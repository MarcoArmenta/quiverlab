"""The `koszul` compute kind end to end through the server/HPC tier (Plan 77 / R36).

Unmarked (extras-gated dir, the Plan-32 ruling): tests/webapp may never carry
`oracle_*` markers, since that would make the audited class counts depend on which
extras are installed.
"""
import json
import pathlib
import tempfile

import pytest

from quiverlab.hpc.spec import parse_compute_item, parse_request
from quiverlab.hpc.spec import run as spec_run
from webapp.server.cache import canonical_key
from webapp.server.schema import ComputeRequest as WebRequest


def _run(compute, relations=("x*x*x",), vertices=(1,), arrows=None, field=None):
    body = {
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": list(vertices),
                    "arrows": arrows or {"x": [1, 1]},
                    "relations": list(relations),
                    "field": field or {"kind": "QQ"}},
        "compute": list(compute),
    }
    with tempfile.TemporaryDirectory() as td:
        out = spec_run(parse_request(body), pathlib.Path(td))
    return body, out


def test_koszul_block_shape():
    _, out = _run(["koszul:0..8"])
    b = out["results"]["koszul"]
    assert b["kind"] == "koszul"
    assert b["n_homogeneous"] == 3
    assert b["k2"]["generator_degrees"] == [1, 2]
    assert b["generation_degrees"]["1"][:4] == [0, 1, 3, 4]   # Berger delta, N=3
    assert "berger_nonquadratic" in b["references"]
    assert b["citations"] and all(len(pair) == 2 for pair in b["citations"])


def test_every_citation_pair_is_fully_formatted():
    # The Plan-77 bib entries must render with their YEAR: an entry whose last field
    # lacks a trailing comma has that field dropped by the parser, which prints
    # "Author (). Title." -- caught live while adding these citations.
    _, out = _run(["koszul"])
    for key, formatted in out["results"]["koszul"]["citations"]:
        assert "()" not in formatted, f"{key} lost its year: {formatted}"
        assert formatted.strip()


def test_default_top_is_eight():
    item = parse_compute_item("koszul")
    assert item.kind == "koszul" and item.hi is None
    _, out = _run(["koszul"])
    assert out["results"]["koszul"]["certified_through_degree"] == 8


def test_range_grammar_matches_the_ext_algebra_sibling():
    item = parse_compute_item("koszul:0..5")
    assert item.kind == "koszul" and item.lo == 0 and item.hi == 5
    _, out = _run(["koszul:0..5"])
    assert out["results"]["koszul"]["certified_through_degree"] == 5


def test_preprojective_A3_is_classified_almost_koszul():
    body = {"schema": 1,
            "algebra": {"kind": "family", "family": "PreprojectiveAlgebra",
                        "params": {"type_or_quiver": "A3"}, "field": {"kind": "QQ"}},
            "compute": ["koszul:0..6"]}
    with tempfile.TemporaryDirectory() as td:
        out = spec_run(parse_request(body), pathlib.Path(td))
    b = out["results"]["koszul"]
    assert b["almost_koszul"]["verdict"] is True
    assert (b["almost_koszul"]["p"], b["almost_koszul"]["q"]) == (2, 2)
    assert b["k2"]["verdict"] is False


def test_adding_koszul_leaves_an_existing_request_key_untouched():
    # `koszul` adds NO request field -- it is a compute-list string -- so a request that
    # does not ask for it keys byte-identically to before Plan 77.
    plain = {"schema": 1,
             "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                         "relations": ["x*x*x"], "field": {"kind": "QQ"}},
             "compute": ["dimension"]}
    with_kz = dict(plain, compute=["dimension", "koszul:0..4"])
    k1 = canonical_key(WebRequest.model_validate(plain).model_dump(by_alias=True), "v")
    k2 = canonical_key(WebRequest.model_validate(with_kz).model_dump(by_alias=True), "v")
    assert k1 != k2                      # the new kind DOES change its own request
    k1b = canonical_key(WebRequest.model_validate(plain).model_dump(by_alias=True), "v")
    assert k1 == k1b                     # ... and leaves the old one byte-stable


def test_presentation_less_refusal_is_typed_and_currently_unreachable():
    # The runner wraps koszul_profile_block in a QuiverlabError catch so a
    # presentation-less algebra becomes an `error` field, never a traceback out of the
    # service. Two honest halves, no skip:
    #  (a) the refusal the catch converts is real and TYPED at the library level, and
    #  (b) today it is UNREACHABLE through the service -- the request schema exposes
    #      only `quiver` and `family` algebras, both of which are presented -- so the
    #      catch is a guard for a future presentation-less family, not dead code
    #      pretending to be covered.
    import typing

    from quiverlab import Algebra
    from quiverlab.errors import QuiverlabError
    from quiverlab.fields import QQ
    from quiverlab.modules.nkoszul import koszul_profile_block
    from webapp.server.schema import AlgebraSpec

    sc = Algebra(QQ, [[[QQ.one()]]], [QQ.one()])
    assert sc.quiver is None
    with pytest.raises(QuiverlabError):
        koszul_profile_block(sc, top=4)

    kinds = set()
    for arg in typing.get_args(AlgebraSpec):
        for member in typing.get_args(arg):
            fields = getattr(member, "model_fields", None)
            if fields and "kind" in fields:
                kinds.update(typing.get_args(fields["kind"].annotation) or ())
    assert kinds == {"quiver", "family"}, kinds


def test_reproduce_snippet_names_the_public_method():
    _, out = _run(["koszul:0..8"])
    assert "A.koszul_profile(8)" in out["reproduce"]
    _, out_default = _run(["koszul"])
    assert "A.koszul_profile(8)" in out_default["reproduce"]   # the default top


def test_tikz_artifact_carries_the_generation_degree_staircase():
    body = {"schema": 1,
            "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                        "relations": ["x*x*x"], "field": {"kind": "QQ"}},
            "compute": ["koszul:0..6"],
            "artifacts": {"tikz": True}}
    with tempfile.TemporaryDirectory() as td:
        out = spec_run(parse_request(body), pathlib.Path(td))
        tex = list(pathlib.Path(td).rglob("*tikz*"))
        src = tex[0].read_text(encoding="utf-8") if tex else ""
    assert out["results"]["koszul"]["n_homogeneous"] == 3
    if src:
        assert "ell = n" in src or r"\ell" in src
