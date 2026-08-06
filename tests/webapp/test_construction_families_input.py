"""Algebra-INPUT reachability: the P44/P46/P48 construction families, ℚ as an
input field, and quiver-with-potential (P44 Jacobian) wired into the no-code
webapp/GUI surface.

UNMARKED (no oracle-class markers) per the Plan-32 extras-gated ruling for
webapp/GUI wiring tests -- they need the webapp import surface and the Pyodide
twin file, not the algebra oracles.

Every new construction family is a webapp ``algebra`` block of ``kind == "family"``
(the P50 integration sweep's scope note: these build an algebra, they are not
compute kinds). The Pyodide twin (``docs/gui/runner.py``) does not build ``family``
blocks (pinned refusal, ``tests/gui/test_runner_build.test_family_kind_points_at_plan_09``);
the cross-tier byte agreement is therefore shown by RECONSTRUCTING the built
algebra as an equivalent ``kind == "quiver"`` block and running the SAME cheap kind
through the twin -- the family produces an algebra whose compute is byte-identical
across the two runners.

ℚ and the potential both live on the ``quiver`` block, which BOTH runners build, so
they get direct ``run_build``/``compute_one`` twin agreement (the ``*_runner_twin``
pattern, ``json.dumps(sort_keys=True)``).
"""
import importlib.util
import json
import pathlib
import tempfile

import pydantic
import pytest

from quiverlab.hpc.spec import SpecError, build_algebra, parse_request
from webapp.server.cache import canonical_key
from webapp.server.catalog import build_catalog
from webapp.server.runner import RunError, run_spec
from webapp.server.schema import ComputeRequest

ROOT = pathlib.Path(__file__).resolve().parents[2]
_V = "0.1.0.dev0"          # frozen library version for the cache-key pins (see below)

# --------------------------------------------------------------------------- #
# Runners
# --------------------------------------------------------------------------- #
_PYO = None


def _pyodide_runner():
    """Import ``docs/gui/runner.py`` (the Pyodide twin) as a fresh module -- the exact
    file the browser worker gets, loaded as ``tests/gui/conftest`` and the P50 sweep do."""
    global _PYO
    if _PYO is None:
        path = ROOT / "docs" / "gui" / "runner.py"
        s = importlib.util.spec_from_file_location("_cfam_gui_runner_twin", path)
        m = importlib.util.module_from_spec(s)
        s.loader.exec_module(m)
        _PYO = m
    return _PYO


def _server(body, kind):
    with tempfile.TemporaryDirectory() as d:
        res = run_spec(ComputeRequest.model_validate(body), d)
    assert "error" not in res, res.get("error")
    return res["results"][kind]


def _pyo(body, item):
    gui = _pyodide_runner()
    assert json.loads(gui.run_build(json.dumps(body)))["ok"], "Pyodide run_build failed"
    out = json.loads(gui.compute_one(item))
    assert out["ok"], out
    return out["block"]


def _canon(block):
    """The block's MATH content, provenance dropped, canonical JSON. The server (raw
    Python, int dict-keys) and the Pyodide twin (already crossed JSON, so string keys)
    compare EQUAL when the mathematics agrees -- the ``*_runner_twin`` idiom."""
    math = {k: v for k, v in block.items() if k not in ("references", "citations")}
    return json.dumps(json.loads(json.dumps(math, default=str)), sort_keys=True)


# --------------------------------------------------------------------------- #
# The five construction families (flattened, catalog-expressible params).
# --------------------------------------------------------------------------- #
_GF5 = {"kind": "GF", "p": 5, "n": 1}
# name -> (params, produced dim over GF(5))
_FAMILIES = {
    # line-graph Brauer tree = kZ_2/J^3, dim 6 (tests/families/test_brauer)
    "BrauerGraphAlgebra": ({"edges": [[0, 1], [1, 2]],
                            "cyclic_order": [[0, [0]], [1, [0, 1]], [2, [1]]],
                            "multiplicities": [[0, 1], [1, 1], [2, 1]]}, 6),
    # A2[S_1] one-point extension: 1 + dim S_1 + dim kA2 = 1 + 1 + 3 = 5
    "OnePointExtension": ({"base": "A2", "module_kind": "simple", "module_vertex": 1}, 5),
    # eAe of kA3 at {1, 2}
    "CornerAlgebra": ({"base": "A3", "vertices": [1, 2]}, 3),
    # (kA3)^op
    "OppositeAlgebra": ({"base": "A3"}, 6),
    # gentle Jacobian of the disc fan -> gentle A3
    "MarkedSurface": ({"preset": "disc_fan_A3"}, 6),
}


def _family_body(name, params, compute):
    return {"schema": 2,
            "algebra": {"kind": "family", "family": name, "params": params, "field": _GF5},
            "compute": compute, "artifacts": {"pdf": False, "tikz": False}}


def _quiver_body(A, compute):
    return {"schema": 2,
            "algebra": {"kind": "quiver", "vertices": list(A.quiver.vertices),
                        "arrows": {a: list(A.quiver.arrows[a]) for a in A.quiver.arrows},
                        "relations": [str(r) for r in (A.relations or [])], "field": _GF5},
            "compute": compute, "artifacts": {"pdf": False, "tikz": False}}


# --------------------------------------------------------------------------- #
# (1) every new family builds + computes one cheap kind end-to-end (server tier)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name", sorted(_FAMILIES))
def test_family_builds_and_computes_end_to_end(name):
    params, dim = _FAMILIES[name]
    block = _server(_family_body(name, params, ["recognizers"]), "recognizers")
    assert isinstance(block, dict) and "error" not in block, block
    # the produced algebra has the expected dimension over GF(5) (also gated by the
    # catalog prefill test, restated here per-family as an explicit contract)
    A = build_algebra({"kind": "family", "family": name, "params": params, "field": _GF5})
    assert A.dim == dim, (name, A.dim)


# --------------------------------------------------------------------------- #
# (2) cross-tier byte agreement: the family (server) == the reconstructed quiver
#     run through the Pyodide twin (json sort_keys equal -- the *_runner_twin idiom)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name", sorted(_FAMILIES))
def test_family_math_byte_agrees_with_pyodide_twin(name):
    params, _ = _FAMILIES[name]
    server_block = _server(_family_body(name, params, ["recognizers"]), "recognizers")
    A = build_algebra({"kind": "family", "family": name, "params": params, "field": _GF5})
    qb = _quiver_body(A, ["recognizers"])
    # the reconstructed quiver is a byte-identical presentation of the family algebra
    # (same cartan/basis), so the twin -- which cannot build families -- computes the
    # SAME mathematics as the family did on the server.
    assert _canon(server_block) == _canon(_pyo(qb, "recognizers")), name
    # ... and the server agrees on the quiver form too (same algebra, same runner)
    assert _canon(server_block) == _canon(_server(qb, "recognizers")), name


# --------------------------------------------------------------------------- #
# (3) catalog: every synthetic family is listed with a 4-language summary,
#     per-param help in every language, and a verified prefill
# --------------------------------------------------------------------------- #
def test_synthetic_families_in_catalog_with_full_i18n():
    from webapp.server.i18n import LANGS
    fams = {f["name"]: f for f in build_catalog()["families"]}
    for name in _FAMILIES:
        assert name in fams, f"{name} missing from the catalog"
        fam = fams[name]
        assert set(fam.get("summary", {})) == set(LANGS), name
        assert set(fam["fields"]) >= {"CC", "GF", "QQ"}, name
        for p in fam["params"]:
            assert set(p.get("help", {})) == set(LANGS), (name, p["name"])


# --------------------------------------------------------------------------- #
# (4) ℚ as an input field
# --------------------------------------------------------------------------- #
def test_qq_field_builds_and_matches_cc_on_ka2():
    # kA2 HH^0..2 over QQ equals the value over CC (they coincide there: [1, 0, 0]).
    ka2 = lambda fk: {"schema": 1,   # noqa: E731
                      "algebra": {"kind": "quiver", "vertices": [1, 2],
                                  "arrows": {"a": [1, 2]}, "relations": [], "field": fk},
                      "compute": ["hh_cohomology:0..2"], "artifacts": {"pdf": False}}
    qq = _server(ka2({"kind": "QQ"}), "hh_cohomology")
    cc = _server(ka2({"kind": "CC"}), "hh_cohomology")
    assert qq["dims"] == cc["dims"] == [1, 0, 0]


def test_qq_field_differs_from_gf_when_it_should():
    # k[x]/(x^2): HH_* is characteristic-sensitive. Over QQ (char 0) HH_n = 1 for n >= 1;
    # over GF(2) it is 2 -- so they DIFFER (cheap witness that QQ is a real field choice,
    # not silently GF/CC).
    dual = lambda fk: {"schema": 1,   # noqa: E731
                       "algebra": {"kind": "quiver", "vertices": [1],
                                   "arrows": {"x": [1, 1]}, "relations": ["x*x"],
                                   "field": fk},
                       "compute": ["hh_homology:0..4"], "artifacts": {"pdf": False}}
    qq = _server(dual({"kind": "QQ"}), "hh_homology")["dims"]
    gf = _server(dual({"kind": "GF", "p": 2, "n": 1}), "hh_homology")["dims"]
    assert qq == [2, 1, 1, 1, 1]
    assert gf == [2, 2, 2, 2, 2]
    assert qq != gf


def test_qq_recognizers_twin():
    body = {"schema": 2,
            "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                        "relations": [], "field": {"kind": "QQ"}},
            "compute": ["recognizers"], "artifacts": {"pdf": False, "tikz": False}}
    assert _canon(_server(body, "recognizers")) == _canon(_pyo(body, "recognizers"))


# --------------------------------------------------------------------------- #
# (5) quiver-with-potential -> P44 Jacobian
# --------------------------------------------------------------------------- #
_TRIANGLE = {"kind": "quiver", "vertices": [1, 2, 3],
             "arrows": {"a": [1, 2], "b": [2, 3], "c": [3, 1]}}


def _pot_body(compute, potential="a*b*c", relations=None, field=None):
    alg = dict(_TRIANGLE, potential=potential, field=(field or _GF5))
    if relations is not None:
        alg["relations"] = relations
    return {"schema": 2, "algebra": alg, "compute": compute,
            "artifacts": {"pdf": False, "tikz": False}}


def test_potential_triangle_is_cyclic_nakayama_kz3_j2():
    # P44 oracle: Jac(3-cycle, abc) = kZ_3/J^2 (cyclic Nakayama, dim 6, identical
    # Cartan). Reuses the numbers from tests/families/test_jacobian.
    from quiverlab import GF, NakayamaAlgebra
    J = build_algebra(dict(_pot_body(["cartan"])["algebra"]))
    N = NakayamaAlgebra(n=3, l=2, cyclic=True, field=GF(5))
    assert J.dim == N.dim == 6
    assert [list(r) for r in J.cartan_matrix()] == [list(r) for r in N.cartan_matrix()]
    assert sorted(str(r) for r in J.relations) == ["a*b", "b*c", "c*a"]


def test_potential_twin_byte_agrees():
    body = _pot_body(["recognizers"])
    assert _canon(_server(body, "recognizers")) == _canon(_pyo(body, "recognizers"))


def test_potential_and_relations_conflict_is_loud_4xx():
    # At the HTTP boundary (pydantic) it is a ValidationError -> typed 4xx; in the
    # shared spec core it is a SpecError. Both are loud, neither a 500.
    bad = _pot_body(["cartan"], potential="a*b*c", relations=["a*b"])
    with pytest.raises(pydantic.ValidationError):
        ComputeRequest.model_validate(bad)
    with pytest.raises(SpecError):
        parse_request(bad)


def test_non_cyclic_potential_refuses_cleanly():
    # A non-cyclic potential term is the library's loud QuiverlabError, surfaced as a
    # clean typed RunError (never a bare 500). "a" alone is not a cycle (source 1 != target 2).
    body = _pot_body(["cartan"], potential="a")
    with pytest.raises(RunError):
        _server(body, "cartan")


# --------------------------------------------------------------------------- #
# (6) canonical-key byte-stability: an ABSENT potential serializes away, so every
#     pre-change request hashes IDENTICALLY (keys computed on clean HEAD before the
#     schema gained the ``potential`` field).
# --------------------------------------------------------------------------- #
_PRECHANGE_KEYS = {
    # triangle quiver, GF(5), recognizers -- same quiver shape as the potential test,
    # but WITHOUT a potential: must hash to its pre-change value.
    "triangle_gf5_recognizers": (
        {"schema": 2, "algebra": {"kind": "quiver", "vertices": [1, 2, 3],
                                  "arrows": {"a": [1, 2], "b": [2, 3], "c": [3, 1]},
                                  "relations": [], "field": {"kind": "GF", "p": 5, "n": 1}},
         "compute": ["recognizers"], "artifacts": {"pdf": False, "tikz": False}},
        "4c8249d821df4a104f6436a463213e13e8662f28907b05d9ab76c16bd88b1ab3"),
    "kA2_cc_hh": (
        {"schema": 1, "algebra": {"kind": "quiver", "vertices": [1, 2],
                                  "arrows": {"a": [1, 2]}, "relations": [],
                                  "field": {"kind": "CC"}},
         "compute": ["hh_cohomology:0..2"], "artifacts": {"pdf": False, "tikz": False}},
        "c0724fde898761428d104c848c4d3a4a57c330897e8a9b4fe181aadd47f627cb"),
    "loop_gf2_cartan": (
        {"schema": 1, "algebra": {"kind": "quiver", "vertices": [1],
                                  "arrows": {"x": [1, 1]}, "relations": ["x*x*x"],
                                  "field": {"kind": "GF", "p": 2, "n": 1}},
         "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False}},
        "20baab2b9008b580f20c6eeafeffe7e0aba004f7e72ed12e90677734460a9272"),
}


@pytest.mark.parametrize("name", sorted(_PRECHANGE_KEYS))
def test_absent_potential_keeps_canonical_key_byte_stable(name):
    body, pre = _PRECHANGE_KEYS[name]
    key = canonical_key(ComputeRequest.model_validate(body).model_dump(by_alias=True), _V)
    assert key == pre, f"cache key for {name!r} drifted -- absent potential must serialize away"
