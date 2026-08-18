"""The `cluster_category` compute kind through the server/HPC tier (Plan 79 / R31, G1).

Unmarked (extras-gated dir, the Plan-32 ruling): tests/webapp may never carry `oracle_*`
markers, since that would make the audited class counts depend on installed extras.
"""
import pathlib
import tempfile

import pytest

from quiverlab.hpc.spec import parse_compute_item, parse_request
from quiverlab.hpc.spec import run as spec_run
from webapp.server.cache import canonical_key
from webapp.server.schema import ComputeRequest as WebRequest
from webapp.server.schema import parse_compute_item as web_parse

_A3 = {"kind": "quiver", "vertices": [1, 2, 3], "arrows": {"a": [1, 2], "b": [2, 3]},
       "relations": [], "field": {"kind": "QQ"}}


def _line(n):
    return {"kind": "quiver", "vertices": list(range(1, n + 1)),
            "arrows": {f"a{i}": [i, i + 1] for i in range(1, n)},
            "relations": [], "field": {"kind": "QQ"}}


def _dynkin_D(n):
    """A genuine D_n: the fork hangs off vertex n-2, not off the end of the chain."""
    arrows = {f"a{i}": [i, i + 1] for i in range(1, n - 1)}
    arrows["b"] = [n - 2, n]
    return {"kind": "quiver", "vertices": list(range(1, n + 1)), "arrows": arrows,
            "relations": [], "field": {"kind": "QQ"}}


def _run(algebra=None, compute="cluster_category"):
    body = {"schema": 1, "algebra": algebra or _A3, "compute": [compute],
            "artifacts": {"pdf": False, "tikz": False}}
    with tempfile.TemporaryDirectory() as td:
        return spec_run(parse_request(body), pathlib.Path(td))


def test_the_kA3_block_reproduces_every_pin():
    b = _run()["results"]["cluster_category"]
    assert b["kind"] == "cluster_category"
    assert b["n"] == 3 and b["dynkin_type"] == "A3" and b["hereditary"] is True
    assert b["num_indec"] == 9
    assert b["num_cluster_tilting"]["count"] == 14
    assert b["num_cluster_tilting"]["certified"] is True
    assert b["cluster_tilted"]["dim"] == 6
    assert b["two_cy"]["ar_formula_holds"] is True
    assert "amiot_cluster_category" in b["references"]


def test_every_citation_pair_is_fully_formatted():
    b = _run()["results"]["cluster_category"]
    assert len(b["citations"]) == len(b["references"])
    for key, formatted in b["citations"]:
        assert "()" not in formatted, f"{key} lost its year: {formatted}"


def test_the_budget_grammar_is_parsed_IDENTICALLY_by_all_three_sites():
    # The special form bypasses the 'name:0..N' degree grammar; if the three sites
    # disagreed about what a request means the tiers would silently diverge.
    for spec, hi in (("cluster_category", None), ("cluster_category:64", 64)):
        a = parse_compute_item(spec)
        b = web_parse(spec)
        assert a.kind == b.kind == "cluster_category"
        assert a.hi == b.hi == hi
        assert a.lo is None and b.lo is None


def test_a_bad_budget_is_refused_by_both_server_grammars():
    from quiverlab.hpc.spec import SpecError
    from webapp.server.schema import SchemaError
    with pytest.raises(SpecError, match="positive integer"):
        parse_compute_item("cluster_category:abc")
    with pytest.raises(SchemaError, match="positive integer"):
        web_parse("cluster_category:abc")


def test_a_non_hereditary_algebra_yields_a_typed_error_block_not_a_500():
    alg = dict(_A3, relations=["a*b"])
    b = _run(algebra=alg)["results"]["cluster_category"]
    assert "error" in b and "hereditary" in b["error"]
    assert b["kind"] == "cluster_category"


def test_a_cyclic_quiver_is_refused_as_a_TYPED_error_before_the_cluster_layer():
    # Through the SERVICE the refusal lands one layer earlier than the cluster gate: a
    # cyclic quiver's path algebra is infinite-dimensional, so kQ cannot be built at all
    # and the runner raises a typed ComputeError -- never a 500. The cluster layer's own
    # "oriented cycle" gate is what a LIBRARY caller passing a Quiver object hits
    # (tests/cluster/test_category_indec.py::test_G2_*), so both layers are covered and
    # neither is assumed to cover the other.
    from quiverlab.hpc.spec import ComputeError
    alg = {"kind": "quiver", "vertices": [1, 2],
           "arrows": {"a": [1, 2], "b": [2, 1]}, "relations": [],
           "field": {"kind": "QQ"}}
    with pytest.raises(ComputeError) as exc:
        _run(algebra=alg)
    assert "infinite-dimensional" in str(exc.value)


def test_a_representation_infinite_input_is_a_POPULATED_refusal_not_an_error_block():
    # The Kronecker: the CATEGORY is a legitimate object, only the enumeration is
    # infinite -- so this must be a populated block with an honest note, NOT an error
    # block that hides the rest of the payload.
    alg = {"kind": "quiver", "vertices": [1, 2],
           "arrows": {"a": [1, 2], "b": [1, 2]}, "relations": [],
           "field": {"kind": "QQ"}}
    b = _run(algebra=alg)["results"]["cluster_category"]
    assert "error" not in b
    assert b["num_indec"] is None
    assert "infinitely many" in b["note"]
    assert b["num_cluster_tilting"]["certified"] is False


def test_reproduce_snippet_names_the_public_method():
    assert "A.cluster_category()" in _run()["reproduce"]


def test_the_estimator_sizes_on_the_CLUSTER_NUMBER_not_the_algebra_dimension():
    # The plan's MAJOR-2 item. E6 has dim kQ = 36 -- comfortably instant under the
    # dim**3 bar model -- but 833 cluster-tilting objects and a minutes-long BFS. The
    # sizing must therefore route the big types OFF the instant tier while leaving the
    # small ones (and every non-cluster request) exactly where they were.
    from webapp.server.config import Config
    from webapp.server.estimator import classify
    cfg = Config.from_env()
    cases = [(_line(3), 6, "instant"), (_line(4), 10, "instant"),
             (_line(5), 15, "queued"), (_line(7), 28, "queued"),
             (_dynkin_D(5), 20, "queued")]
    for alg, dim, expected in cases:
        req = WebRequest.model_validate(
            {"schema": 1, "algebra": alg, "compute": ["cluster_category"]})
        assert classify(dim, req, cfg)["tier"] == expected, (alg["vertices"], expected)
        # ... and the SAME algebra asking for something else is untouched
        plain = WebRequest.model_validate(
            {"schema": 1, "algebra": alg, "compute": ["dimension"]})
        assert classify(dim, plain, cfg)["tier"] == "instant"


def test_adding_cluster_category_leaves_an_existing_request_key_untouched():
    plain = {"schema": 1, "algebra": _A3, "compute": ["dimension"]}
    k1 = canonical_key(WebRequest.model_validate(plain).model_dump(by_alias=True), "v")
    k2 = canonical_key(WebRequest.model_validate(plain).model_dump(by_alias=True), "v")
    assert k1 == k2
    withcc = dict(plain, compute=["dimension", "cluster_category"])
    assert canonical_key(
        WebRequest.model_validate(withcc).model_dump(by_alias=True), "v") != k1
