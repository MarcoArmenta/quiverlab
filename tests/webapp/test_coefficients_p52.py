"""Plan 52 coefficient block end-to-end (both runners). A schema-3 request carries
a Hochschild coefficient bimodule; the hh block gains a ``coefficients`` provenance
string; a coefficients-LESS request keeps its byte-identical canonical key
(conditional versioning); and a big explicit coefficient sizes the job off instant.

No oracle-class marker (tests/webapp/ collects only with the [web] extra; contract
tests for the two runners' block shapes, per the Plan-32 extras-gated ruling)."""
import json
import tempfile
import pathlib
import importlib.util

import pytest
import quiverlab as ql

from quiverlab.hpc import spec

_A2 = {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
       "relations": [], "field": {"kind": "CC"}}


def _server(req):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]


def _pyodide(req):
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_p52", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    out = {}
    for item in req["compute"]:
        r = json.loads(mod.compute_one(item))
        assert r["ok"], r
        out[item.split(":")[0]] = r["block"]
    return out


def _req_dual():
    return {"schema": 3, "algebra": _A2, "compute": ["hh_cohomology:0..4"],
            "coefficients": {"builtin": {"kind": "dual"}}}


def _expected_dims_dual():
    A = ql.Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=ql.CC)
    return A.hochschild_cohomology(4, engine="bar", coefficients=ql.Bimodule.dual(A)).dims


def test_coefficient_block_shape():
    out = _server(_req_dual())
    blk = out["hh_cohomology"]
    assert blk["coefficients"] == "D(A)"
    assert blk["dims"] == _expected_dims_dual()


def test_twin_parity():
    sv = _server(_req_dual())["hh_cohomology"]
    py = _pyodide(_req_dual())["hh_cohomology"]
    # the coefficient contract fields must agree across the two runners
    assert sv["coefficients"] == py["coefficients"] == "D(A)"
    assert sv["dims"] == py["dims"]
    assert sv["kind"] == py["kind"]


def test_absent_coefficients_keeps_canonical_key():
    from webapp.server.schema import ComputeRequest
    from webapp.server.cache import canonical_key
    plain = {"schema": 2, "algebra": _A2, "compute": ["hh_cohomology:0..4"]}
    r = ComputeRequest(**plain)
    d = r.model_dump(by_alias=True)
    # conditional versioning: the absent coefficients block is POPPED, so the
    # canonical form (and hence the cache key) is byte-identical to pre-Plan-52.
    assert "coefficients" not in d
    k = canonical_key(d, "x")
    # a hand-built dict WITHOUT the coefficients key must produce the SAME key
    manual = {"schema": 2, "algebra": r.model_dump(by_alias=True)["algebra"],
              "compute": ["hh_cohomology:0..4"],
              "artifacts": d["artifacts"]}
    assert canonical_key(manual, "x") == k


def test_absent_coefficients_v3_still_pops():
    # a v3 request with NO coefficients block still drops the key (pop-list), so a
    # gratuitous version bump alone does not re-key.
    from webapp.server.schema import ComputeRequest
    plain = {"schema": 3, "algebra": _A2, "compute": ["hh_cohomology:0..4"]}
    d = ComputeRequest(**plain).model_dump(by_alias=True)
    assert "coefficients" not in d


def _req_big_explicit_coeff(dim_M):
    return {"schema": 3, "algebra": {"kind": "quiver", "vertices": [1],
                                     "arrows": {"x": [1, 1]}, "relations": ["x*x"],
                                     "field": {"kind": "GF", "p": 5}},
            "compute": ["hh_cohomology:0..2"],
            "coefficients": {"dim": dim_M,
                             "left_maps": {"x": [[0]]}, "right_maps": {"x": [[0]]}}}


def test_large_explicit_coefficient_does_not_route_instant():
    # C4: an explicit coefficient of large dim over a SMALL algebra must size the
    # job (bar cost is quadratic in dim_M) -> NOT the instant tier.
    from webapp.server.estimator import sizing_dim
    from webapp.server.schema import ComputeRequest
    req = ComputeRequest(**_req_big_explicit_coeff(dim_M=60))
    assert sizing_dim(3, req) >= 60             # coefficient dominates the small algebra


def test_coefficients_require_schema_3():
    from webapp.server.schema import ComputeRequest
    from pydantic import ValidationError
    bad = {"schema": 2, "algebra": _A2, "compute": ["hh_cohomology:0..4"],
           "coefficients": {"builtin": {"kind": "dual"}}}
    with pytest.raises(ValidationError):
        ComputeRequest(**bad)


def test_coefficients_only_on_hh_kinds():
    from webapp.server.schema import ComputeRequest
    from pydantic import ValidationError
    bad = {"schema": 3, "algebra": _A2, "compute": ["center"],
           "coefficients": {"builtin": {"kind": "dual"}}}
    with pytest.raises(ValidationError):
        ComputeRequest(**bad)
