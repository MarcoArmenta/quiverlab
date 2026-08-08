"""radical_filtration + ar_invariants algebra kinds (Plan 57): served by hpc.spec,
mirrored byte-identically by the Pyodide twin (docs/gui/runner.py), with the honest
off-scope banner. Unmarked: webapp-contract + cross-runner infrastructure, not a
mathematical oracle class (the Plan-32 extras-gated ruling -- the math oracles live
in tests/modules/)."""
import importlib.util
import json
import os
import pathlib
import tempfile

import pytest

from quiverlab.hpc import spec

_ROOT = pathlib.Path(__file__).resolve().parents[2]

_kA3 = {"kind": "quiver", "vertices": [1, 2, 3], "arrows": {"a1": [1, 2], "a2": [2, 3]},
        "relations": [], "field": {"kind": "QQ"}}
_kxx = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
        "relations": ["x*x"], "field": {"kind": "GF", "p": 5, "n": 1}}
# The non-uniform cyclic Nakayama [3,2,2] over GF(2): knittable + NON-self-injective,
# but char 2 <= dim 3 trips the trace-form rad-End char guard -> a loud typed error.
_n322_gf2 = {"kind": "quiver", "vertices": [1, 2, 3],
             "arrows": {"a1": [1, 2], "a2": [2, 3], "a3": [3, 1]},
             "relations": ["a1*a2*a3", "a2*a3", "a3*a1"],
             "field": {"kind": "GF", "p": 2, "n": 1}}


def _server_block(req, kind):
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"][kind]


def _twin_block(req, compute_key):
    path = os.path.join(_ROOT, "docs", "gui", "runner.py")
    s = importlib.util.spec_from_file_location("gui_twin_p57_test", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    r = json.loads(mod.compute_one(compute_key))
    assert r["ok"], r
    return r["block"]


def test_radical_filtration_block_shape():
    b = _server_block({"schema": 1, "algebra": _kA3, "compute": ["radical_filtration"]},
                      "radical_filtration")
    assert b["kind"] == "radical_filtration"
    assert b["complete"] is True and b["status"] == "complete"
    assert b["nilpotency_index"] == 3
    assert b["rad_infinity_zero"] is True
    assert b["generalized_standard"] is True
    assert b["layer_profile"] == [9, 3]
    # citation pairs are (bibtex_key, formatted) -- the ar_quiver/ARS1995 precedent.
    assert "ChaioLiu2013" in [k for k, _ in b["citations"]]


def test_ar_invariants_block_shape():
    b = _server_block({"schema": 1, "algebra": _kA3, "compute": ["ar_invariants"]},
                      "ar_invariants")
    assert b["kind"] == "ar_invariants"
    assert b["representation_directed"] is True
    assert "regular" not in b["partition_counts"]
    assert isinstance(b["degrees"], dict) and len(b["degrees"]) == 6
    assert "Liu1992degrees" in [k for k, _ in b["citations"]]


def test_self_injective_off_scope_banner():
    # cyclic Nakayama kZ_3/J^2 is self-injective (a single loop x^2 here): status
    # 'unsupported', no fake index.
    for kind in ("radical_filtration", "ar_invariants"):
        b = _server_block({"schema": 1, "algebra": _kxx, "compute": [kind]}, kind)
        assert b["complete"] is False and b["status"] == "unsupported"
        assert b["nilpotency_index"] is None


@pytest.mark.parametrize("req,key,kind", [
    ({"schema": 1, "algebra": _kA3, "compute": ["radical_filtration"]},
     "radical_filtration", "radical_filtration"),
    ({"schema": 1, "algebra": _kA3, "compute": ["radical_filtration:64"]},
     "radical_filtration:64", "radical_filtration"),
    ({"schema": 1, "algebra": _kA3, "compute": ["ar_invariants"]},
     "ar_invariants", "ar_invariants"),
    ({"schema": 1, "algebra": _kxx, "compute": ["radical_filtration"]},
     "radical_filtration", "radical_filtration"),
])
def test_twin_parity(req, key, kind):
    sb = _server_block(req, kind)
    tb = _twin_block(req, key)
    assert sb.get("kind") == kind
    assert json.dumps(sb, sort_keys=True, default=str) == \
           json.dumps(tb, sort_keys=True, default=str)


@pytest.mark.parametrize("kind", ["radical_filtration", "ar_invariants"])
def test_char_le_dim_yields_a_clean_typed_error_block(kind):
    # char 2 <= dim 3 on [3,2,2]/GF(2): the trace-form rad-End refusal must surface as
    # a clean typed {"kind":..., "error":...} block through the dispatch -- never a 500
    # or an uncaught raise.
    b = _server_block({"schema": 1, "algebra": _n322_gf2, "compute": [kind]}, kind)
    assert b["kind"] == kind
    assert "error" in b and "char 2" in b["error"]
    # a char-scope refusal never emits a fake verdict
    assert "nilpotency_index" not in b or b.get("nilpotency_index") is None


def test_budget_is_not_a_degree_range():
    # radical_filtration / ar_invariants carry a module budget, not a degree range;
    # the 'name:0..N' form is rejected at parse time (like ar_quiver / tau_tilting).
    with pytest.raises(spec.SpecError):
        spec.parse_compute_item("radical_filtration:0..4")
    assert spec.parse_compute_item("radical_filtration:512").hi == 512
    assert spec.parse_compute_item("ar_invariants").hi is None
    # and radical_filtration must NOT collide with radical_filtration_ss (a range kind)
    assert spec.parse_compute_item("radical_filtration_ss:0..4").kind == \
        "radical_filtration_ss"
