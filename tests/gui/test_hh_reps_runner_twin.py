"""docs/gui/runner.py serves the PLAIN hh_cohomology / hh_homology kinds
SHAPE-IDENTICALLY to the server tier, INCLUDING the Plan 35 wave-3d explicit
representatives (basis_classes / chain_basis / differentials / inner_dims). The Pyodide
twin and the wheel's spec core cannot import each other, so the contract is that the
block one emits is key-for-key equal to the other's -- and both carry the reps.
"""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

# k<x>/(x^3) over GF(2): the F_p fast dims path -> the GF(p) bar reps route.
_BODY = {"schema": 2,
         "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                     "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
         "compute": ["hh_cohomology:0..3", "hh_homology:0..3"],
         "artifacts": {"pdf": False, "tikz": False}}
# gui block kinds are HHTable's kind strings ("HH^" / "HH_"); server results are keyed
# by the compute name.
_MAP = {"hh_cohomology": "HH^", "hh_homology": "HH_"}
_COMPARE = ("dims", "engine", "basis_classes", "chain_basis", "differentials",
            "inner_dims")


def _load_gui_runner():
    spec = importlib.util.spec_from_file_location("gui_runner", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _gui_blocks(body):
    gui = _load_gui_runner()
    assert json.loads(gui.run_build(json.dumps(body)))["ok"]
    blocks = {}
    for spec in body["compute"]:
        out = json.loads(gui.compute_one(spec))
        assert out["ok"], out
        blocks[out["block"]["kind"]] = out["block"]
    return blocks


def _server_blocks(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"]


def test_hh_blocks_match_spec_runner(tmp_path):
    gui = _gui_blocks(_BODY)
    ref = _server_blocks(_BODY, tmp_path)
    for name, gkind in _MAP.items():
        for key in _COMPARE:
            assert gui[gkind].get(key) == ref[name].get(key), (name, key)


def test_both_runners_carry_hh_reps(tmp_path):
    """Both runners emit basis_classes / chain_basis / differentials / inner_dims (the
    key-for-key equality passes even if BOTH dropped them, so pin PRESENCE)."""
    gui = _gui_blocks(_BODY)
    ref = _server_blocks(_BODY, tmp_path)
    for name, gkind in _MAP.items():
        for field in ("basis_classes", "chain_basis", "differentials", "inner_dims"):
            assert field in gui[gkind], (name, field, "gui")
            assert field in ref[name], (name, field, "server")
        assert sorted(gui[gkind]["basis_classes"]) == ["0", "1", "2", "3"], name
