"""docs/gui/runner.py serves the cyclic_homology kind SHAPE-IDENTICALLY to the server
tier, INCLUDING the Plan 35 wave-3b explicit HC representatives (basis_classes /
chain_basis / differentials / column_structure). The Pyodide twin and the wheel's spec
core cannot import each other, so the contract is that the block one emits is
key-for-key equal to the other's -- and both carry the reps + the column structure.
"""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

_BODY = {"schema": 2,
         "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                     "relations": ["x*x"], "field": {"kind": "GF", "p": 7, "n": 1}},
         "compute": ["cyclic_homology:0..3"],
         "artifacts": {"pdf": False, "tikz": False}}

_COMPARE = ("kind", "top", "dims", "engine", "references", "basis_classes",
            "chain_basis", "differentials", "column_structure")


def _load_gui_runner():
    spec = importlib.util.spec_from_file_location("gui_runner_cyc", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _gui_block():
    gui = _load_gui_runner()
    assert json.loads(gui.run_build(json.dumps(_BODY)))["ok"]
    out = json.loads(gui.compute_one("cyclic_homology:0..3"))
    assert out["ok"], out
    return out["block"]


def _server_block(tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(_BODY), tmp_path)
    return ref["results"]["cyclic_homology"]


def test_cyclic_blocks_match_spec_runner(tmp_path):
    gui = _gui_block()
    srv = _server_block(tmp_path)
    for key in _COMPARE:
        assert gui.get(key) == srv.get(key), "cyclic block mismatch on %r" % key
    # the reps really are present (not both-None-therefore-equal)
    assert gui.get("basis_classes") and gui.get("column_structure")
