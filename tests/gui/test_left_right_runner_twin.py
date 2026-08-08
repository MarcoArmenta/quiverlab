"""The Pyodide GUI runner (docs/gui/runner.py) serves the `left_right_parts` kind (Plan 55)
SHAPE-IDENTICALLY to the server core (quiverlab.hpc.spec): the SAME library block builder
(modules.left_right.left_right_parts_block) + `references`->citations, so the module-category
atlas (both parts + complement + Ext-injectives + supports) cannot drift between the browser
GUI and the server report."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

_BODY = {"schema": 1,
         "algebra": {"kind": "quiver", "vertices": [1, 2, 3],
                     "arrows": {"a": [1, 2], "b": [2, 3]},
                     "relations": [], "field": {"kind": "QQ"}},
         "compute": ["left_right_parts:256"], "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner_lr_twin", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("left_right_parts:256"))
    assert out["ok"], out
    return out["block"]


def test_gui_left_right_parts_block_shape():
    b = _gui_block(_BODY)
    assert b["complete"] and b["n"] == 3 and b["universe_size"] == 6
    assert len(b["left"]) == len(b["right"]) == 6
    assert b["complement"] == []
    assert set(b["left_support"]["vertices"]) == {1, 2, 3}
    assert "act_left_right" in b["references"]


def test_gui_matches_server(tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    server = run_spec(ComputeRequest.model_validate(_BODY),
                      tmp_path)["results"]["left_right_parts"]
    gui = _gui_block(_BODY)
    assert json.dumps(server, sort_keys=True) == json.dumps(gui, sort_keys=True)
