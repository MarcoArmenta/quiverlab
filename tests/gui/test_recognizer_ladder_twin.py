"""The Pyodide GUI runner (docs/gui/runner.py) serves the `recognizer_ladder` kind (Plan 61)
SHAPE-IDENTICALLY to the server core (quiverlab.hpc.spec): the SAME library block builder
(modules.recognizers_ladder.recognizer_ladder_block) + `references`->citations, so the
quasi-tilted/shod/weakly-shod/laura/ada ladder + the ada/HH^1 simple-connectedness block
cannot drift between the browser GUI and the server report."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

_BODY = {"schema": 1,
         "algebra": {"kind": "quiver", "vertices": [1, 2, 3],
                     "arrows": {"a": [1, 2], "b": [2, 3]},
                     "relations": [], "field": {"kind": "QQ"}},
         "compute": ["recognizer_ladder:256"], "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner_rl_twin", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("recognizer_ladder:256"))
    assert out["ok"], out
    return out["block"]


def test_gui_recognizer_ladder_block_shape():
    b = _gui_block(_BODY)
    assert b["complete"] and b["n"] == 3 and b["universe_size"] == 6 and b["gldim"] == 1
    assert all(b["ladder"][k]["verdict"] for k in
               ("quasi_tilted", "shod", "weakly_shod", "laura", "ada"))
    assert b["complement"] == []
    assert "hrs_quasitilted" in b["references"]


def test_gui_matches_server(tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    server = run_spec(ComputeRequest.model_validate(_BODY),
                      tmp_path)["results"]["recognizer_ladder"]
    gui = _gui_block(_BODY)
    assert json.dumps(server, sort_keys=True) == json.dumps(gui, sort_keys=True)
