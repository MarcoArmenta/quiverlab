"""The Pyodide GUI runner (docs/gui/runner.py) serves the `strings` kind (Plan 46)
SHAPE-IDENTICALLY to the server core (quiverlab.hpc.spec): the SAME library block
builder (strings.block.strings_block) + `references`->citations, so census / bands /
rep-type / AG cannot drift between the browser GUI and the server report."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

_BODY = {"schema": 1,
         "algebra": {"kind": "quiver", "vertices": [1, 2, 3],
                     "arrows": {"a": [1, 2], "b": [2, 3]},
                     "relations": ["a*b"], "field": {"kind": "GF", "p": 7, "n": 1}},
         "compute": ["strings"], "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner_strings_twin", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("strings"))
    assert out["ok"], out
    return out["block"]


def test_gui_strings_block_shape():
    b = _gui_block(_BODY)
    assert b["recognizers"]["is_gentle"] is True
    assert b["strings"]["status"] == "complete" and b["strings"]["count"] == 5
    assert b["bands"]["exist"] is False and b["rep_type"] == "finite"
    assert b["ag_invariant"] == [[4, 2]]                        # AAG value for kA3/(ab)
    assert "butler_ringel" in b["references"]


def test_gui_matches_server(tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    server = run_spec(ComputeRequest.model_validate(_BODY), tmp_path)["results"]["strings"]
    gui = _gui_block(_BODY)
    assert json.dumps(server, sort_keys=True) == json.dumps(gui, sort_keys=True)
