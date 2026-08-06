"""Plan 45 / C4 -- both runners (the wheel's hpc.spec core and the Pyodide GUI twin) emit
the ``tau_tilting`` block byte-for-byte identical, so the browser's live wall-and-chamber
picture cannot drift from the server report / cluster run. Algebra-level kind (no module
block); the budget rides in the compute string ('tau_tilting:512')."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

_BODY = {"schema": 1,
         "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                     "relations": [], "field": {"kind": "GF", "p": 7, "n": 1}},
         "compute": ["tau_tilting:512"],
         "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner_tt", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("tau_tilting:512"))
    assert out["ok"], out
    return out["block"]


def _server_block(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"]["tau_tilting"]


def test_tau_tilting_block_matches_across_runners(tmp_path):
    gui = _gui_block(_BODY)
    ref = _server_block(_BODY, tmp_path)
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))


def test_block_is_the_full_ka2_run(tmp_path):
    ref = _server_block(_BODY, tmp_path)
    assert ref["complete"] and ref["num_pairs"] == 5 and ref["n"] == 2
    assert ref["counts"] == {"pairs": 5, "torsion": 5, "silting": 5, "semibricks": 5}
    assert len(ref["fan"]["chambers"]) == 5
