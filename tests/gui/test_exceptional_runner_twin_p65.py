"""Plan 65 / R27+R28 -- both runners (the wheel's hpc.spec core and the Pyodide GUI twin)
emit the ``exceptional_sequences`` block byte-for-byte identical, so the browser cannot
drift from the server report / cluster run. Algebra-level kind (no module block); the
enumeration budget rides in the compute string ('exceptional_sequences:512')."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

_KA3 = {"schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1, 2, 3],
                    "arrows": {"a": [1, 2], "b": [2, 3]}, "relations": [],
                    "field": {"kind": "GF", "p": 7, "n": 1}},
        "compute": ["exceptional_sequences:512"],
        "artifacts": {"pdf": False, "tikz": False}}

# non-hereditary but tau-tilting-finite: k(1<->2)/rad^2 (all length-2 paths zero).
_NONHERED = {"schema": 1,
             "algebra": {"kind": "quiver", "vertices": [1, 2],
                         "arrows": {"a": [1, 2], "b": [2, 1]},
                         "relations": ["a*b", "b*a"],
                         "field": {"kind": "GF", "p": 7, "n": 1}},
             "compute": ["exceptional_sequences:512"],
             "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner_exc", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("exceptional_sequences:512"))
    assert out["ok"], out
    return out["block"]


def _server_block(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"]["exceptional_sequences"]


def test_exceptional_block_matches_across_runners_hereditary(tmp_path):
    gui = _gui_block(_KA3)
    ref = _server_block(_KA3, tmp_path)
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))


def test_exceptional_block_matches_across_runners_nonhereditary(tmp_path):
    gui = _gui_block(_NONHERED)
    ref = _server_block(_NONHERED, tmp_path)
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))
