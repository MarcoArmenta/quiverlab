"""Plan 72 -- both runners (the wheel's hpc.spec core and the Pyodide GUI twin)
emit the ``split_extension`` / ``arrow_removal`` blocks byte-for-byte identical, so
the browser cannot drift from the server report / cluster run. Algebra-level kinds
(no module block); the top-degree budget rides in the compute string
('split_extension:4' / 'arrow_removal:4'). Unmarked per the Plan-32 extras-gated
ruling (a cross-runner contract, not an oracle)."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

# T(kA2) needs char > dim = 6; GF(7) is fine.
_KA2 = {"schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1, 2],
                    "arrows": {"a1": [1, 2]}, "relations": [],
                    "field": {"kind": "GF", "p": 7, "n": 1}},
        "compute": ["split_extension:4"],
        "artifacts": {"pdf": False, "tikz": False}}

# P1: a loop x (x^2=0) at 1 plus an inert bridge c:1->2.
_P1 = {"schema": 1,
       "algebra": {"kind": "quiver", "vertices": [1, 2],
                   "arrows": {"x": [1, 1], "c": [1, 2]}, "relations": ["x^2"],
                   "field": {"kind": "GF", "p": 7, "n": 1}},
       "compute": ["arrow_removal:4"],
       "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body, compute):
    spec = importlib.util.spec_from_file_location("gui_runner_p72", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one(compute))
    assert out["ok"], out
    return out["block"]


def _server_block(body, kind, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"][kind]


def test_split_extension_block_matches_across_runners(tmp_path):
    gui = _gui_block(_KA2, "split_extension:4")
    ref = _server_block(_KA2, "split_extension", tmp_path)
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))


def test_arrow_removal_block_matches_across_runners(tmp_path):
    gui = _gui_block(_P1, "arrow_removal:4")
    ref = _server_block(_P1, "arrow_removal", tmp_path)
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))
