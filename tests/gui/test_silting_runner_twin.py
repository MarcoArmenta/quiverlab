"""Plan 67 -- both runners (the wheel's hpc.spec core and the Pyodide GUI twin) emit the
``silting`` block byte-for-byte identical, so the browser's live silting verifier /
mutation / exploration picture cannot drift from the server report / cluster run.
Algebra-level kind (no module block); the radius,budget rides in the compute string
('silting:2,64')."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

# local k[x]/(x^2): the exploration closes complete at radius 0 (Thm 2.26).
_BODY_LOCAL = {"schema": 1,
               "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                           "relations": ["x*x"], "field": {"kind": "GF", "p": 32003, "n": 1}},
               "compute": ["silting:2,64"],
               "artifacts": {"pdf": False, "tikz": False}}

# kA2: infinite silting quiver -> honest truncation.
_BODY_KA2 = {"schema": 1,
             "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                         "relations": [], "field": {"kind": "GF", "p": 7, "n": 1}},
             "compute": ["silting:2,64"],
             "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner_silting", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("silting:2,64"))
    assert out["ok"], out
    return out["block"]


def _server_block(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"]["silting"]


def test_silting_block_matches_across_runners_local(tmp_path):
    gui = _gui_block(_BODY_LOCAL)
    ref = _server_block(_BODY_LOCAL, tmp_path)
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))


def test_silting_block_matches_across_runners_ka2(tmp_path):
    gui = _gui_block(_BODY_KA2)
    ref = _server_block(_BODY_KA2, tmp_path)
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))
