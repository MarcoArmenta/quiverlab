"""The Pyodide GUI runner (docs/gui/runner.py) serves the Plan-63 `wall_chamber` kind
SHAPE-IDENTICALLY to the server core (quiverlab.hpc.spec): the SAME library block builder
(tautilting.wallchamber.wall_chamber_structure) + `references`->citations, so the chambers,
the D(B) walls, the counts and the fan render cannot drift between the browser GUI and the
server report."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

_KA2 = {"schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                    "relations": [], "field": {"kind": "QQ"}},
        "compute": ["wall_chamber:512"], "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body, spec_item):
    spec = importlib.util.spec_from_file_location("gui_runner_wall_chamber_twin", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one(spec_item))
    assert out["ok"], out
    return out["block"]


def test_gui_wall_chamber_block_shape():
    b = _gui_block(_KA2, "wall_chamber:512")
    assert b["complete"] and b["num_chambers"] == 5 and b["num_walls"] == 3
    assert b["render"] == "fan2d"
    rays = [w for w in b["walls"] if not w["is_full_hyperplane"]]
    assert len(rays) == 1                                    # D(P1) a ray
    assert "brustle_smith_treffinger" in b["references"]


def test_gui_matches_server(tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    server = run_spec(ComputeRequest.model_validate(_KA2), tmp_path)["results"]["wall_chamber"]
    gui = _gui_block(_KA2, "wall_chamber:512")
    assert (json.dumps(server, sort_keys=True, default=str)
            == json.dumps(gui, sort_keys=True, default=str))
