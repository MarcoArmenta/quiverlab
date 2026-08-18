"""The Pyodide GUI runner (docs/gui/runner.py) serves the `skew_gentle` kind (Plan 68)
SHAPE-IDENTICALLY to the server core (quiverlab.hpc.spec): the SAME library block builder
(skewgentle.block.skew_gentle_block) reads the triple off the split algebra's construction
marker + `references`->citations, so the skew-gentle block (recognizer + split shape + dim
law + classification counts + support tau-tilting + rep-type certificate) cannot drift
between the browser GUI and the server report."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

_BODY = {"schema": 1,
         "algebra": {"kind": "family", "family": "SkewGentleAlgebra",
                     "params": {"vertices": [1, 2], "arrows": {"a": [1, 2]},
                                "relations": [], "special": [2]},
                     "field": {"kind": "QQ"}},
         "compute": ["skew_gentle"], "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner_sg_twin", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("skew_gentle"))
    assert out["ok"], out
    return out["block"]


def test_gui_skew_gentle_block_shape():
    b = _gui_block(_BODY)
    assert b["is_skew_gentle"] is True
    assert b["dim_law"]["ok"] is True and b["dim_law"]["split_dim"] == 5
    assert b["rank"] == 3                                    # |Q_0| + |Sp|
    assert b["rep_type"]["rep_finite"] is True
    assert "he_zhou_zhu" in b["references"]


def test_gui_matches_server(tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    server = run_spec(ComputeRequest.model_validate(_BODY),
                      tmp_path)["results"]["skew_gentle"]
    gui = _gui_block(_BODY)
    assert json.dumps(server, sort_keys=True) == json.dumps(gui, sort_keys=True)
