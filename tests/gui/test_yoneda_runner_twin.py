"""Plan 35 wave 3c -- both runners (the wheel's hpc spec core and the Pyodide GUI twin)
emit the Ext block's Yoneda ``interpretation`` key-for-key equal. The interpretation is
captured by ``ext_dims(..., interpret=True)`` in both tiers, so a class's constructed
exact sequence (middle module, connecting maps, exactness facts) cannot drift between
the server report and the browser GUI.
"""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

_BODY = {"schema": 2,
         "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                     "relations": [], "field": {"kind": "GF", "p": 7, "n": 1}},
         "module": {"builtin": {"kind": "simple", "vertex": 1}},
         "ext_target": {"builtin": {"kind": "simple", "vertex": 2}},
         "compute": ["ext:0..3"],
         "artifacts": {"pdf": False, "tikz": False}}


def _gui_ext_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("ext:0..3"))
    assert out["ok"], out
    return out["block"]


def _server_ext_block(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"]["ext"]


def test_ext_interpretation_matches_across_runners(tmp_path):
    gui = _gui_ext_block(_BODY)
    ref = _server_ext_block(_BODY, tmp_path)
    assert "interpretation" in gui and "interpretation" in ref
    assert gui["interpretation"] == ref["interpretation"]


def test_interpretation_is_the_certified_baer_extension(tmp_path):
    """The shared payload is the constructed 0 -> S_2 -> P_1 -> S_1 -> 0, certified, with
    the middle module named as the projective cover."""
    ref = _server_ext_block(_BODY, tmp_path)
    seqs = ref["interpretation"]["sequences"]
    assert list(seqs) == ["1"]                       # only Ext^1 is nonzero for kA_2
    s = seqs["1"][0]
    assert s["certified"] is True and s["kind"] == "baer"
    mid = [m for m in s["modules"] if m["role"] == "middle"][0]
    assert mid["standard"] == {"kind": "projective", "vertex": "1"}
