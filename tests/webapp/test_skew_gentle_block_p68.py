"""The `skew_gentle` algebra-level scalar kind (Plan 68): served by hpc.spec from a
triple (family: SkewGentleAlgebra), mirrored byte-identically by the Pyodide twin
(docs/gui/runner.py). Unmarked (extras-gated dir per the Plan-32 ruling)."""
import importlib.util
import json
import pathlib

from webapp.server.runner import run_spec
from webapp.server.schema import ComputeRequest

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")


def _body(vertices, arrows, relations, special):
    return {"schema": 1,
            "algebra": {"kind": "family", "family": "SkewGentleAlgebra",
                        "params": {"vertices": vertices, "arrows": arrows,
                                   "relations": relations, "special": special},
                        "field": {"kind": "QQ"}},
            "compute": ["skew_gentle"], "artifacts": {"pdf": False, "tikz": False}}


def _server_block(body, tmp_path):
    return run_spec(ComputeRequest.model_validate(body), tmp_path)["results"]["skew_gentle"]


def _gui_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner_sg_block", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("skew_gentle"))
    assert out["ok"], out
    return out["block"]


def test_skew_gentle_block_shape(tmp_path):
    # triple: Q = 1 --a--> 2, Sp = {2}, over QQ, compute ["skew_gentle"].
    block = _server_block(_body([1, 2], {"a": [1, 2]}, [], [2]), tmp_path)
    assert block["is_skew_gentle"] is True
    assert block["dim_law"]["ok"] is True
    assert block["dim_law"]["split_dim"] == block["dim_law"]["assoc_gentle_dim"] == 5
    assert block["rank"] == 3                                # |Q_0| + |Sp|
    assert block["rep_type"]["rep_finite"] is True
    assert "he_zhou_zhu" in block["references"]              # registry key
    assert "HeZhouZhu2020" in [k for k, _ in block["citations"]]   # resolved BibTeX id


def test_sp_empty_block_matches_gentle(tmp_path):
    # Sp = empty triple -> the split byte-reduces to the plain gentle A(Q,I); the block
    # reports rank = |Q_0| and the same dim law.
    block = _server_block(_body([1, 2, 3], {"a": [1, 2], "b": [2, 3]}, ["a*b"], []),
                          tmp_path)
    assert block["is_skew_gentle"] is True
    assert block["dim_law"]["split_dim"] == block["dim_law"]["assoc_gentle_dim"]
    assert block["rank"] == 3                                # |Q_0| + 0
    assert block["classification"]["num_special"] == 0


def test_twin_parity(tmp_path):
    body = _body([1, 2], {"a": [1, 2]}, [], [2])
    server = _server_block(body, tmp_path)
    gui = _gui_block(body)
    assert json.dumps(server, sort_keys=True) == json.dumps(gui, sort_keys=True)
