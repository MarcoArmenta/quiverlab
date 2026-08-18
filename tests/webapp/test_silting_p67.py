"""The silting algebra-level compute kind (Plan 67): served by hpc.spec (via the webapp
runner), mirrored byte-identically by the Pyodide twin. Local k[x]/(x^2) -> complete; kA2
-> honest truncation; the regular verdict is silting; the AI citation rides along."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")


def _body(algebra, compute):
    return {"schema": 1, "algebra": algebra, "compute": compute,
            "artifacts": {"pdf": False, "tikz": False}}


_LOCAL = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
          "relations": ["x*x"], "field": {"kind": "GF", "p": 32003, "n": 1}}
_KA2 = {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
        "relations": [], "field": {"kind": "GF", "p": 7, "n": 1}}


def _server_block(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    return run_spec(ComputeRequest.model_validate(body), tmp_path)["results"]["silting"]


def test_silting_block_local_complete(tmp_path):
    block = _server_block(_body(_LOCAL, ["silting:2,64"]), tmp_path)
    assert block["regular"]["is_silting"] is True
    assert block["exploration"]["status"] == "complete"
    assert block["exploration"]["finite_class"] == "local"
    # the AI ground-truth reference rides along (registry key in references, resolved to a
    # (bibtex_key, formatted) pair in citations).
    assert "aihara_iyama_silting" in block["references"]
    assert "AiharaIyama2012" in [k for k, _ in block["citations"]]


def test_silting_block_kA2_truncates(tmp_path):
    block = _server_block(_body(_KA2, ["silting:2,64"]), tmp_path)
    assert block["exploration"]["status"] in ("radius", "budget")
    assert block["exploration"]["finite_class"] is None
    assert block["regular"]["is_silting"] is True          # A is tilting (=> silting)


def test_twin_parity(tmp_path):
    body = _body(_LOCAL, ["silting:2,64"])
    ref = _server_block(body, tmp_path)
    spec = importlib.util.spec_from_file_location("gui_runner_p67", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    gui = json.loads(mod.compute_one("silting:2,64"))["block"]
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))
