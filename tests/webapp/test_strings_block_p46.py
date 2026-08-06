"""The `strings` algebra-only scalar kind (Plan 46): served by hpc.spec / the webapp
runner, mirrored by the Pyodide twin (docs/gui/runner.py), byte-identical block.

PLAN CORRECTION (documented): the plan sketch used the 2-cycle kQ/(ab,ba) as the
"reports infinite" example, but that algebra is self-injective Nakayama = rep-FINITE
(no bands); the genuine band example is the Kronecker quiver. And gentle kA3/(ab) has
5 string modules (S1,S2,S3 + [1;2] + [2;3]), not 6 -- 6 is kA3 WITHOUT the relation."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

_FIELD = {"kind": "GF", "p": 7, "n": 1}


def _gentle_a3_body():
    return {"schema": 1,
            "algebra": {"kind": "quiver", "vertices": [1, 2, 3],
                        "arrows": {"a": [1, 2], "b": [2, 3]},
                        "relations": ["a*b"], "field": _FIELD},
            "compute": ["strings"], "artifacts": {"pdf": False, "tikz": False}}


def _kronecker_body():
    return {"schema": 1,
            "algebra": {"kind": "quiver", "vertices": [1, 2],
                        "arrows": {"a": [1, 2], "b": [1, 2]},
                        "relations": [], "field": _FIELD},
            "compute": ["strings"], "artifacts": {"pdf": False, "tikz": False}}


def _server_block(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"]["strings"]


def _gui_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner_strings", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("strings"))
    assert out["ok"], out
    return out["block"]


def test_strings_block_shape_gentle_a3(tmp_path):
    b = _server_block(_gentle_a3_body(), tmp_path)
    assert b["recognizers"]["is_gentle"] is True
    assert b["strings"]["status"] == "complete" and b["strings"]["count"] == 5
    assert b["bands"]["exist"] is False and b["rep_type"] == "finite"
    assert b["ag_invariant"] is not None                       # gentle => AG present
    assert "avella_geiss" in b["references"]                    # AG citation carried
    assert b["citations"]                                       # resolved pairs present


def test_strings_block_kronecker_reports_infinite(tmp_path):
    b = _server_block(_kronecker_body(), tmp_path)
    assert b["bands"]["exist"] is True
    assert b["strings"]["status"] == "budget"
    assert b["rep_type"] == "infinite"


def test_twin_parity(tmp_path):
    body = _gentle_a3_body()
    server = _server_block(body, tmp_path)
    gui = _gui_block(body)
    assert json.dumps(server, sort_keys=True) == json.dumps(gui, sort_keys=True)
