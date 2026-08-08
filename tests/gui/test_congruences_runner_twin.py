"""Plan 64 / R26 -- both runners (the wheel's hpc.spec core and the Pyodide GUI twin) emit
the ``congruences`` block byte-for-byte identical, so the browser's live torsion-lattice /
Con / wide picture cannot drift from the server report / cluster run. Algebra-level kind (no
module block); the pair budget rides in the compute string ('congruences:512')."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

_BODY = {"schema": 1,
         "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                     "relations": [], "field": {"kind": "QQ"}},
         "compute": ["congruences:512"],
         "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner_cong", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("congruences:512"))
    assert out["ok"], out
    return out["block"]


def _server_block(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"]["congruences"]


def test_congruences_block_matches_across_runners(tmp_path):
    gui = _gui_block(_BODY)
    ref = _server_block(_BODY, tmp_path)
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))


def test_block_is_the_full_ka2_run(tmp_path):
    ref = _server_block(_BODY, tmp_path)
    assert ref["complete"] and ref["lattice"]["size"] == 5
    assert ref["congruences"]["size"] == 5 and ref["wide"]["size"] == 5
    assert ref["lattice"]["is_semidistributive"] is True
    assert ref["lattice"]["is_modular"] is False
