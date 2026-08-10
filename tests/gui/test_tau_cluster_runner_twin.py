"""Plan 66 / R29 -- both runners (the wheel's hpc.spec core and the Pyodide GUI twin) emit the
``tau_cluster`` block byte-for-byte identical, so the browser cannot drift from the server report
/ cluster run. Algebra-level kind (no module block); the pair budget rides in the compute string
('tau_cluster:512'). Shared builder: tautilting.cluster_morphism.tau_cluster_block."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

# kA2 (1->2, no relations) over QQ -- the fast pin (3 gens, 1 atom, face vector (5,11,5)).
_KA2 = {"schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                    "relations": [], "field": {"kind": "QQ"}},
        "compute": ["tau_cluster:512"],
        "artifacts": {"pdf": False, "tikz": False}}

# cyclic self-injective Nakayama kZ3/rad^2 -- the K(pi,1) + 3+3 discriminator pin.
_KZ3 = {"schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1, 2, 3],
                    "arrows": {"a": [1, 2], "b": [2, 3], "c": [3, 1]},
                    "relations": ["a*b", "b*c", "c*a"], "field": {"kind": "QQ"}},
        "compute": ["tau_cluster:512"],
        "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner_tcl", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("tau_cluster:512"))
    assert out["ok"], out
    return out["block"]


def _server_block(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"]["tau_cluster"]


def test_tau_cluster_block_matches_across_runners_kA2(tmp_path):
    gui = _gui_block(_KA2)
    ref = _server_block(_KA2, tmp_path)
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))


def test_tau_cluster_block_matches_across_runners_nakayama(tmp_path):
    gui = _gui_block(_KZ3)
    ref = _server_block(_KZ3, tmp_path)
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))
