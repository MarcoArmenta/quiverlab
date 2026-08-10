"""Plan 69 / R33 -- both runners (the wheel's hpc.spec core and the Pyodide GUI twin)
emit the ``barcode`` block byte-for-byte identical (they share the library core builder
quiverlab.modules.barcode.barcode_block), so the browser cannot drift from the server
report / cluster run. Module-side kind: the A_5 filtration H_0 module (Plan-26 block maps)."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

_A5 = {"schema": 2,
       "algebra": {"kind": "quiver", "vertices": [1, 2, 3, 4, 5],
                   "arrows": {"a": [1, 2], "b": [2, 3], "c": [3, 4], "d": [4, 5]},
                   "relations": [], "field": {"kind": "QQ"}},
       "module": {"dims": {"1": 1, "2": 2, "3": 1, "4": 2, "5": 1},
                  "maps": {"a": [[1], [0]], "b": [[1, 1]], "c": [[1], [0]], "d": [[1, 1]]}},
       "compute": ["barcode"],
       "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body):
    spec = importlib.util.spec_from_file_location("gui_runner_barcode", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one("barcode"))
    assert out["ok"], out
    return out["block"]


def _server_block(body, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"]["barcode"]


def test_barcode_block_matches_across_runners(tmp_path):
    gui = _gui_block(_A5)
    ref = _server_block(_A5, tmp_path)
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))
