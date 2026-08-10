"""Plan 74 -- both runners (the wheel's hpc.spec core and the Pyodide GUI twin)
emit the ``skew_group_hh`` block byte-for-byte identical, and both build the
``SkewGroupAlgebra`` construction family identically, so the browser cannot drift
from the server report / cluster run. Unmarked per the Plan-32 extras-gated ruling
(a cross-runner contract, not an oracle)."""
import importlib.util
import json
import pathlib

RUNNER_PATH = (pathlib.Path(__file__).resolve().parents[2]
               / "docs" / "gui" / "runner.py")

# Flagship: Z/2 on k[x]/(x^2), sigma(x) = -x, over QQ.
_FLAGSHIP = {
    "schema": 1,
    "algebra": {"kind": "family", "family": "SkewGroupAlgebra",
                "params": {"vertices": [1], "arrows": {"x": [1, 1]},
                           "relations": ["x*x"],
                           "generators": [{"vertex_perm": {"1": 1},
                                           "arrow_perm": {"x": "x"},
                                           "arrow_scalars": {"x": -1}}]},
                "field": {"kind": "QQ"}},
    "compute": ["skew_group_hh:3"],
    "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body, compute):
    spec = importlib.util.spec_from_file_location("gui_runner_p74", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one(compute))
    assert out["ok"], out
    return out["block"]


def _server_block(body, kind, tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    ref = run_spec(ComputeRequest.model_validate(body), tmp_path)
    return ref["results"][kind]


def test_skew_group_hh_block_matches_across_runners(tmp_path):
    gui = _gui_block(_FLAGSHIP, "skew_group_hh:3")
    ref = _server_block(_FLAGSHIP, "skew_group_hh", tmp_path)
    assert gui["dims"] == [1, 1, 1, 1]
    assert gui["direct_dims"] == [1, 1, 1, 1]
    assert gui["agrees"] is True
    assert "stefan_hopf_galois" in [c[0] for c in gui["citations"]] or \
           "stefan_hopf_galois" in gui["references"]
    assert (json.dumps(gui, sort_keys=True, default=str)
            == json.dumps(ref, sort_keys=True, default=str))
