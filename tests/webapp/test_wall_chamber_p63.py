"""The wall_chamber algebra-level compute kind (Plan 63 / R25): served by hpc.spec, mirrored
by the Pyodide twin (docs/gui/runner.py), both runners byte-identical via the shared builder
tautilting.wallchamber.wall_chamber_structure. kA2 -> 5 chambers / 3 walls, complete, render
fan2d; the wild 2-Kronecker trips status='budget' (honest bounded region). Unmarked
(extras-gated dir), per the Plan-32 cross-runner ruling."""
import json
import pathlib
import tempfile


def _kA2_request(budget=512):
    return {
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                    "relations": [], "field": {"kind": "QQ"}},
        "compute": [f"wall_chamber:{budget}"],
        "artifacts": {"pdf": False, "tikz": False},
    }


def test_wall_chamber_block_shape(tmp_path):
    from quiverlab.hpc.spec import parse_request, run
    req = parse_request(_kA2_request(512))
    out = run(req, tmp_path)
    b = out["results"]["wall_chamber"]
    assert b["complete"] and b["num_chambers"] == 5 and b["num_walls"] == 3
    assert b["render"] == "fan2d"
    assert "brustle_smith_treffinger" in b["references"]
    assert b["citations"]                                # BST/DIJ/Asai/King/KT resolved


def _load_gui_runner():
    import importlib.util
    p = str(pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py")
    spec = importlib.util.spec_from_file_location("gui_runner_p63", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_twin_parity():
    # run the same request through the server (hpc.spec) and the Pyodide twin
    # (docs/gui/runner.py); assert the wall_chamber block is byte-identical.
    from quiverlab.hpc.spec import parse_request, run
    req_dict = _kA2_request(512)
    with tempfile.TemporaryDirectory() as d:
        server = run(parse_request(req_dict), d)["results"]["wall_chamber"]
    gr = _load_gui_runner()
    assert json.loads(gr.run_build(json.dumps(req_dict)))["ok"]
    twin = json.loads(gr.compute_one("wall_chamber:512"))["block"]
    assert (json.dumps(server, sort_keys=True, default=str)
            == json.dumps(twin, sort_keys=True, default=str))


def test_estimator_budget_is_not_a_degree():
    # the wall_chamber budget (512) must NOT be read as homological degree 512 by the tier
    # classifier -- otherwise the demo would misroute / be rejected through /api/compute.
    from webapp.server.estimator import _max_degree
    from webapp.server.schema import ComputeRequest
    req = ComputeRequest.model_validate(_kA2_request(512))
    assert _max_degree(req) == 0


def test_wild_budget_status(tmp_path):
    # 2-Kronecker with a small budget -> complete False, status "budget", truncation set,
    # counts null, no crash (the honest bounded region). Budget 8: the 2-Kronecker's
    # exchange_graph is super-linear (~6s at 10), so a small budget keeps this fast.
    from quiverlab.hpc.spec import parse_request, run
    req = parse_request({
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1, 2],
                    "arrows": {"a": [1, 2], "b": [1, 2]}, "relations": [],
                    "field": {"kind": "QQ"}},
        "compute": ["wall_chamber:8"],
        "artifacts": {"pdf": False, "tikz": False}})
    b = run(req, tmp_path)["results"]["wall_chamber"]
    assert b["complete"] is False and b["status"] == "budget"
    assert b["counts"] is None and b["green_count"] is None
    assert b["truncation"] and b["num_chambers"] > 0
