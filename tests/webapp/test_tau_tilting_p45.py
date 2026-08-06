"""The tau_tilting algebra-level compute kind (Plan 45 / C4): served by hpc.spec, mirrored
by the Pyodide twin (docs/gui/runner.py). kA2 -> 5 pairs, complete, n=2 fan with 5
chambers, the AIR four-way counts all 5; the wild 2-Kronecker trips status='budget'.
Unmarked (extras-gated dir), per the Plan-32 cross-runner ruling."""
import json
import pathlib
import tempfile


def _kA2_request(budget=512):
    return {
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
                    "relations": [], "field": {"kind": "GF", "p": 7, "n": 1}},
        "compute": [f"tau_tilting:{budget}"],
        "artifacts": {"pdf": False, "tikz": False},
    }


def test_tau_tilting_block_shape(tmp_path):
    from quiverlab.hpc.spec import parse_request, run
    req = parse_request(_kA2_request(512))
    out = run(req, tmp_path)
    b = out["results"]["tau_tilting"]
    assert b["complete"] and b["num_pairs"] == 5 and b["n"] == 2
    assert b["fan"] is not None and len(b["fan"]["chambers"]) == 5
    assert b["counts"] == {"pairs": 5, "torsion": 5, "silting": 5, "semibricks": 5}
    assert b["green_count"] == 2
    assert b["citations"]                                # AIR/DIJ/King resolved


def _load_gui_runner():
    import importlib.util
    p = ("/Users/marco/Desktop/HomologicalNetworks/quiverlab/.claude/worktrees/"
         "agent-af2763c265c1dfad0/docs/gui/runner.py")
    p = str(pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py")
    spec = importlib.util.spec_from_file_location("gui_runner_p45", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_twin_parity():
    # run the same request through the server (hpc.spec) and the Pyodide twin
    # (docs/gui/runner.py); assert the tau_tilting block is byte-identical.
    from quiverlab.hpc.spec import parse_request, run
    req_dict = _kA2_request(512)
    with tempfile.TemporaryDirectory() as d:
        server = run(parse_request(req_dict), d)["results"]["tau_tilting"]
    gr = _load_gui_runner()
    assert json.loads(gr.run_build(json.dumps(req_dict)))["ok"]
    twin = json.loads(gr.compute_one("tau_tilting:512"))["block"]
    assert (json.dumps(server, sort_keys=True, default=str)
            == json.dumps(twin, sort_keys=True, default=str))


def test_estimator_budget_is_not_a_degree():
    # the tau_tilting budget (512) must NOT be read as homological degree 512 by the tier
    # classifier -- otherwise the flagship live demo would misroute / be rejected through
    # /api/compute. _max_degree must ignore the tau_tilting budget.
    from webapp.server.estimator import _max_degree
    from webapp.server.schema import ComputeRequest
    req = ComputeRequest.model_validate(_kA2_request(512))
    assert _max_degree(req) == 0


def test_wild_budget_status(tmp_path):
    # 2-Kronecker with a small budget -> status "budget", complete False, no crash.
    from quiverlab.hpc.spec import parse_request, run
    req = parse_request({
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1, 2],
                    "arrows": {"a": [1, 2], "b": [1, 2]}, "relations": [],
                    "field": {"kind": "GF", "p": 32003, "n": 1}},
        "compute": ["tau_tilting:8"],
        "artifacts": {"pdf": False, "tikz": False}})
    b = run(req, tmp_path)["results"]["tau_tilting"]
    assert b["complete"] is False and b["status"] == "budget"
    assert b["counts"] is None and b["fan"] is None
