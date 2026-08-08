"""The left_right_parts algebra-level compute kind (Plan 55 / R15): served by hpc.spec,
mirrored byte-for-byte by the Pyodide twin (docs/gui/runner.py). kA3 -> L_A = R_A = ind A
(6 indecs), empty complement, both supports = A; the rad^2=0 A5 -> complement = [S3] (the
non-empty ada laura datum), e_lambda = {1,2,3}, e_rho = {3,4,5}. Unmarked (extras-gated dir),
per the Plan-32 cross-runner ruling."""
import json
import pathlib
import tempfile


def _request(vertices, arrows, relations, budget=256):
    return {
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": vertices, "arrows": arrows,
                    "relations": relations, "field": {"kind": "QQ"}},
        "compute": [f"left_right_parts:{budget}"],
        "artifacts": {"pdf": False, "tikz": False},
    }


def _kA3_request(budget=256):
    return _request([1, 2, 3], {"a": [1, 2], "b": [2, 3]}, [], budget)


def _radsq_a5_request(budget=256):
    # 1 <- 2 <- 3 <- 4 <- 5 bound by rad^2 = 0 (all length-2 paths).
    return _request([1, 2, 3, 4, 5],
                    {"a1": [2, 1], "a2": [3, 2], "a3": [4, 3], "a4": [5, 4]},
                    ["a2*a1", "a3*a2", "a4*a3"], budget)


def _load_gui_runner():
    import importlib.util
    p = str(pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py")
    spec = importlib.util.spec_from_file_location("gui_runner_p55", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_left_right_parts_block_shape(tmp_path):
    from quiverlab.hpc.spec import parse_request, run
    out = run(parse_request(_kA3_request(256)), tmp_path)
    b = out["results"]["left_right_parts"]
    assert b["complete"] and b["n"] == 3 and b["universe_size"] == 6
    assert len(b["left"]) == len(b["right"]) == 6 and b["complement"] == []
    assert set(b["left_support"]["vertices"]) == {1, 2, 3}
    assert b["citations"]                       # act/aclv/survey resolved


def test_nonempty_complement_block(tmp_path):
    # the rad^2=0 A5: complement carries exactly S3 (dim-vector {3:1}) -- the
    # ada-with-non-empty-complement datum P61 consumes.
    from quiverlab.hpc.spec import parse_request, run
    b = run(parse_request(_radsq_a5_request(256)), tmp_path)["results"]["left_right_parts"]
    assert b["complete"] and b["universe_size"] == 9
    assert len(b["left"]) == 4 and len(b["right"]) == 4
    assert [r["name"] for r in b["complement"]] == ["S_3"]
    assert b["complement"][0]["dimvec"] == {"1": 0, "2": 0, "3": 1, "4": 0, "5": 0}
    assert set(b["left_support"]["vertices"]) == {1, 2, 3}
    assert set(b["right_support"]["vertices"]) == {3, 4, 5}


def test_twin_parity():
    # same request through docs/gui/runner.py; byte-identical block (both runners call the
    # shared left_right_parts_block builder).
    from quiverlab.hpc.spec import parse_request, run
    req_dict = _kA3_request(256)
    with tempfile.TemporaryDirectory() as d:
        server = run(parse_request(req_dict), d)["results"]["left_right_parts"]
    gr = _load_gui_runner()
    assert json.loads(gr.run_build(json.dumps(req_dict)))["ok"]
    twin = json.loads(gr.compute_one("left_right_parts:256"))["block"]
    assert (json.dumps(server, sort_keys=True, default=str)
            == json.dumps(twin, sort_keys=True, default=str))


def test_estimator_budget_is_not_a_degree():
    # the module budget (256) must NOT be read as homological degree 256 by the tier
    # classifier -- _max_degree must ignore it (like ar_quiver / tau_tilting).
    from webapp.server.estimator import _max_degree
    from webapp.server.schema import ComputeRequest
    req = ComputeRequest.model_validate(_kA3_request(256))
    assert _max_degree(req) == 0


def test_selfinjective_status(tmp_path):
    # k[x]/(x^3) -> status "unsupported", complete False, no crash, empty lists / None supports.
    from quiverlab.hpc.spec import parse_request, run
    req = {
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                    "relations": ["x*x*x"], "field": {"kind": "QQ"}},
        "compute": ["left_right_parts:256"],
        "artifacts": {"pdf": False, "tikz": False}}
    b = run(parse_request(req), tmp_path)["results"]["left_right_parts"]
    assert b["complete"] is False and b["status"] == "unsupported"
    assert b["left"] == [] and b["right"] == [] and b["left_support"] is None
