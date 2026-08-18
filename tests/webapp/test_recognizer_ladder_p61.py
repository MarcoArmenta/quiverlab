"""The recognizer_ladder algebra-level compute kind (Plan 61 / R18): served by hpc.spec,
mirrored byte-for-byte by the Pyodide twin (docs/gui/runner.py). kA3 -> all five verdicts
True, empty complement, gl.dim 1 (hereditary); the rad^2=0 A5 -> (qt,shod,wshod,laura,ada) =
(F,F,T,T,T) with complement = [S_3] and the shod S_3 witness; self-injective -> status
"unsupported" (no partial ladder); over CC the ada tree fixture emits the Theorem-B verdict
(simply connected). Unmarked (extras-gated dir), per the Plan-32 cross-runner ruling."""
import json
import pathlib
import tempfile


def _request(vertices, arrows, relations, field=None, budget=256):
    return {
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": vertices, "arrows": arrows,
                    "relations": relations, "field": field or {"kind": "QQ"}},
        "compute": [f"recognizer_ladder:{budget}"],
        "artifacts": {"pdf": False, "tikz": False},
    }


def _kA3_request(budget=256):
    return _request([1, 2, 3], {"a": [1, 2], "b": [2, 3]}, [], budget=budget)


def _radsq_a5_request(field=None, budget=256):
    # 1 <- 2 <- 3 <- 4 <- 5 bound by rad^2 = 0 (all length-2 paths).
    return _request([1, 2, 3, 4, 5],
                    {"a1": [2, 1], "a2": [3, 2], "a3": [4, 3], "a4": [5, 4]},
                    ["a2*a1", "a3*a2", "a4*a3"], field=field, budget=budget)


def _load_gui_runner():
    import importlib.util
    p = str(pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py")
    spec = importlib.util.spec_from_file_location("gui_runner_p61", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _verdicts(b):
    return {k: b["ladder"][k]["verdict"] for k in
            ("quasi_tilted", "shod", "weakly_shod", "ada")} | {
        "laura": b["ladder"]["laura"]["verdict"]}


def test_recognizer_ladder_block_shape_kA3(tmp_path):
    from quiverlab.hpc.spec import parse_request, run
    b = run(parse_request(_kA3_request(256)), tmp_path)["results"]["recognizer_ladder"]
    assert b["complete"] and b["n"] == 3 and b["universe_size"] == 6 and b["gldim"] == 1
    assert _verdicts(b) == {"quasi_tilted": True, "shod": True, "weakly_shod": True,
                            "laura": True, "ada": True}
    assert b["complement"] == []
    # ada True over QQ (not algebraically closed): dim HH^1 reported, no SC verdict.
    sc = b["simple_connectedness"]
    assert sc["applicable"] is False and sc["hh1_dim"] == 0 and sc["verdict"] is None
    assert b["citations"]                       # HRS/CL/AC/ACLV/Smith/BFT/survey resolved


def test_full_five_verdict_ladder(tmp_path):
    # ACLV Example 2.2(b): the rad^2=0 A5 is (F,F,T,T,T), complement = [S_3], gl.dim 4.
    from quiverlab.hpc.spec import parse_request, run
    b = run(parse_request(_radsq_a5_request()), tmp_path)["results"]["recognizer_ladder"]
    assert b["complete"] and b["universe_size"] == 9 and b["gldim"] == 4
    assert _verdicts(b) == {"quasi_tilted": False, "shod": False, "weakly_shod": True,
                            "laura": True, "ada": True}
    assert [r["name"] for r in b["complement"]] == ["S_3"]
    assert b["ladder"]["shod"]["witness"]["module"] == "S_3"


def test_theorem_b_verdict_over_CC(tmp_path):
    # the ada tree fixture over CC: applicable True, HH^1 = 0, simply connected (Theorem B).
    from quiverlab.hpc.spec import parse_request, run
    b = run(parse_request(_radsq_a5_request(field={"kind": "CC"})),
            tmp_path)["results"]["recognizer_ladder"]
    sc = b["simple_connectedness"]
    assert sc["applicable"] is True and sc["field_algebraically_closed"] is True
    assert sc["hh1_dim"] == 0 and sc["verdict"] is True


def test_twin_parity(tmp_path):
    # same request through docs/gui/runner.py; byte-identical block (both runners call the
    # shared recognizer_ladder_block builder).
    from quiverlab.hpc.spec import parse_request, run
    req_dict = _kA3_request(256)
    with tempfile.TemporaryDirectory() as d:
        server = run(parse_request(req_dict), d)["results"]["recognizer_ladder"]
    gr = _load_gui_runner()
    assert json.loads(gr.run_build(json.dumps(req_dict)))["ok"]
    twin = json.loads(gr.compute_one("recognizer_ladder:256"))["block"]
    assert (json.dumps(server, sort_keys=True, default=str)
            == json.dumps(twin, sort_keys=True, default=str))


def test_estimator_budget_is_not_a_degree():
    # the module budget (256) must NOT be read as homological degree 256 by the tier
    # classifier -- _max_degree must ignore it (like ar_quiver / left_right_parts).
    from webapp.server.estimator import _max_degree
    from webapp.server.schema import ComputeRequest
    req = ComputeRequest.model_validate(_kA3_request(256))
    assert _max_degree(req) == 0


def test_selfinjective_status(tmp_path):
    # k[x]/(x^3) -> status "unsupported", complete False, empty ladder, no crash.
    from quiverlab.hpc.spec import parse_request, run
    req = {
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                    "relations": ["x*x*x"], "field": {"kind": "QQ"}},
        "compute": ["recognizer_ladder:256"],
        "artifacts": {"pdf": False, "tikz": False}}
    b = run(parse_request(req), tmp_path)["results"]["recognizer_ladder"]
    assert b["complete"] is False and b["status"] == "unsupported"
    assert b["ladder"] == {} and b["complement"] == []
