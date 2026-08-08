"""The tilted_check algebra-level compute kind (Plan 60 / R17): served by hpc.spec, mirrored
byte-for-byte by the Pyodide twin (docs/gui/runner.py). kA3 -> tilted (hereditary, type A3,
certified projective slice); kA3/rad^2 -> tilted (faithful_section_found, slice {S_2,P_2,P_3},
type A3, End_A(S) = kA3 recovered); kZ3/J^2 -> not tilted (self_injective, complete). The knit
budget rides the compute string (budget_sections stays internal); the estimator must not read it
as a homological degree. Unmarked (extras-gated dir), per the Plan-32 cross-runner ruling."""
import json
import pathlib
import tempfile


def _request(vertices, arrows, relations, budget=256):
    return {
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": vertices, "arrows": arrows,
                    "relations": relations, "field": {"kind": "QQ"}},
        "compute": [f"tilted_check:{budget}"],
        "artifacts": {"pdf": False, "tikz": False},
    }


def _kA3_request(budget=256):                       # hereditary kA3 = 1->2->3
    return _request([1, 2, 3], {"a": [1, 2], "b": [2, 3]}, [], budget)


def _kA3_radsq_request(budget=256):                 # kA3/rad^2 (3->2->1, length-2 path zero)
    return _request([1, 2, 3], {"a1": [2, 1], "a2": [3, 2]}, ["a2*a1"], budget)


def _kZ3_J2_request(budget=256):                    # the 3-cycle rad^2=0 = kZ3/J^2 (self-injective)
    return _request([1, 2, 3], {"a": [1, 2], "b": [2, 3], "c": [3, 1]},
                    ["a*b", "b*c", "c*a"], budget)


def _load_gui_runner():
    import importlib.util
    p = str(pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py")
    spec = importlib.util.spec_from_file_location("gui_runner_p60", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_tilted_check_hereditary_block(tmp_path):
    from quiverlab.hpc.spec import parse_request, run
    b = run(parse_request(_kA3_request(256)), tmp_path)["results"]["tilted_check"]
    assert b["complete"] and b["verdict"] == "tilted" and b["reason"] == "hereditary"
    assert b["hereditary_type"] == "A_3" and b["n"] == 3
    assert len(b["slice"]) == 3                      # certified projective slice A_A
    assert b["citations"]                            # liu/happel-ringel/assem resolved


def test_tilted_check_faithful_section_block(tmp_path):
    # kA3/rad^2: verdict tilted, reason faithful_section_found, slice named, type A3, H recovered.
    from quiverlab.hpc.spec import parse_request, run
    b = run(parse_request(_kA3_radsq_request(256)), tmp_path)["results"]["tilted_check"]
    assert b["complete"] and b["verdict"] == "tilted"
    assert b["reason"] == "faithful_section_found"
    assert {r["name"] for r in b["slice"]} == {"S_2", "P_2", "P_3"}
    assert b["hereditary_type"] == "A_3"
    assert b["hereditary_algebra"] is not None
    assert len(b["hereditary_algebra"]["vertices"]) == 3
    assert not b["hereditary_algebra"]["relations"]  # kA3 is hereditary
    assert b["reconstruction"]["type"] == "A_3"


def test_tilted_check_not_tilted_selfinjective(tmp_path):
    # kZ3/J^2 (self-injective): verdict not_tilted, reason self_injective, complete True (Gate S,
    # a REAL verdict -- the knit refuses self-injective, the theorem gate does not).
    from quiverlab.hpc.spec import parse_request, run
    b = run(parse_request(_kZ3_J2_request(256)), tmp_path)["results"]["tilted_check"]
    assert b["complete"] and b["verdict"] == "not_tilted" and b["reason"] == "self_injective"


def test_twin_parity():
    # same request through docs/gui/runner.py; byte-identical block (both runners call the
    # shared tilted_check_block builder).
    from quiverlab.hpc.spec import parse_request, run
    for req_dict in (_kA3_request(256), _kA3_radsq_request(256), _kZ3_J2_request(256)):
        with tempfile.TemporaryDirectory() as d:
            server = run(parse_request(req_dict), d)["results"]["tilted_check"]
        gr = _load_gui_runner()
        assert json.loads(gr.run_build(json.dumps(req_dict)))["ok"]
        twin = json.loads(gr.compute_one("tilted_check:256"))["block"]
        assert (json.dumps(server, sort_keys=True, default=str)
                == json.dumps(twin, sort_keys=True, default=str))


def test_estimator_budget_is_not_a_degree():
    # the knit budget (256) must NOT be read as homological degree 256 by the tier classifier
    # -- _max_degree must ignore it (like ar_quiver / tau_tilting / left_right_parts).
    from webapp.server.estimator import _max_degree
    from webapp.server.schema import ComputeRequest
    req = ComputeRequest.model_validate(_kA3_request(256))
    assert _max_degree(req) == 0
