import json

LOOP_GF2 = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
    "compute": [], "artifacts": {"pdf": False, "tikz": True},
}


def _build(runner, req):
    return json.loads(runner.run_build(json.dumps(req)))


def test_build_loop_gf2(runner):
    out = _build(runner, LOOP_GF2)
    assert out["ok"] is True
    assert out["dim"] == 3 and out["n_vertices"] == 1 and out["n_arrows"] == 1
    assert out["algebra"].startswith("Algebra of dimension 3")


def test_build_gf_prime_power(runner):
    req = json.loads(json.dumps(LOOP_GF2))
    req["algebra"]["field"] = {"kind": "GF", "p": 2, "n": 2}   # GF(4)
    req["algebra"]["relations"] = ["x*x"]
    out = _build(runner, req)
    assert out["ok"] is True and out["dim"] == 2


def test_family_kind_points_at_plan_09(runner):
    req = {"schema": 1, "algebra": {"kind": "family", "family": "truncated_polynomial",
                                    "params": {"n": 2}, "field": {"kind": "CC"}},
           "compute": []}
    out = _build(runner, req)
    assert out["ok"] is False
    assert out["error"]["type"] == "RequestError"
    assert "Plan 09" in out["error"]["message"]


def test_bad_schema_version_rejected(runner):
    req = json.loads(json.dumps(LOOP_GF2))
    req["schema"] = 2
    out = _build(runner, req)
    assert out["ok"] is False and out["error"]["type"] == "RequestError"


def test_field_error_surfaces_verbatim(runner):
    req = json.loads(json.dumps(LOOP_GF2))
    req["algebra"]["field"] = {"kind": "GF", "p": 6, "n": 1}   # 6 is not prime
    out = _build(runner, req)
    assert out["ok"] is False
    assert out["error"]["type"] == "FieldError"
    assert out["error"]["message"]                       # library text, verbatim


def test_relation_error_surfaces_verbatim(runner):
    req = json.loads(json.dumps(LOOP_GF2))
    req["algebra"]["relations"] = ["x*nosucharrow"]
    out = _build(runner, req)
    assert out["ok"] is False
    assert out["error"]["type"].endswith("Error")
    assert out["error"]["type"] != "InternalError"       # a NAMED library error


def test_malformed_arrows_rejected(runner):
    req = json.loads(json.dumps(LOOP_GF2))
    req["algebra"]["arrows"] = {"x": [1]}                # not a [source, target] pair
    out = _build(runner, req)
    assert out["ok"] is False and out["error"]["type"] == "RequestError"


def test_verbose_forced_off(runner):
    import quiverlab
    quiverlab.verbose = True
    _build(runner, LOOP_GF2)
    assert quiverlab.verbose is False    # the GUI must never write trace files
