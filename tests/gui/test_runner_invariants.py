import json

LOOP_GF2 = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
    "compute": [], "artifacts": {"pdf": True, "tikz": True},
}
KRONECKER_CC = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1, 2],
                "arrows": {"a": [1, 2], "b": [1, 2]},
                "relations": [], "field": {"kind": "CC"}},
    "compute": [], "artifacts": {"pdf": True, "tikz": True},
}


def _ready(runner, req):
    out = json.loads(runner.run_build(json.dumps(req)))
    assert out["ok"], out
    return out


def _one(runner, spec):
    return json.loads(runner.compute_one(spec))


def test_goldens_loop_gf2(runner):
    _ready(runner, LOOP_GF2)
    hh = _one(runner, "hh_cohomology:0..2")
    assert hh["ok"] and hh["block"]["dims"] == [3, 2, 2]
    assert hh["block"]["kind"] == "HH^" and hh["block"]["top"] == 2
    assert hh["block"]["citations"] and len(hh["block"]["citations"][0]) == 2
    assert _one(runner, "cartan")["block"]["matrix"] == [[3]]
    cox = _one(runner, "coxeter_polynomial")["block"]
    assert cox["latex"] == "t + 1"
    g = _one(runner, "global_dimension")["block"]
    assert g["text"].startswith(">=") and g["exact"] is False
    assert _one(runner, "center")["block"]["dim"] == 3


def test_goldens_kronecker_cc(runner):
    _ready(runner, KRONECKER_CC)
    assert _one(runner, "hh_cohomology:0..3")["block"]["dims"] == [1, 3, 0, 0]
    assert _one(runner, "hh_homology:0..2")["block"]["dims"] == [2, 0, 0]
    assert _one(runner, "cartan")["block"]["matrix"] == [[1, 2], [0, 1]]
    assert _one(runner, "coxeter_polynomial")["block"]["latex"] == "t^{2} - 2 t + 1"
    g = _one(runner, "global_dimension")["block"]
    assert g["exact"] is True and g["value"] == 1
    assert _one(runner, "center")["block"]["dim"] == 1


def test_cartan_latex_is_pmatrix(runner):
    _ready(runner, KRONECKER_CC)
    latex = _one(runner, "cartan")["block"]["latex"]
    assert latex == r"\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}"


def test_dimension_kind_matches_server(runner):
    # Fix F1-4: the GUI runner now serves `dimension` = A.dim, at parity with the
    # server/HPC runner (quiverlab.hpc.spec._dispatch). Same `value` semantics.
    _ready(runner, LOOP_GF2)          # k[x]/(x^3): dim 3
    out = _one(runner, "dimension")
    assert out["ok"] and out["block"]["value"] == 3
    assert "citations" in out["block"]
    _ready(runner, KRONECKER_CC)      # 2-Kronecker path algebra: dim 4
    assert _one(runner, "dimension")["block"]["value"] == 4


def test_python_snippet_serves_dimension_kind(runner):
    # Correction #3: `dimension` was added to compute_one but MISSING from the
    # snippet `calls` dict, so python_snippet() KeyError'd on `calls['dimension']`.
    # Now it maps to A.dim -- exercise the whole build -> compute -> snippet path.
    req = dict(LOOP_GF2, compute=["dimension"])
    _ready(runner, req)
    assert _one(runner, "dimension")["ok"]
    snippet = runner.python_snippet()          # must not KeyError
    assert "print(A.dim)" in snippet


def test_failed_rebuild_clears_stale_tor_target(runner):
    # Fix F1-5: run_build's state reset now clears tor_target too, so a stale Tor
    # second module cannot survive a FAILED rebuild into the next computation.
    ok = dict(LOOP_GF2, tor_target={"builtin": {"kind": "simple", "vertex": 1}})
    assert json.loads(runner.run_build(json.dumps(ok)))["ok"]
    assert runner._state["tor_target"] is not None
    bad = {"schema": 1, "algebra": {"kind": "quiver", "vertices": [1],
           "arrows": {"x": [1, 1]}, "relations": [], "field": {"kind": "BOGUS"}}}
    assert not json.loads(runner.run_build(json.dumps(bad)))["ok"]
    assert runner._state["tor_target"] is None


def test_results_and_events_accumulate(runner):
    _ready(runner, KRONECKER_CC)
    _one(runner, "hh_cohomology:0..1")
    _one(runner, "cartan")
    assert [r["invariant"] for r in runner._state["results"]] == [
        "hh_cohomology:0..1", "cartan"]
    assert runner._state["events"]          # HH filled the shared event sink


def test_unknown_invariant_and_missing_range(runner):
    _ready(runner, LOOP_GF2)
    out = _one(runner, "determinant")
    assert not out["ok"] and out["error"]["type"] == "RequestError"
    out = _one(runner, "hh_cohomology")     # range is mandatory for HH
    assert not out["ok"] and out["error"]["type"] == "RequestError"


def test_degree_cap_enforced(runner):
    _ready(runner, LOOP_GF2)
    out = _one(runner, "hh_cohomology:0..11")
    assert not out["ok"] and "cap" in out["error"]["message"]


def test_compute_before_build_is_loud(runner):
    out = _one(runner, "cartan")
    assert not out["ok"] and out["error"]["type"] == "RequestError"
