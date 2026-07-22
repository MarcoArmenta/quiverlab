import json

LOOP_GF2 = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
    "compute": ["hh_cohomology:0..4", "cartan"],
    "artifacts": {"pdf": False, "tikz": True},
}
KRONECKER_CC = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1, 2],
                "arrows": {"a": [1, 2], "b": [1, 2]},
                "relations": [], "field": {"kind": "CC"}},
    "compute": ["hh_cohomology:0..4"],
    "artifacts": {"pdf": False, "tikz": True},
}


def _est(runner, req, factor=1.0):
    out = json.loads(runner.run_build(json.dumps(req)))
    assert out["ok"], out
    return json.loads(runner.estimate(factor))


def test_bar_guard_cells_pinned(runner):
    # The library's own DepthLimitError printed 8748 x 2916 for Kronecker d^6.
    assert runner._bar_guard_cells(4, 6) == 8748 * 2916
    # m*(m-1)^n basis arithmetic, squared-ish guard:
    assert runner._bar_guard_cells(3, 2) == (3 * 2 ** 3) * (3 * 2 ** 2)


def test_cap_degrees_match_engine(runner):
    assert runner._cap_degree(4, 10) == 6     # Kronecker/CC caps at d^6 (observed)
    assert runner._cap_degree(3, 10) == 9     # loop m=3 runs top 8, caps at 9
    assert runner._cap_degree(9, 10) == 3     # dim-9 zoo caps at 3
    assert runner._cap_degree(3, 8) is None   # inside the cap: no cap
    assert runner._cap_degree(2, 10) is None  # m=2: 4*1^k never exceeds the cap


def test_estimate_buckets_and_breakdown(runner):
    out = _est(runner, KRONECKER_CC, factor=1.0)
    assert out["ok"] and out["dim"] == 4 and out["cap_degree"] is None
    assert out["bucket"] in ("seconds", "minute")     # ~2 s native at factor 1
    assert [b["invariant"] for b in out["breakdown"]] == ["hh_cohomology:0..4"]
    assert abs(sum(b["units"] for b in out["breakdown"]) - out["units"]) < 1e-9
    # Huge factor pushes the same request into the long bucket.
    slow = json.loads(runner.estimate(1000.0))
    assert slow["bucket"] == "long" and "Cancel" in slow["label"]


def test_estimate_cap_bucket(runner):
    req = dict(KRONECKER_CC, compute=["hh_cohomology:0..8"])
    out = _est(runner, req)
    assert out["bucket"] == "cap" and out["cap_degree"] == 6
    assert "degree 6" in out["label"]


def test_monotonicity(runner):
    def units(req):
        return _est(runner, req)["units"]
    base = dict(LOOP_GF2, compute=["hh_cohomology:0..4"])
    deeper = dict(LOOP_GF2, compute=["hh_cohomology:0..6"])
    assert units(deeper) > units(base)                  # more degrees costs more
    cc = json.loads(json.dumps(base))
    cc["algebra"]["field"] = {"kind": "CC"}
    assert units(cc) > units(base)                      # bar route dearer than fast
    more = dict(base, compute=["hh_cohomology:0..4", "center", "cartan"])
    assert units(more) > units(base)                    # extra invariants add units


def test_estimate_before_build_is_loud(runner):
    out = json.loads(runner.estimate(1.0))
    assert not out["ok"] and out["error"]["type"] == "RequestError"


def test_bucket_for_seconds_mapping(runner):
    for s, want in ((1, "seconds"), (30, "minute"), (200, "minutes"), (900, "long")):
        got = json.loads(runner.bucket_for_seconds(s))
        assert got["bucket"] == want, (s, got)
    assert "estimated:" in json.loads(runner.bucket_for_seconds(30))["label"]


def test_calibrate_sane_and_stateless(runner):
    import json as _json
    before = dict(runner._state)
    out = _json.loads(runner.calibrate())
    assert out["seconds"] > 0 and out["units"] > 0 and out["factor"] > 0
    assert out["seconds"] < 30                      # native: well under a second
    assert abs(out["factor"] - out["seconds"] / out["units"]) < 1e-9
    assert runner._state == before                  # calibration must not clobber
