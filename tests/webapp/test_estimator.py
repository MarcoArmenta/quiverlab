from webapp.server.config import Config
from webapp.server.estimator import classify, decide_tier, estimate_ops
from webapp.server.schema import ComputeRequest


def _req(field, compute):
    return ComputeRequest.model_validate({
        "schema": 1,
        "algebra": {"kind": "family", "family": "QuantumCI",
                    "params": {"n": 3}, "field": field},
        "compute": compute, "artifacts": {"pdf": False, "tikz": False}})


def test_gf_is_cheaper_than_cc():
    assert estimate_ops(10, 6, "GF") < estimate_ops(10, 6, "CC")


def test_small_gf_is_instant(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    req = _req({"kind": "GF", "p": 5, "n": 1}, ["hh_cohomology:0..4"])
    assert decide_tier(4, req, cfg) == "instant"


def test_deep_request_is_queued(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    # Degree over instant cap but within the anonymous queued band (≤20).
    req = _req({"kind": "GF", "p": 5, "n": 1}, ["hh_cohomology:0..12"])
    assert decide_tier(50, req, cfg) == "queued"


def test_high_degree_forces_queue_even_if_cheap(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    req = _req({"kind": "GF", "p": 2, "n": 1}, ["hh_cohomology:0..20"])
    assert decide_tier(2, req, cfg) == "queued"  # exceeds instant_max_degree


def test_over_anonymous_caps_needs_big_when_enabled(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_SMTP_HOST": "relay", "QLWEB_SMTP_FROM": "q@e.org"})
    req = _req({"kind": "CC"}, ["hh_cohomology:0..30"])   # deep + slow CC field
    assert decide_tier(300, req, cfg) == "big"


def test_over_anonymous_caps_rejected_when_big_disabled(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})   # no SMTP → big disabled
    req = _req({"kind": "CC"}, ["hh_cohomology:0..30"])
    assert decide_tier(300, req, cfg) == "reject"


def test_beyond_big_caps_is_reject(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_SMTP_HOST": "relay", "QLWEB_SMTP_FROM": "q@e.org"})
    req = _req({"kind": "CC"}, ["hh_cohomology:0..200"])   # over big_max_degree
    assert decide_tier(5000, req, cfg) == "reject"


def test_classify_carries_estimate_and_reason(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})   # no SMTP
    info = classify(300, _req({"kind": "CC"}, ["hh_cohomology:0..30"]), cfg)
    assert info["tier"] == "reject" and info["reason"] == "big_disabled"
    assert info["estimate"]["cells"] > 0 and info["estimate"]["minutes"] >= 1
    beyond = classify(5000, _req({"kind": "CC"}, ["hh_cohomology:0..200"]), cfg)
    assert beyond["reason"] == "beyond_big_cap"
