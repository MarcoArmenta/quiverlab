from webapp.server.config import Config
from webapp.server.estimator import classify, decide_tier, estimate_ops, sizing_dim
from webapp.server.schema import ComputeRequest


def _req(field, compute, pdf=False, tikz=False):
    return ComputeRequest.model_validate({
        "schema": 1,
        "algebra": {"kind": "family", "family": "QuantumCI",
                    "params": {"n": 3}, "field": field},
        "compute": compute, "artifacts": {"pdf": pdf, "tikz": tikz}})


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


# --------------------------------------------------------------------------- #
# Report artifacts (artifacts.pdf) need the queued tier's persistent artifact dir:
# the instant tier discards its dir (and runs capture_reps=False), so an instant
# report request would silently return no report. A would-be-instant request that
# asks for the report is therefore upgraded to queued; big/reject are unaffected.
# --------------------------------------------------------------------------- #

def test_small_no_report_stays_instant(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    req = _req({"kind": "GF", "p": 5, "n": 1}, ["hh_cohomology:0..4"], pdf=False)
    assert decide_tier(4, req, cfg) == "instant"          # unchanged


def test_small_report_request_upgrades_to_queued(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    req = _req({"kind": "GF", "p": 5, "n": 1}, ["hh_cohomology:0..4"], pdf=True)
    info = classify(4, req, cfg)                           # same dims as instant test
    assert info["tier"] == "queued"                        # ...but the report forces queued
    assert info["reason"] == "report_artifacts"


def test_tikz_alone_stays_instant(tmp_path):
    # TikZ is written to the same discarded dir but is NOT a trigger: the canvas GUI
    # requests it on every compute, so gating on it would defeat the instant tier.
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    req = _req({"kind": "GF", "p": 5, "n": 1}, ["hh_cohomology:0..4"], tikz=True)
    assert decide_tier(4, req, cfg) == "instant"


def test_report_flag_does_not_disturb_big_or_reject(tmp_path):
    # The pdf flag only ever downgrades instant->queued; it never rescues or reroutes
    # a big/reject request. Same big (SMTP on) and reject (SMTP off) verdicts as the
    # no-report cases above.
    big_cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                               "QLWEB_SMTP_HOST": "relay", "QLWEB_SMTP_FROM": "q@e.org"})
    big = _req({"kind": "CC"}, ["hh_cohomology:0..30"], pdf=True)
    assert decide_tier(300, big, big_cfg) == "big"
    reject_cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})   # no SMTP
    reject = _req({"kind": "CC"}, ["hh_cohomology:0..30"], pdf=True)
    assert decide_tier(300, reject, reject_cfg) == "reject"
    beyond = _req({"kind": "CC"}, ["hh_cohomology:0..200"], pdf=True)
    assert classify(5000, beyond, big_cfg)["reason"] == "beyond_big_cap"


# --------------------------------------------------------------------------- #
# Module-aware sizing (Plan 26): a big module must size the job even over a small
# algebra, exactly like an oversized family request.
# --------------------------------------------------------------------------- #

def _module_req(compute, module, ext_target=None):
    body = {"schema": 2,
            "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                        "relations": ["x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
            "compute": compute, "artifacts": {"pdf": False, "tikz": False},
            "module": module}
    if ext_target is not None:
        body["ext_target"] = ext_target
    return ComputeRequest.model_validate(body)


def test_sizing_dim_uses_module_when_larger():
    req = _module_req(["dimension_vector"],
                      {"dims": {"1": 50}, "maps": {"x": [[0] * 50 for _ in range(50)]}})
    assert sizing_dim(2, req) == 50            # module (50) dwarfs the algebra (2)


def test_sizing_dim_ignores_builtin_and_falls_back():
    # A builtin pick-list is bounded by the algebra, so it adds nothing extra.
    req = _module_req(["dimension_vector"], {"builtin": {"kind": "simple", "vertex": 1}})
    assert sizing_dim(2, req) == 2


def test_sizing_dim_unchanged_without_module():
    req = _req({"kind": "GF", "p": 2, "n": 1}, ["cartan"])
    assert sizing_dim(4, req) == 4             # every existing request classifies as before


def test_big_module_routes_off_the_instant_tier(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    small = _module_req(["dimension_vector"],
                        {"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]}})
    assert classify(sizing_dim(2, small), small, cfg)["tier"] == "instant"
    big = _module_req(["ext:0..6"],
                      {"dims": {"1": 400}, "maps": {"x": [[0] * 400 for _ in range(400)]}},
                      ext_target={"builtin": {"kind": "simple", "vertex": 1}})
    assert classify(sizing_dim(2, big), big, cfg)["tier"] != "instant"


# --------------------------------------------------------------------------- #
# Plan 69 (H2): a `barcode` request on a COMMUTATIVE LADDER runs a full AR knit
# (knit-heavy), so it is upgraded instant->queued (reason="knit_heavy"); a plain
# A_n/zigzag barcode only decomposes (module-sized) and stays instant-eligible.
# --------------------------------------------------------------------------- #
def _barcode_req(vertices, arrows, relations, module):
    return ComputeRequest.model_validate({
        "schema": 2,
        "algebra": {"kind": "quiver", "vertices": vertices, "arrows": arrows,
                    "relations": relations, "field": {"kind": "QQ"}},
        "compute": ["barcode"], "artifacts": {"pdf": False, "tikz": False},
        "module": module})


def _cl3_barcode_req():
    # is_commutative_ladder recognizes by GRAPH SHAPE only, so the commutativity
    # relations are not needed to trigger the knit_heavy upgrade -- the CL(3) box-product
    # graph is enough (relations=[] keeps the request a plain string-relation schema).
    from quiverlab.families.commutative_ladder import CommutativeLadder
    from quiverlab.fields import QQ
    Q = CommutativeLadder(3, field=QQ).quiver
    return _barcode_req(list(Q.vertices),
                        {name: list(st) for name, st in Q.arrows.items()}, [],
                        {"builtin": {"kind": "simple", "vertex": Q.vertices[0]}})


def test_cl_barcode_is_knit_heavy_not_instant(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    req = _cl3_barcode_req()
    out = classify(sizing_dim(18, req), req, cfg)      # CL(3) dim 18 (would-be-instant)
    assert out["tier"] == "queued" and out["reason"] == "knit_heavy"


def test_an_barcode_stays_instant_eligible(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    req = _barcode_req([1, 2, 3, 4, 5],
                       {"a": [1, 2], "b": [2, 3], "c": [3, 4], "d": [4, 5]}, [],
                       {"dims": {"1": 1, "2": 2, "3": 1, "4": 2, "5": 1},
                        "maps": {"a": [[1], [0]], "b": [[1, 1]],
                                 "c": [[1], [0]], "d": [[1, 1]]}})
    assert decide_tier(sizing_dim(15, req), req, cfg) == "instant"
