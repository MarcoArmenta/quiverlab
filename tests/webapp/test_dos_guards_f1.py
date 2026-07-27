"""Fix worker F1 -- sync/instant-tier memory-DoS guards.

Covers the four gaps that let an oversized module (esp. a Tor second module) drive
a multi-GB dense-matrix allocation on the instant tier:

  * ``sizing_dim`` now counts ``tor_target`` (Fix 1a) -- a big Tor target routes off
    the instant tier just like a big module / ext_target.
  * an absolute parse-time cap on total module dimension (Fix 1c) refuses a
    pathological module BEFORE ``_full_matrices`` allocates the n x n action
    matrices -- proven WITHOUT allocating GB (a huge ``dims`` dict is tiny; the
    matrices are never built).
  * the instant path is rate-limited (Fix 1d) with the same per-IP/queue gate the
    queued path uses.
"""
import pytest
from fastapi.testclient import TestClient

from webapp.server.app import create_app
from webapp.server.config import Config
from webapp.server.estimator import classify, sizing_dim
from webapp.server.instant import InstantRateLimiter
from webapp.server.runner import RunError, run_spec
from webapp.server.schema import ComputeRequest

from quiverlab.hpc import spec as hpc_spec

_ALG = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
        "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}}


def _tor_req(module, tor_target):
    return ComputeRequest.model_validate(
        {"schema": 2, "algebra": _ALG, "compute": ["tor:0..4"],
         "artifacts": {"pdf": False, "tikz": False},
         "module": module, "tor_target": tor_target})


# --------------------------------------------------------------------------- #
# Fix 1a -- sizing_dim counts tor_target
# --------------------------------------------------------------------------- #

def test_sizing_dim_uses_tor_target_when_larger():
    # A tiny module M (builtin, bounded by the algebra) but a huge Tor target N:
    # the job must size on N, not the dim-3 algebra.
    req = _tor_req({"builtin": {"kind": "simple", "vertex": 1}},
                   {"dims": {"1": 5000}})           # side omitted -> defaults left
    assert sizing_dim(3, req) == 5000


def test_big_tor_target_routes_off_the_instant_tier(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    small = _tor_req({"builtin": {"kind": "simple", "vertex": 1}},
                     {"builtin": {"kind": "simple", "vertex": 1, "side": "left"}})
    assert classify(sizing_dim(3, small), small, cfg)["tier"] == "instant"
    # The regression the critic confirmed: WITHOUT counting tor_target this sized 3
    # and mis-classified "instant". With the fix it sizes 400 and routes off.
    big = _tor_req({"builtin": {"kind": "simple", "vertex": 1}},
                   {"dims": {"1": 400}})
    assert classify(sizing_dim(3, big), big, cfg)["tier"] != "instant"


# --------------------------------------------------------------------------- #
# Fix 1c -- absolute parse-time module-dimension cap (no GB allocated)
# --------------------------------------------------------------------------- #

def test_oversized_module_refused_at_parse_with_typed_error(tmp_path):
    # dims sum well over the 2048 cap; the request dict is tiny (no matrices), and
    # the cap fires in parse_request -> _parse_module BEFORE _full_matrices could
    # allocate the n x n action matrix. So this asserts the refusal without ever
    # allocating multi-GB.
    req = ComputeRequest.model_validate(
        {"schema": 2, "algebra": _ALG, "compute": ["dimension_vector"],
         "artifacts": {"pdf": False, "tikz": False},
         "module": {"dims": {"1": 100000}}})
    with pytest.raises(RunError) as ei:
        run_spec(req, tmp_path)
    assert ei.value.error_type == "SpecError"
    assert "exceeds" in ei.value.message and "100000" in ei.value.message


def _parse_body(total):
    return {"schema": 2, "algebra": _ALG, "compute": ["dimension_vector"],
            "artifacts": {"pdf": False, "tikz": False},
            "module": {"dims": {"1": total}}}


def test_cap_boundary_is_exact():
    # Exactly AT the cap parses (nothing legitimate is refused); one over is refused.
    # Checked at the parse layer so neither case builds a 2048 x 2048 matrix.
    req = hpc_spec.parse_request(_parse_body(hpc_spec._MAX_MODULE_DIM))
    assert req.module.dims == {"1": hpc_spec._MAX_MODULE_DIM}
    with pytest.raises(hpc_spec.SpecError):
        hpc_spec.parse_request(_parse_body(hpc_spec._MAX_MODULE_DIM + 1))


def test_cap_is_above_the_largest_test_module():
    # Guard the choice: the cap sits well above the largest module any test/legit
    # request constructs (a dim-400 module, only in a classification test).
    assert hpc_spec._MAX_MODULE_DIM >= 2048 > 400


# --------------------------------------------------------------------------- #
# Correction #1 -- the cap RESPECTS the large-compute opt-in (hpc.allow_large).
# The webapp never sets allow_large, so it stays capped; the offline/cluster CLI
# opts in and computes big local/cluster modules unrefused.
# --------------------------------------------------------------------------- #

def _mod_body(dim, hpc=None):
    body = {"schema": 2, "algebra": _ALG, "compute": ["dimension_vector"],
            "artifacts": {"pdf": False, "tikz": False},
            "module": {"dims": {"1": dim}}}
    if hpc is not None:
        body["hpc"] = hpc
    return body


def test_allow_large_bypasses_module_cap():
    # A dim-3000 module (well over the 2048 cap) with the CLI's large-compute opt-in
    # PARSES -- no cap refusal (the critic's proven regression: it was wrongly refused
    # with "run locally: pip install quiverlab" despite allow_large + max_mem_bytes).
    req = hpc_spec.parse_request(
        _mod_body(3000, hpc={"allow_large": True, "max_mem_bytes": 200 * 1024 ** 3}))
    assert req.module.dims == {"1": 3000}
    assert req.hpc is not None and req.hpc.allow_large is True


def test_allow_large_bypasses_cap_for_tor_and_ext_targets():
    # The bypass covers the second-module targets too (ext_target / tor_target),
    # since those also drive dense-matrix allocation.
    body = {"schema": 2, "algebra": _ALG, "compute": ["tor:0..2"],
            "artifacts": {"pdf": False, "tikz": False},
            "module": {"builtin": {"kind": "simple", "vertex": 1}},
            "tor_target": {"dims": {"1": 5000}, "side": "left"},
            "hpc": {"allow_large": True}}
    req = hpc_spec.parse_request(body)
    assert req.tor_target.dims == {"1": 5000}


def test_webapp_shape_still_capped_without_allow_large():
    # The WEBAPP path (no hpc block at all) is unchanged: a dim-3000 module is refused
    # with the typed SpecError that the webapp maps to a 4xx.
    with pytest.raises(hpc_spec.SpecError) as ei:
        hpc_spec.parse_request(_mod_body(3000))
    assert "exceeds" in str(ei.value) and "3000" in str(ei.value)


def test_webapp_run_spec_refuses_dim3000_module_typed(tmp_path):
    # End-to-end webapp worker path (run_spec): the critic's dim-3000 module is
    # refused at parse time with a typed SpecError (-> 4xx), WITHOUT allocating the
    # matrix. The webapp never sets allow_large, so the cap always applies here.
    req = ComputeRequest.model_validate(
        {"schema": 2, "algebra": _ALG, "compute": ["dimension_vector"],
         "artifacts": {"pdf": False, "tikz": False},
         "module": {"dims": {"1": 3000}}})
    with pytest.raises(RunError) as ei:
        run_spec(req, tmp_path)
    assert ei.value.error_type == "SpecError" and "3000" in ei.value.message


def test_hpc_block_without_allow_large_stays_capped():
    # An hpc block that does NOT opt into large compute (allow_large defaults False)
    # keeps the cap -- max_mem_bytes alone is a deepen guard, not a cap bypass.
    with pytest.raises(hpc_spec.SpecError):
        hpc_spec.parse_request(_mod_body(3000, hpc={"max_mem_bytes": 200 * 1024 ** 3}))


def test_webapp_never_sets_allow_large():
    # Belt-and-braces: the webapp's pydantic request has no hpc field, so its
    # model_dump carries no 'hpc' key -> parse_request sees allow_large=False.
    req = ComputeRequest.model_validate(
        {"schema": 2, "algebra": _ALG, "compute": ["dimension_vector"],
         "artifacts": {"pdf": False, "tikz": False},
         "module": {"dims": {"1": 2}}})
    assert "hpc" not in req.model_dump(by_alias=True)


# --------------------------------------------------------------------------- #
# Correction #2 -- a real per-IP flood throttle on the INSTANT path (the queued
# limiter never sees instant traffic, so instant requests must be throttled
# separately or an IP can fire unbounded spawned interpreters).
# --------------------------------------------------------------------------- #

def test_instant_rate_limiter_sliding_window():
    # Unit test the limiter with a fake clock: max 2 per 10s window, per IP hash.
    now = {"t": 100.0}
    lim = InstantRateLimiter(2, 10, clock=lambda: now["t"])
    assert lim.allow("ipA") is True          # 1st
    assert lim.allow("ipA") is True          # 2nd
    assert lim.allow("ipA") is False         # 3rd -> throttled
    assert lim.allow("ipB") is True          # a different IP is independent
    now["t"] += 11                            # slide past the window
    assert lim.allow("ipA") is True          # budget refreshed


def test_instant_rate_limiter_disabled_when_nonpositive():
    lim = InstantRateLimiter(0, 60)
    assert all(lim.allow("ip") for _ in range(1000))   # never throttles


def test_instant_flood_is_throttled(tmp_path):
    # Integration: cap the instant tier at ONE per window. The first instant request
    # runs (200); the second from the SAME IP is throttled with the RateLimited 429 --
    # WITHOUT this the critic fired 250 instant requests and got 250x 200 (each a
    # fresh spawned interpreter = CPU/spawn DoS).
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_INSTANT_RATE_MAX": "1",
                           "QLWEB_INSTANT_RATE_WINDOW_SECONDS": "3600"})
    client = TestClient(create_app(cfg))
    r1 = client.post("/api/compute", json=_gf_body(["hh_cohomology:0..3"]))
    assert r1.status_code == 200 and r1.json()["tier"] == "instant", r1.text
    r2 = client.post("/api/compute", json=_gf_body(["hh_cohomology:0..3"]))
    assert r2.status_code == 429
    assert r2.json()["error_type"] == "RateLimited"


def test_instant_throttle_does_not_charge_queued_limiter(tmp_path):
    # The instant throttle must NOT double-count against the queued per-IP limiter:
    # with a generous instant budget, a normal instant request still succeeds and no
    # job row is created (instant enqueues nothing).
    import sqlite3
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))
    r = client.post("/api/compute", json=_gf_body(["hh_cohomology:0..3"]))
    assert r.status_code == 200 and r.json()["tier"] == "instant"
    rows = sqlite3.connect(cfg.db_path).execute("SELECT COUNT(*) FROM jobs").fetchone()
    assert rows[0] == 0          # instant created no queued job


# --------------------------------------------------------------------------- #
# Fix 1d -- the instant path is rate-limited
# --------------------------------------------------------------------------- #

def _gf_body(compute):
    return {"schema": 1,
            "algebra": {"kind": "family", "family": "QuantumCI",
                        "params": {"q": 1}, "field": {"kind": "GF", "p": 2, "n": 1}},
            "compute": compute, "artifacts": {"pdf": False, "tikz": False}}


def test_instant_path_is_rate_limited(tmp_path):
    # A full global queue must block the instant tier too (it spawns a capped child
    # per request). GLOBAL_QUEUE_MAX=1 with a pre-filled queue -> the instant
    # request gets the same 429 the queued path returns.
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_GLOBAL_QUEUE_MAX": "1"})
    client = TestClient(create_app(cfg))
    # Queue one job (degree over the instant cap -> queued) to fill the cap.
    r_q = client.post("/api/compute", json=_gf_body(["hh_cohomology:0..12"]))
    assert r_q.status_code == 202 and r_q.json()["tier"] == "queued"
    # Now an instant-classified request is refused with 429 (queue full).
    r_i = client.post("/api/compute", json=_gf_body(["hh_cohomology:0..3"]))
    assert r_i.status_code == 429
    assert r_i.json()["error_type"] == "RateLimited"


def test_instant_path_ok_when_within_limits(tmp_path):
    # Sanity: with headroom, the instant tier still serves normally (the rate-limit
    # gate does not spuriously refuse a first request).
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    client = TestClient(create_app(cfg))
    r = client.post("/api/compute", json=_gf_body(["hh_cohomology:0..3"]))
    assert r.status_code == 200 and r.json()["tier"] == "instant"
