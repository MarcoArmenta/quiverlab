"""Wave 3 (v1.0.1): server-tier hardening.

Four defects the audit reproduced against the released v1.0.0 server:

1. The algebra was BUILT synchronously in the request thread before any rate
   limiter, and `validate_family` capped nothing, so one unauthenticated request
   could burn unbounded CPU in the app process. The `big` and `reject` branches
   returned without touching a limiter at all.
2. `/api/gui/probe` and `/api/gui/random-module` are `async def`, so that
   synchronous build blocked the whole event loop, not one threadpool slot.
3. `purge_pending_big()` had zero callers and `feedback` had no purge at all, so
   two tables of user-supplied PII grew without bound.
4. `/docs` + `/openapi.json` were public, `ip_hash_salt` was unguarded, and a
   deployed boot would silently accept `QLWEB_JOB_WALL_SECONDS=0` (which disarms
   both the wall kill and RLIMIT_CPU -- correct offline, never on a shared host).
"""
import os
import tempfile

import pytest

from webapp.server import config as cfgmod


def _count(cfg, table):
    """Row count read straight from the store's SQLite file (no test-only API)."""
    import sqlite3
    with sqlite3.connect(cfg.db_path) as c:
        return c.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]


def _env(**over):
    base = {
        "QLWEB_DATA_DIR": tempfile.mkdtemp(prefix="w3-"),
        "QLWEB_TOKEN_SECRET": "a-real-secret-not-the-default",
        "QLWEB_IP_HASH_SALT": "a-real-salt",
    }
    base.update(over)
    return base


# --- 3a: a cheap pre-build size guard ---------------------------------------

def test_oversized_family_is_refused_before_any_build():
    """The guard must refuse from the DECLARED size, without constructing the
    algebra -- that is the whole point (the build is what costs 32s)."""
    from webapp.server.catalog import assert_within_build_budget
    cfg = cfgmod.Config.from_env(_env())
    with pytest.raises(Exception) as exc:
        assert_within_build_budget({"family": "ExteriorAlgebra", "params": {"n": 12}}, cfg)
    assert "too large" in str(exc.value).lower() or "budget" in str(exc.value).lower()


def test_ordinary_family_passes_the_build_budget():
    from webapp.server.catalog import assert_within_build_budget
    cfg = cfgmod.Config.from_env(_env())
    assert_within_build_budget({"family": "ExteriorAlgebra", "params": {"n": 3}}, cfg)


def test_build_budget_disabled_when_zero():
    """The offline desktop app sets 0 -- the user's own machine, no refusals."""
    from webapp.server.catalog import assert_within_build_budget
    cfg = cfgmod.Config.from_env(_env(QLWEB_BUILD_MAX_DIM="0"))
    assert_within_build_budget({"family": "ExteriorAlgebra", "params": {"n": 20}}, cfg)


# --- 3c: both PII tables are actually purged --------------------------------

def test_sweep_purges_unverified_big_job_emails():
    """`purge_pending_big` existed but had ZERO callers, so a plaintext email from
    a big-job request that was never confirmed lived forever."""
    from webapp.server.store import JobStore
    from webapp.worker.sweeper import sweep_once
    cfg = cfgmod.Config.from_env(_env())
    store = JobStore(cfg.db_path)
    store.init_schema()
    store.create_pending_big(spec={"x": 1}, email="someone@example.com",
                             email_hash="h")
    assert _count(cfg, "pending_big") == 1
    sweep_once(store, cfg, "2030-01-01T00:00:00Z")
    assert _count(cfg, "pending_big") == 0, "stale unverified email survived the sweep"


def test_sweep_purges_stale_feedback_contacts():
    """feedback.contact is user-supplied and frequently an email; there was no
    DELETE for the feedback table anywhere in the codebase."""
    from webapp.server.store import JobStore
    from webapp.worker.sweeper import sweep_once
    cfg = cfgmod.Config.from_env(_env())
    store = JobStore(cfg.db_path)
    store.init_schema()
    store.create_feedback(category="bug", message="m", contact="a@b.c",
                          ip="h", job_ref=None)
    assert _count(cfg, "feedback") == 1
    sweep_once(store, cfg, "2030-01-01T00:00:00Z")
    assert _count(cfg, "feedback") == 0, "stale feedback PII survived the sweep"


# --- 3d: deployment defaults ------------------------------------------------

def test_openapi_and_docs_are_not_public_on_the_deployed_app():
    from fastapi.testclient import TestClient
    from webapp.server.app import create_app
    cfg = cfgmod.Config.from_env(_env())
    c = TestClient(create_app(cfg))
    for path in ("/docs", "/redoc", "/openapi.json"):
        assert c.get(path).status_code == 404, f"{path} is publicly served"


def test_default_ip_hash_salt_is_reported_loudly(caplog):
    """assert_production_secrets guarded token_secret but said nothing about the IP
    salt, so a boot bypassing compose silently hashed every client IP with a salt
    published in this repo. It is a WARNING and not a refusal on purpose: no
    in-process signal separates a real deployment from a test simulating one, and
    docker-compose's ${VAR:?} is the layer that can and does refuse."""
    env = _env()
    env.pop("QLWEB_IP_HASH_SALT")
    with caplog.at_level("WARNING"):
        cfgmod.assert_production_secrets(cfgmod.Config.from_env(env))
    assert any("IP_HASH_SALT" in r.message for r in caplog.records), caplog.text


def test_unbounded_wall_is_reported_loudly(caplog):
    """`QLWEB_JOB_WALL_SECONDS=0` disarms the parent deadline kill AND the child
    RLIMIT_CPU -- right for the offline app, never for a shared host."""
    cfg = cfgmod.Config.from_env(_env(QLWEB_JOB_WALL_SECONDS="0"))
    with caplog.at_level("WARNING"):
        cfgmod.assert_production_secrets(cfg)
    assert any("WALL_SECONDS" in r.message for r in caplog.records), caplog.text


def test_a_real_secret_and_wall_warn_about_nothing(caplog):
    with caplog.at_level("WARNING"):
        cfgmod.assert_production_secrets(cfgmod.Config.from_env(_env()))
    assert not [r for r in caplog.records if "QLWEB_" in r.message], caplog.text


def test_offline_config_still_allows_an_unbounded_wall():
    """The offline app must keep its no-time-limit behaviour untouched."""
    from webapp.server.offline import build_offline_config
    cfg = build_offline_config(tempfile.mkdtemp(prefix="w3-off-"),
                               {"cores": 8, "mem_bytes": 16 * 1024 ** 3},
                               env={"QLWEB_TOKEN_SECRET": "x", "QLWEB_IP_HASH_SALT": "y"})
    assert cfg.job_wall_seconds == 0


# --- 3a/3b: the guard is wired into the routes, before the build -------------

def _client(**over):
    from fastapi.testclient import TestClient
    from webapp.server.app import create_app
    return TestClient(create_app(cfgmod.Config.from_env(_env(**over))))


def _big_req():
    return {"algebra": {"kind": "family", "family": "ExteriorAlgebra",
                        "params": {"n": 14}, "field": {"kind": "GF", "p": 2}},
            "compute": ["hh_cohomology:0..2"]}


def test_oversized_compute_is_refused_cheaply_not_built():
    """Before v1.0.1 this drove an unbounded synchronous build in the app process.
    It must now come back fast, as a typed 4xx, without constructing anything."""
    import time
    c = _client()
    t0 = time.time()
    r = c.post("/api/compute", json=_big_req())
    elapsed = time.time() - t0
    assert r.status_code == 422, r.status_code
    assert r.json()["error_type"] == "TooLarge", r.json()
    assert elapsed < 5.0, f"refusal took {elapsed:.1f}s -- it built the algebra"


def test_oversized_probe_is_refused_cheaply():
    """/api/gui/probe is `async def`, so a synchronous build there froze the whole
    event loop, not just one threadpool slot."""
    import time
    c = _client()
    t0 = time.time()
    r = c.post("/api/gui/probe", json=_big_req())
    elapsed = time.time() - t0
    assert r.status_code < 500, r.status_code
    assert elapsed < 5.0, f"probe took {elapsed:.1f}s -- it built the algebra"


