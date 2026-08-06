"""Feature 2 -- local-only "cancel the running job".

Marco: on the "you already have a job running" rate-limit, the OFFLINE desktop
app offers to cancel that running job; the deployed cloud server must not.

The cancel route (``POST /api/cancel-running``) is registered ONLY by
``webapp.server.offline.create_offline_app`` -- so the deployed ``create_app``
never carries it and a deployed request is a clean 404 (the SECURITY PIN below,
the whole point of the feature). The worker parent loop
(``webapp/worker/worker.py::run_one_job``) polls a ``cancel_requests`` signal and
kills a running child promptly, marking the job ``failed`` -- a terminal state the
read-only /draw poller (``webapp/static/gui/worker.js``) recognises.

Offline app is built without a live server (``create_offline_app`` +
``TestClient``), as in tests/webapp/test_offline.py. The worker-kills-the-child
integration uses the deterministic fake-context pattern from
tests/webapp/test_offline_no_time_limit.py (no real spawn child, no timing race).
"""
import queue

import pytest
from fastapi.testclient import TestClient

from webapp.server import offline as offline_mod
from webapp.server.app import create_app
from webapp.server.config import Config
from webapp.server.offline import create_offline_app
from webapp.server.security import hash_ip
from webapp.server.store import CANCEL_ERROR, JobStore

FAKE_RES = {"cores": 4, "mem_bytes": 8 * 1024 ** 3, "gpus": 0, "gpus_used": False}

CALLER_IP = "5.6.7.8"                       # sent as X-Forwarded-For (last hop)
HDR = {"x-forwarded-for": CALLER_IP}


@pytest.fixture()
def offline_app(monkeypatch, tmp_path):
    """A laptop-tuned offline app + its cfg, isolated from ambient seed/data env."""
    monkeypatch.setattr(offline_mod, "detect_resources", lambda: dict(FAKE_RES))
    monkeypatch.delenv("QUIVERLAB_SEED_CACHE", raising=False)
    monkeypatch.delenv("QUIVERLAB_DATA", raising=False)
    app, cfg, _res = create_offline_app(data_dir=tmp_path, env={})
    return app, cfg


def _store(cfg):
    s = JobStore(cfg.db_path)
    s.init_schema()
    return s


def _caller_hash(cfg):
    return hash_ip(CALLER_IP, cfg.ip_hash_salt)


# --------------------------------------------------------------------------- #
# The route exists on the offline app + client detection flag
# --------------------------------------------------------------------------- #
def test_offline_app_registers_the_cancel_route(offline_app):
    app, _cfg = offline_app
    paths = {getattr(r, "path", None) for r in app.routes}
    assert "/api/cancel-running" in paths


def test_offline_body_flag_lets_the_client_detect_offline(offline_app, tmp_path):
    # gui.js's offlineApp() reads <body data-offline>; base.html sets it from
    # app.state.offline. True offline, false on the deployed server.
    app, _cfg = offline_app
    assert 'data-offline="true"' in TestClient(app).get("/draw").text
    dep = create_app(Config.from_env({"QLWEB_DATA_DIR": str(tmp_path / "dep")}))
    assert 'data-offline="false"' in TestClient(dep).get("/draw").text


# --------------------------------------------------------------------------- #
# Cancelling the caller's job: running (signalled), pending (finalised), noop
# --------------------------------------------------------------------------- #
def test_cancel_running_job_records_the_signal(offline_app):
    app, cfg = offline_app
    store = _store(cfg)
    jid = store.create_job({"schema": 1}, ip=_caller_hash(cfg))
    store.mark_running(jid)                          # as if the worker claimed it
    r = TestClient(app).post("/api/cancel-running", json={}, headers=HDR)
    assert r.status_code == 200, r.text
    assert r.json() == {"cancelled": True, "job_id": jid}
    # A running job stays running until the worker acts; the SIGNAL is recorded so
    # the worker parent loop will kill the child on its next tick.
    assert store.cancel_requested(jid) is True
    assert store.get_job(jid).status == "running"


def test_cancel_pending_job_is_finalised_directly(offline_app):
    app, cfg = offline_app
    store = _store(cfg)
    jid = store.create_job({"schema": 1}, ip=_caller_hash(cfg))   # pending, unclaimed
    r = TestClient(app).post("/api/cancel-running", json={}, headers=HDR)
    assert r.status_code == 200 and r.json() == {"cancelled": True, "job_id": jid}
    job = store.get_job(jid)
    assert job.status == "failed" and job.error == CANCEL_ERROR   # terminal at once
    assert store.cancel_requested(jid) is False                  # no worker needed
    assert store.count_running_for_ip(_caller_hash(cfg)) == 0     # slot freed


def test_cancel_is_a_noop_when_nothing_is_running(offline_app):
    app, _cfg = offline_app
    r = TestClient(app).post("/api/cancel-running", json={}, headers=HDR)
    assert r.status_code == 200                       # 200 no-op, not 409
    assert r.json() == {"cancelled": False, "job_id": None}


def test_cancel_only_touches_the_callers_own_job(offline_app):
    app, cfg = offline_app
    store = _store(cfg)
    other = store.create_job({"schema": 1}, ip=hash_ip("9.9.9.9", cfg.ip_hash_salt))
    store.mark_running(other)
    r = TestClient(app).post("/api/cancel-running", json={}, headers=HDR)
    assert r.json() == {"cancelled": False, "job_id": None}       # caller has nothing
    assert store.get_job(other).status == "running"               # untouched
    assert store.cancel_requested(other) is False


# --------------------------------------------------------------------------- #
# The worker actually kills the running child and reaches a terminal state
# --------------------------------------------------------------------------- #
def test_worker_kills_the_running_child_on_cancel(offline_app, monkeypatch):
    """End to end through run_one_job with the deterministic fake-context probe
    (tests/webapp/test_offline_no_time_limit.py pattern): a child that would run
    forever is killed on the FIRST loop tick once the cancel signal is present,
    the job ends ``failed``, and the signal is torn down. No orphaned running row,
    no timing race."""
    from webapp.worker import worker as W

    _app, cfg = offline_app
    store = _store(cfg)
    jid = store.create_job({"schema": 1}, ip=_caller_hash(cfg))
    store.mark_running(jid)
    assert store.request_cancel_for_ip(_caller_hash(cfg),
                                       "2026-01-01T00:00:00Z") == jid
    assert store.cancel_requested(jid) is True

    class _Proc:
        def __init__(self, *a, **k):
            self.started = self.terminated = self.killed = False
            self._alive = True

        def start(self):
            self.started = True

        def is_alive(self):
            return self._alive

        def join(self, timeout=None):
            pass

        def terminate(self):
            self.terminated = True
            self._alive = False               # SIGTERM ends it

        def kill(self):
            self.killed = True
            self._alive = False

    procs = []

    class _Q:
        def get_nowait(self):                 # never reached: we cancel first
            raise queue.Empty

    class _Ctx:
        def Queue(self):
            return _Q()

        def Process(self, *a, **k):
            p = _Proc()
            procs.append(p)
            return p

    monkeypatch.setattr(W.mp, "get_context", lambda *a, **k: _Ctx())
    W.run_one_job(store, cfg, store.get_job(jid))

    assert procs, "run_one_job never started a child"
    assert procs[0].started is True
    assert procs[0].terminated is True        # the child was killed...
    assert procs[0].is_alive() is False       # ...and is gone
    job = store.get_job(jid)
    assert job.status == "failed"             # terminal state worker.js recognises
    assert job.error == CANCEL_ERROR
    assert store.cancel_requested(jid) is False   # signal torn down after acting


# --------------------------------------------------------------------------- #
# THE SECURITY PIN: the deployed server has NO such route (404, not 403/405)
# --------------------------------------------------------------------------- #
def test_deployed_server_has_no_cancel_route_404(tmp_path):
    """The whole point of the feature: cancelling does NOT exist on the deployed
    cloud server. Not 403 (would imply the route exists but is forbidden), not
    405 (would imply a colliding route pattern) -- a clean 404: the route is
    absent. This is why the path is /api/cancel-running, not
    /api/jobs/cancel-running (the latter collides with GET /api/jobs/{job_id} and
    would 405)."""
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    r = TestClient(create_app(cfg)).post("/api/cancel-running", json={})
    assert r.status_code == 404, (
        "deployed create_app MUST NOT expose the offline cancel route")
    # Belt and suspenders: it is genuinely absent from the route table.
    paths = {getattr(route, "path", None) for route in create_app(cfg).routes}
    assert "/api/cancel-running" not in paths
