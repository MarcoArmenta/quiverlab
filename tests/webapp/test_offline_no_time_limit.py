"""The offline desktop app does not stop a computation on the clock.

Marco (2026-07-30): "ComputeError: DepthLimitError: job exceeded the wall-time
cap -- we don't want this in the app, because they can leave their computer
overnight or something."

The 15-minute cap and the "too big, use the email tier" refusal exist because the
DEPLOYED server is a shared public service: the cap is DoS protection, the big-job
tier gates cost behind email verification. Neither applies on the user's own
laptop, computing for the user alone. So the offline tier:

  * runs a job to completion, however long it takes (wall = 0);
  * queues and RUNS anything the GUI can express, rather than refusing it;
  * still caps MEMORY -- taking the whole desktop down with the job is never what
    the user wanted;
  * and treats quitting the app as the cancel button, rather than restarting the
    interrupted job on every launch (which, with no wall cap, would run forever).

The deployed server's own defaults must be UNCHANGED by all of this.
"""
import pathlib

import pytest

from webapp.server.config import Config
from webapp.server.estimator import classify
from webapp.server.offline import build_offline_config, runtime_caps, _banner_lines
from webapp.server.schema import ComputeRequest

_ROOT = pathlib.Path(__file__).resolve().parents[2]
_RESOURCES = {"cores": 8, "mem_bytes": 16 * 1024 ** 3}


def _offline(tmp_path, **env):
    base = {"QLWEB_TOKEN_SECRET": "x", "QLWEB_IP_HASH_SALT": "y", **env}
    return build_offline_config(tmp_path, _RESOURCES, env=base)


def _req(dim_degree=10):
    return ComputeRequest.model_validate({
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1],
                    "arrows": {"x": [1, 1]}, "relations": ["x*x*x"],
                    "field": {"kind": "CC"}},
        "compute": ["hh_cohomology:0..%d" % dim_degree],
    })


# --------------------------------------------------------------------------- #
# The offline config itself
# --------------------------------------------------------------------------- #
def test_offline_has_no_wall_cap_but_keeps_the_memory_ceiling(tmp_path):
    cfg = _offline(tmp_path)
    assert cfg.job_wall_seconds == 0                    # 0 = run to completion
    assert cfg.job_mem_bytes == 16 * 1024 ** 3 * 4 // 5  # ...memory still bounded


def test_offline_never_refuses_a_computation_as_too_big(tmp_path):
    """On the deployed server a large request lands in "big" (email tier) or is
    rejected outright; offline it must simply queue and run."""
    cfg = _offline(tmp_path)
    assert not cfg.big_jobs_enabled                     # no SMTP relay, as before
    for dim in (3, 500, 50_000):
        for degree in (2, 10, 20):
            got = classify(dim, _req(degree), cfg)["tier"]
            assert got in ("instant", "queued"), (dim, degree, got)


def test_the_deployed_server_defaults_are_untouched():
    """The cap is real DoS protection on a shared service -- lifting it offline
    must not lift it for the public deployment."""
    cfg = Config.from_env({"QLWEB_TOKEN_SECRET": "x", "QLWEB_IP_HASH_SALT": "y"})
    assert cfg.job_wall_seconds == 900
    assert cfg.queued_ops_threshold == 500_000_000
    assert cfg.queued_max_degree == 20
    # ...and a request past those still routes to the gated tiers.
    assert classify(50_000, _req(20), cfg)["tier"] == "reject"


def test_an_explicit_env_override_still_wins(tmp_path):
    """Every offline default is a setdefault: a user who WANTS a cap can set one."""
    cfg = _offline(tmp_path, QLWEB_JOB_WALL_SECONDS="60")
    assert cfg.job_wall_seconds == 60


# --------------------------------------------------------------------------- #
# The worker honours wall <= 0
# --------------------------------------------------------------------------- #
# POSIX rlimits do not exist on Windows: `_apply_caps` logs and returns early
# there, so a test that asserts a cap WAS applied can only run where one can be.
_HAS_RLIMITS = True
try:
    import resource as _resource                       # noqa: F401
except ImportError:                                    # pragma: no cover - Windows
    _HAS_RLIMITS = False


def test_no_cpu_rlimit_is_applied_when_the_wall_is_disabled(monkeypatch):
    """wall <= 0 must not set RLIMIT_CPU -- that cap is the in-child half of what
    produced Marco's message. Asserted on every platform: nothing may be applied."""
    from webapp.worker import worker as W
    applied = []
    monkeypatch.setattr(W, "_try_setrlimit",
                        lambda which, cap, name: applied.append((name, cap)))
    W._apply_caps(0, 123)
    assert not any(name == "RLIMIT_CPU" for name, _ in applied)


@pytest.mark.skipif(not _HAS_RLIMITS, reason="no POSIX rlimits on this platform")
def test_a_positive_wall_still_applies_the_cpu_rlimit(monkeypatch):
    """The other side of the branch: disabling the cap offline must not disable it
    for a caller that asks for one."""
    from webapp.worker import worker as W
    applied = []
    monkeypatch.setattr(W, "_try_setrlimit",
                        lambda which, cap, name: applied.append((name, cap)))
    W._apply_caps(900, 123)
    assert ("RLIMIT_CPU", 900) in applied


def test_the_parent_kill_is_disarmed_when_the_wall_is_disabled():
    """The wall-time kill is what produced Marco's message; with wall <= 0 the
    parent must have NO deadline at all."""
    src = (_ROOT / "webapp" / "worker" / "worker.py").read_text(encoding="utf-8")
    assert "deadline = (time.monotonic() + wall + 5) if wall > 0 else None" in src
    assert "deadline is None or time.monotonic() < deadline" in src


def test_an_explicit_zero_on_the_row_is_not_swallowed_by_the_config(tmp_path):
    """`job.wall_seconds or cfg.job_wall_seconds` treated an explicit 0 (no limit)
    as "unset" and silently substituted the config cap. 0 is meaningful now."""
    from webapp.server.store import JobStore

    cfg = _offline(tmp_path, QLWEB_JOB_WALL_SECONDS="900")   # a cap in the config
    store = JobStore(cfg.db_path)
    store.init_schema()
    job_id = store.create_job({"schema": 1}, ip="local", wall_seconds=0,
                              mem_bytes=0)
    job = store.get_job(job_id)
    wall = cfg.job_wall_seconds if job.wall_seconds is None else job.wall_seconds
    assert wall == 0, "an explicit no-limit row must not inherit the config cap"


def test_a_long_job_is_not_killed(tmp_path, monkeypatch):
    """End to end through run_one_job: a child that outlives any plausible cap
    finishes normally when the wall is disabled."""
    import time as _time
    from webapp.worker import worker as W
    from webapp.server.store import JobStore

    cfg = _offline(tmp_path)
    store = JobStore(cfg.db_path)
    store.init_schema()
    job_id = store.create_job({"schema": 1}, ip="local")
    job = store.get_job(job_id)
    assert job.wall_seconds is None                     # anonymous row: config decides

    # A "child" that runs well past the old 15-minute cap, compressed: patch the
    # clock so the parent's loop would have fired a deadline kill long ago.
    t0 = _time.monotonic()
    monkeypatch.setattr(W.time, "monotonic", lambda: t0 + 10 ** 6)

    class _Proc:
        def __init__(self, *a, **k):
            self.calls = 0
            self.killed = False

        def start(self):
            pass

        def is_alive(self):
            self.calls += 1
            return self.calls <= 2                      # exits on its own, later

        def join(self, timeout=None):
            pass

        def terminate(self):
            self.killed = True

        def kill(self):
            self.killed = True

    procs = []

    class _Ctx:
        def Queue(self):
            class _Q:
                def get_nowait(self_inner):
                    return ("ok", None)
            return _Q()

        def Process(self, *a, **k):
            p = _Proc()
            procs.append(p)
            return p

    monkeypatch.setattr(W.mp, "get_context", lambda kind: _Ctx())
    W.run_one_job(store, cfg, job)
    assert not procs[0].killed, "the child was killed despite the wall being off"
    assert store.get_job(job_id).status == "done"


# --------------------------------------------------------------------------- #
# Quitting is the cancel button
# --------------------------------------------------------------------------- #
def test_an_interrupted_job_is_failed_offline_not_restarted(tmp_path):
    """With no wall cap, requeueing a job the user ended by quitting would run it
    again on every launch, forever. Offline it is marked failed, honestly."""
    import threading
    from webapp.server import offline as O
    from webapp.server.store import JobStore

    cfg = _offline(tmp_path)
    store = JobStore(cfg.db_path)
    store.init_schema()
    job_id = store.create_job({"schema": 1}, ip="local")
    store.mark_running(job_id)                          # as if the app was quit

    stop = threading.Event()
    stop.set()                                          # one pass, then exit
    O._worker_loop(cfg, stop)

    job = store.get_job(job_id)
    assert job.status == "failed"
    assert "app was closed" in job.error
    assert "again" in job.error                         # ...and says what to do


# --------------------------------------------------------------------------- #
# The banner must not advertise a cap that no longer exists
# --------------------------------------------------------------------------- #
def test_the_banner_says_no_time_limit(tmp_path):
    cfg = _offline(tmp_path)
    caps = runtime_caps(cfg, _RESOURCES)
    text = "\n".join(_banner_lines(8000, cfg, caps, open_hint=True))
    assert "no time limit" in text
    assert "wall 0s" not in text                        # never a nonsense "0s" cap


def test_the_client_polls_until_the_job_actually_ends():
    """The page used to stop polling after 30 minutes -- on the offline app that
    was a lie about a job that was still running."""
    src = (_ROOT / "webapp" / "static" / "gui" / "worker.js").read_text(encoding="utf-8")
    assert "30 * 60 * 1000" not in src
    assert "timed out waiting for the queued job" not in src
    assert "for (;;)" in src                            # until terminal, not a clock
