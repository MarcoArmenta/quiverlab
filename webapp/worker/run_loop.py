"""Long-running worker entrypoint. One container runs ``cfg.worker_processes``
concurrent poll loops (spawn context), so the single worker service scales
across cores without extra containers. Each loop polls for jobs; loop 0 also
sweeps retention hourly."""
from __future__ import annotations

import logging
import multiprocessing as mp
import signal
import sys
import time
from datetime import datetime, timezone

from webapp.server.config import Config, get_config
from webapp.server.store import JobStore
from webapp.worker.sweeper import sweep_cache_once, sweep_once
from webapp.worker.worker import worker_tick

_log = logging.getLogger("quiverlab_web.worker")

# Parent join grace on shutdown: after signalling stop, ``main`` waits this long
# for each loop to finish its in-flight job and exit before terminating
# stragglers. Kept just BELOW the compose ``stop_grace_period`` (job wall + 30s)
# so ``main`` reaps its own children and exits cleanly before Docker's hard
# SIGKILL. Derived from the anonymous job wall cap, which bounds a normal
# in-flight job (the child also enforces its own wall/CPU caps).
_STOP_GRACE_SLACK_SECONDS = 20


def exclude_tiers_for(index: int, fleet_size: int) -> tuple[str, ...] | None:
    """The tier reservation for poll loop ``index`` in a fleet of ``fleet_size``.

    Loop 0 refuses the ``big`` tier (so one loop is always free for instant/queued
    work) ONLY when the fleet has a second loop to run big jobs; every other loop,
    and any loop in a single-loop fleet (the offline embedded worker, where the big
    tier is off anyway), claims any tier -- byte-identical to the pre-M1 behaviour.
    Returns the ``exclude_tiers`` tuple to pass to ``worker_tick`` / ``claim_next``,
    or ``None`` for 'claim anything'."""
    return ("big",) if (index == 0 and fleet_size >= 2) else None


def _make_mailer(cfg: Config):
    """Return the real SMTP mailer for big-job completion emails, or ``None``
    when the big-job tier is off (no outbound relay configured).

    The import of ``webapp.server.mail`` is deliberately lazy. That module ships
    with the big-jobs email tier (spec §17); the worker loop must import cleanly
    and run the anonymous tier even before it exists. This mirrors ``worker.py``,
    which also imports the mail module lazily and only when a big job with an
    email address actually completes -- Task 8 carries no hard dependency on the
    big-jobs tier. When SMTP is unconfigured (the default and every test), the
    branch below never runs and ``mail`` is never imported."""
    if not cfg.big_jobs_enabled:
        return None
    from webapp.server.mail import smtp_mailer
    return smtp_mailer(cfg)


def _loop(cfg: Config, index: int, stop, fleet_size: int) -> None:
    """One poll loop. Claims and runs a single job per tick; between ticks it
    checks the shared ``stop`` event so a shutdown finishes the in-flight job
    (``worker_tick`` runs it to completion) and then exits cleanly.

    Retention sweeping runs on ``index == 0`` only: one sweeper is enough, and N
    loops all racing to ``rmtree`` the same expired dirs would be wasteful (and
    log spurious errors). ``sweep_once`` itself passes an explicit ``now_iso``.

    **Tier reservation (the M1 guarantee).** The claim queue is a single shared
    FIFO, so without this loop 0 would happily run a big job and, with a fleet of
    2, two verified-email big jobs could occupy every loop and freeze the
    anonymous instant/queued tier behind a 24 h wall. When the fleet has >= 2
    loops, loop 0 EXCLUDES the ``big`` tier (`exclude_tiers=('big',)`), so one
    loop is always free to serve instant/queued work while big jobs run. A
    single-loop fleet (the offline embedded worker, where the big tier is off
    anyway) excludes nothing and claims everything -- byte-identical to before."""
    store = JobStore(cfg.db_path)
    store.init_schema()
    mailer = _make_mailer(cfg)
    # Loop 0 is the reserved instant/queued loop when the fleet is big enough to
    # spare it; every other loop (and any loop in a 1-loop fleet) claims any tier.
    exclude_tiers = exclude_tiers_for(index, fleet_size)
    last_sweep = 0.0
    while not stop.is_set():
        did = worker_tick(store, cfg, mailer=mailer, exclude_tiers=exclude_tiers)
        if index == 0:
            now = time.time()
            if now - last_sweep > 3600:
                # Cache sweep FIRST (lifts stale-version / LRU pins), then retention
                # (reclaims any now-unpinned job past the cutoff) -- Plan 25.
                sweep_cache_once(store, cfg)
                sweep_once(store, cfg,
                           datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))
                last_sweep = now
        if not did and not stop.is_set():
            time.sleep(1)


def start_workers(cfg: Config, count: int | None = None, stop=None):
    """Spawn ``count`` (default ``cfg.worker_processes``) poll loops via the
    spawn context; return ``(procs, stop)`` -- the live processes and the shared
    stop event that signals a graceful drain. A fresh event is created when the
    caller does not supply one. Each loop claims and runs one job at a time; N
    loops occupy N cores (the job child pins every math runtime to a single
    thread).

    The loops are NOT daemonic: ``worker.run_one_job`` runs each claimed job in
    its own resource-capped ``spawn`` child, and Python forbids a daemonic
    process from having children ("daemonic processes are not allowed to have
    children"). Shutdown is cooperative via ``stop`` -- see ``main``."""
    count = count if count is not None else cfg.worker_processes
    ctx = mp.get_context("spawn")
    if stop is None:
        stop = ctx.Event()
    # Pass the ACTUAL fleet size (``count``, which may be an explicit override, not
    # cfg.worker_processes) so ``_loop`` decides its tier reservation on the real
    # number of loops -- loop 0 only reserves itself for instant/queued when there
    # is a second loop to run the big tier.
    procs = [ctx.Process(target=_loop, args=(cfg, i, stop, count), daemon=False)
             for i in range(count)]
    for p in procs:
        p.start()
    return procs, stop


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    cfg = get_config()
    store = JobStore(cfg.db_path)
    store.init_schema()

    # Startup requeue (single-writer moment, BEFORE any loop claims): adopt jobs
    # stranded in `running` by a previous worker death back into the queue.
    # `restart: unless-stopped` + this requeue make an ungraceful cut recoverable.
    requeued = store.requeue_stale_running()
    if requeued:
        _log.warning("requeued %d job(s) stranded in 'running' at startup: %s",
                     len(requeued), requeued)

    procs, stop = start_workers(cfg)

    def _shutdown(signum, _frame) -> None:
        # `docker compose stop`/`kill -TERM` sends SIGTERM to this PID 1. Set the
        # stop flag; each loop finishes its in-flight job, then exits at the next
        # tick. The supervision loop below then joins them within the grace.
        _log.info("signal %d received; draining in-flight jobs then stopping",
                  signum)
        stop.set()

    signal.signal(signal.SIGTERM, _shutdown)
    signal.signal(signal.SIGINT, _shutdown)

    # Supervise until a stop is requested or a loop dies unexpectedly. An
    # unexpected death (crash, OOM-kill) makes main exit non-zero so
    # `restart: unless-stopped` restarts the whole fleet -- the startup requeue
    # above makes that safe. This is the simplest honest supervision.
    exit_code = 0
    while not stop.is_set():
        dead = [p for p in procs if not p.is_alive()]
        if dead:
            _log.error("worker loop(s) exited unexpectedly (pids %s); exiting "
                       "non-zero to trigger a fleet restart",
                       [p.pid for p in dead])
            exit_code = 1
            break
        time.sleep(1)

    # Graceful stop: signal all loops, join within a bounded grace (one job wall
    # + slack, just under the compose stop_grace_period), then terminate any
    # straggler and reap it. An ungraceful straggler kill is still recoverable
    # via the startup requeue on the next boot.
    stop.set()
    deadline = time.monotonic() + cfg.job_wall_seconds + _STOP_GRACE_SLACK_SECONDS
    for p in procs:
        p.join(max(0.0, deadline - time.monotonic()))
    for p in procs:
        if p.is_alive():
            _log.warning("worker loop %s did not drain within the grace; "
                         "terminating", p.pid)
            p.terminate()
    for p in procs:
        p.join()
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
