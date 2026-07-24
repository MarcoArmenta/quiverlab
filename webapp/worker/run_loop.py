"""Long-running worker entrypoint. One container runs ``cfg.worker_processes``
concurrent poll loops (spawn context), so the single worker service scales
across cores without extra containers. Each loop polls for jobs and sweeps
retention hourly."""
from __future__ import annotations

import multiprocessing as mp
import signal
import time
from datetime import datetime, timezone

from webapp.server.config import Config, get_config
from webapp.server.store import JobStore
from webapp.worker.sweeper import sweep_once
from webapp.worker.worker import worker_tick


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


def _loop(cfg: Config) -> None:
    store = JobStore(cfg.db_path)
    store.init_schema()
    mailer = _make_mailer(cfg)
    last_sweep = 0.0
    while True:
        did = worker_tick(store, cfg, mailer=mailer)
        now = time.time()
        if now - last_sweep > 3600:
            sweep_once(store, cfg,
                       datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))
            last_sweep = now
        if not did:
            time.sleep(1)


def start_workers(cfg: Config, count: int | None = None) -> list:
    """Spawn ``count`` (default ``cfg.worker_processes``) poll loops via the
    spawn context and return the live processes. Each loop claims and runs one
    job at a time; N loops occupy N cores (the job child pins every math runtime
    to a single thread).

    The loops are NOT daemonic: ``worker.run_one_job`` runs each claimed job in
    its own resource-capped ``spawn`` child, and Python forbids a daemonic
    process from having children ("daemonic processes are not allowed to have
    children"). Lifecycle is managed explicitly instead -- ``main`` terminates
    them on SIGTERM, and callers that embed ``start_workers`` (the test) join
    them themselves."""
    count = count if count is not None else cfg.worker_processes
    ctx = mp.get_context("spawn")
    procs = [ctx.Process(target=_loop, args=(cfg,), daemon=False) for _ in range(count)]
    for p in procs:
        p.start()
    return procs


def main() -> None:
    cfg = get_config()
    JobStore(cfg.db_path).init_schema()
    procs = start_workers(cfg)

    def _shutdown(signum, _frame) -> None:
        # Graceful container stop: `docker compose stop` sends SIGTERM to this
        # PID 1. Terminate each loop (they die on the default SIGTERM), then the
        # joins below return and this process exits 0. A job child running at
        # that instant is bounded by its own wall/CPU caps.
        for p in procs:
            p.terminate()

    signal.signal(signal.SIGTERM, _shutdown)
    signal.signal(signal.SIGINT, _shutdown)
    for p in procs:
        p.join()


if __name__ == "__main__":
    main()
