"""Plan 28 -- the offline laptop GUI.

``serve_offline`` runs the Plan-09 webapp as a single local process with an
embedded worker: fully offline by construction (server-side compute, vendored
KaTeX, no CDN, no SMTP, no big-job tier). Everything lands under a user data dir
(``$QUIVERLAB_DATA`` or ``~/.quiverlab``); a shipped seed cache is copied in on
first run so the precomputed examples replay instantly.

Import boundary: this module lives in ``webapp/`` and pulls FastAPI/uvicorn, so it
is only importable where the ``[web]`` extra is installed -- the ``quiverlab.hpc``
``gui`` verb imports it lazily, inside the verb, for exactly that reason.

Host-resource awareness (Plan 28 addendum): the offline defaults (worker memory
cap, worker count) and the banner/footer display derive from
``quiverlab.hpc.resources.detect_resources`` when that module is present, and from
a local stdlib probe until it lands. Env/config overrides always win.
"""
from __future__ import annotations

import logging
import os
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path

from webapp.server.app import create_app
from webapp.server.config import Config
from webapp.server.estimator import human_bytes
from webapp.server.store import JobStore

_log = logging.getLogger("quiverlab_web.offline")


# --------------------------------------------------------------------------- #
# Host-resource detection (delegates to the wheel's stdlib-only helper when it
# lands; a local stdlib probe until then). A module-level seam so tests can
# monkeypatch ``webapp.server.offline.detect_resources``.
# --------------------------------------------------------------------------- #

def _detect_mem_bytes() -> int:
    """Total physical RAM in bytes via POSIX ``sysconf`` (Linux + macOS); an
    honest 4 GiB floor when the platform will not say."""
    try:
        return int(os.sysconf("SC_PAGE_SIZE")) * int(os.sysconf("SC_PHYS_PAGES"))
    except (ValueError, OSError, AttributeError):
        return 4 * 1024 ** 3


def _fallback_resources() -> dict:
    """Stdlib-only host resources, matching the sibling helper's key contract:
    ``cores``, ``mem_bytes``, ``gpus`` (0 -- we do not probe accelerators here),
    ``gpus_used`` (always False; quiverlab's engines are exact CPU computation)."""
    return {"cores": os.cpu_count() or 1, "mem_bytes": _detect_mem_bytes(),
            "gpus": 0, "gpus_used": False}


def detect_resources() -> dict:
    """Detected host resources for offline defaults + display. Delegates to
    ``quiverlab.hpc.resources.detect_resources`` (sibling A, stdlib-only,
    SLURM/cgroup-aware) when importable; otherwise a local stdlib probe. The keys
    are normalised so a partial payload never crashes the caller."""
    try:
        from quiverlab.hpc.resources import detect_resources as _dr
    except Exception:                       # module not present yet -> local probe
        return _fallback_resources()
    try:
        res = dict(_dr())
    except Exception:                       # defensive: never let detection abort startup
        _log.warning("detect_resources() failed; using the local stdlib probe",
                     exc_info=True)
        return _fallback_resources()
    return {"cores": int(res.get("cores") or 1),
            "mem_bytes": int(res.get("mem_bytes") or _detect_mem_bytes()),
            "gpus": int(res.get("gpus") or 0),
            "gpus_used": bool(res.get("gpus_used", False))}


# --------------------------------------------------------------------------- #
# Data dir + seed cache
# --------------------------------------------------------------------------- #

def resolve_data_dir(data_dir=None) -> Path:
    """The user data dir: the explicit argument, else ``$QUIVERLAB_DATA``, else
    ``~/.quiverlab``. ``~`` is expanded."""
    if data_dir:
        return Path(data_dir).expanduser()
    env = os.environ.get("QUIVERLAB_DATA")
    if env:
        return Path(env).expanduser()
    return Path.home() / ".quiverlab"


def find_seed_db() -> Path | None:
    """Locate the shipped seed cache DB: ``$QUIVERLAB_SEED_CACHE`` first, then the
    image's ``/opt/quiverlab/seed-cache.db``. Returns None when neither exists."""
    for cand in (os.environ.get("QUIVERLAB_SEED_CACHE"), "/opt/quiverlab/seed-cache.db"):
        if cand:
            p = Path(cand).expanduser()
            if p.is_file():
                return p
    return None


def seed_first_run(cfg: Config, seed_db: Path | None) -> bool:
    """On first run (no DB yet), copy the shipped seed bundle into the data dir:
    the SQLite file to ``cfg.db_path`` and its sibling ``artifacts/`` tree to
    ``cfg.artifacts_dir`` (the same bundle shape ``container/seed_cache.py``
    writes). Returns whether a seed was installed.

    Version-keyed by construction: Plan-25 cache keys embed the library version,
    so a seed built by a different quiverlab simply never hits (its rows age out
    of the LRU) -- no bespoke invalidation lives here. A no-op when there is no
    seed, or when the DB already exists (a returning user keeps their own data)."""
    if seed_db is None or not seed_db.is_file():
        return False
    if cfg.db_path.exists():
        return False
    cfg.data_dir.mkdir(parents=True, exist_ok=True)
    cfg.artifacts_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(seed_db, cfg.db_path)
    seed_artifacts = seed_db.parent / "artifacts"
    if seed_artifacts.is_dir():
        shutil.copytree(seed_artifacts, cfg.artifacts_dir, dirs_exist_ok=True)
    return True


# --------------------------------------------------------------------------- #
# Offline config + caps
# --------------------------------------------------------------------------- #

def build_offline_config(data_dir, resources: dict, *, env=None) -> Config:
    """A ``Config`` tuned for a laptop: data dir under the user home, one worker
    loop, and a worker memory cap derived from the DETECTED RAM (four fifths of
    it, leaving headroom for the OS + browser). Big-job/SMTP tiers stay off (no
    relay is configured -> ``big_jobs_enabled`` is False). Every default is a
    ``setdefault``, so an explicit ``QLWEB_*`` env var always wins."""
    base = dict(os.environ if env is None else env)
    base.setdefault("QLWEB_DATA_DIR", str(data_dir))
    mem = int(resources.get("mem_bytes") or 0)
    if mem > 0:
        base.setdefault("QLWEB_JOB_MEM_BYTES", str(mem * 4 // 5))
    base.setdefault("QLWEB_WORKER_PROCESSES", "1")   # a laptop, not a fleet
    return Config.from_env(base)


def runtime_caps(cfg: Config, resources: dict) -> dict:
    """The numbers the banner + GUI footer show: detected host resources alongside
    the configured worker caps under which the offline user actually computes."""
    mem = int(resources.get("mem_bytes") or 0)
    return {
        "cores": int(resources.get("cores") or 1),
        "ram_bytes": mem,
        "ram_human": human_bytes(mem),
        "worker_mem_bytes": cfg.job_mem_bytes,
        "worker_mem_human": human_bytes(cfg.job_mem_bytes),
        "worker_wall_seconds": cfg.job_wall_seconds,
        "workers": cfg.worker_processes,
        "gpus": int(resources.get("gpus") or 0),
        "gpus_used": bool(resources.get("gpus_used", False)),
    }


# --------------------------------------------------------------------------- #
# App factory (testable without a live server) + the blocking serve entrypoint
# --------------------------------------------------------------------------- #

def create_offline_app(data_dir=None, *, env=None, mailer=None):
    """Build the offline FastAPI app and its ``Config`` WITHOUT starting a server
    or a worker. Resolves the data dir, installs the seed cache on first run,
    builds a laptop-tuned config, and stamps ``app.state`` with the offline flag +
    detected resources + caps (the pages read these to show the footer). Returns
    ``(app, cfg, resources)``."""
    data_dir = resolve_data_dir(data_dir)
    resources = detect_resources()
    cfg = build_offline_config(data_dir, resources, env=env)
    seeded = seed_first_run(cfg, find_seed_db())
    app = create_app(cfg, mailer=mailer)
    app.state.offline = True
    app.state.resources = resources
    app.state.caps = runtime_caps(cfg, resources)
    app.state.seeded = seeded
    return app, cfg, resources


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _worker_loop(cfg: Config, stop) -> None:
    """A single in-process poll loop (background thread): claim + run one job per
    tick in a resource-capped child (``worker_tick`` -> ``run_one_job``), sweeping
    the cache + retention hourly. Cooperative shutdown via the ``stop`` event."""
    from webapp.worker.sweeper import sweep_cache_once, sweep_once
    from webapp.worker.worker import worker_tick

    store = JobStore(cfg.db_path)
    store.init_schema()
    # Adopt any jobs a previous session stranded in 'running' before we start
    # claiming (single-writer moment): this is the only loop.
    try:
        requeued = store.requeue_stale_running()
        if requeued:
            _log.info("offline: requeued %d stranded job(s)", len(requeued))
    except Exception:                        # never let a startup hiccup kill the loop
        _log.warning("offline: startup requeue failed", exc_info=True)
    last_sweep = 0.0
    while not stop.is_set():
        try:
            did = worker_tick(store, cfg)
        except Exception:
            _log.exception("offline: worker tick failed")
            did = False
        now = time.time()
        if now - last_sweep > 3600:
            try:
                sweep_cache_once(store, cfg)
                sweep_once(store, cfg, _now_iso())
            except Exception:
                _log.warning("offline: sweep failed", exc_info=True)
            last_sweep = now
        if not did and not stop.is_set():
            stop.wait(1)                     # idle backoff, promptly interruptible


def _banner_lines(port: int, cfg: Config, caps: dict, open_hint: bool) -> list[str]:
    lines = ["quiverlab -- offline GUI"]
    lines.append(f"  open  http://localhost:{port}" if open_hint
                 else f"  serving on http://127.0.0.1:{port}")
    lines.append(f"  data dir:     {cfg.data_dir}")
    lines.append(f"  host:         {caps['cores']} core(s), {caps['ram_human']} RAM detected")
    lines.append(f"  worker caps:  memory {caps['worker_mem_human']}, "
                 f"wall {caps['worker_wall_seconds']}s, {caps['workers']} worker(s)")
    if caps.get("gpus"):
        lines.append(f"  {caps['gpus']} GPU(s) detected -- not used: quiverlab's "
                     "engines are exact CPU computation")
    lines.append("  big-job/email tier: disabled (all compute is local)")
    return lines


def serve_offline(port: int = 8000, data_dir=None, open_hint: bool = True) -> None:
    """Serve the offline GUI on ``127.0.0.1:port`` (loopback only) with an
    embedded worker thread; block until interrupted. This is the target of the
    ``quiverlab-hpc gui`` verb -- imported lazily there so the ``run``/``render``
    paths carry no webapp dependency."""
    import threading

    import uvicorn

    app, cfg, resources = create_offline_app(data_dir)
    caps = app.state.caps
    # Thread hint from detected cores for any in-process compute (the per-job
    # child re-pins its own thread caps); overrides win via setdefault.
    threads = str(max(1, min(int(resources.get("cores") or 1), 8)))
    os.environ.setdefault("NUMBA_NUM_THREADS", threads)
    os.environ.setdefault("OMP_NUM_THREADS", threads)

    print("\n".join(_banner_lines(port, cfg, caps, open_hint)), flush=True)

    stop = threading.Event()
    worker = threading.Thread(target=_worker_loop, args=(cfg, stop),
                              name="ql-offline-worker", daemon=True)
    worker.start()
    try:
        # uvicorn installs its own SIGINT/SIGTERM handlers on the main thread and
        # returns after a graceful shutdown; we then drain the worker.
        uvicorn.run(app, host="127.0.0.1", port=port, log_level="warning")
    finally:
        stop.set()
        worker.join(timeout=cfg.job_wall_seconds + 5)
