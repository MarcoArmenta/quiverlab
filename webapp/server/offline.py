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

# Stand-in for "no ceiling" in the offline tier's queued thresholds. A plain int
# (the config fields are ints), far above anything the GUI can express, so
# `classify` never reaches the big/reject branches on a local machine.
_OFFLINE_UNBOUNDED = 1 << 62

def build_offline_config(data_dir, resources: dict, *, env=None) -> Config:
    """A ``Config`` tuned for a laptop: data dir under the user home, one worker
    loop, and a worker memory cap derived from the DETECTED RAM (four fifths of
    it, leaving headroom for the OS + browser). Big-job/SMTP tiers stay off (no
    relay is configured -> ``big_jobs_enabled`` is False). Every default is a
    ``setdefault``, so an explicit ``QLWEB_*`` env var always wins.

    NO TIME LIMIT and NO SIZE REFUSAL (Marco 2026-07-30). The deployed server caps
    a job at 15 minutes and refuses anything past the queued thresholds, because
    it is a SHARED public service: the cap is DoS protection and the big-job tier
    gates cost behind email. None of that applies here. This is the user's own
    machine, computing for the user alone -- they may reasonably start a real
    computation and leave it running overnight, and "job exceeded the wall-time
    cap" after 15 minutes is simply the wrong answer. So the offline app disables
    the wall cap (0 = unlimited) and lifts the queued-tier thresholds so every
    request the GUI can express is QUEUED and actually run.

    The MEMORY ceiling stays: exhausting the machine's RAM would take the whole
    desktop down with it, which is never what the user wanted."""
    base = dict(os.environ if env is None else env)
    base.setdefault("QLWEB_DATA_DIR", str(data_dir))
    mem = int(resources.get("mem_bytes") or 0)
    if mem > 0:
        base.setdefault("QLWEB_JOB_MEM_BYTES", str(mem * 4 // 5))
    base.setdefault("QLWEB_WORKER_PROCESSES", "1")   # a laptop, not a fleet
    base.setdefault("QLWEB_JOB_WALL_SECONDS", "0")   # 0 = run until it finishes
    # Nothing the GUI can express should be refused as "too big" on the user's own
    # hardware; past the instant threshold it simply queues and runs.
    base.setdefault("QLWEB_QUEUED_OPS_THRESHOLD", str(_OFFLINE_UNBOUNDED))
    base.setdefault("QLWEB_QUEUED_MAX_DEGREE", str(_OFFLINE_UNBOUNDED))
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


_SHUTDOWN_GRACE_SECONDS = 5


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
    # claiming (single-writer moment): this is the only loop. Offline they are
    # marked FAILED, not requeued: with no wall cap (see build_offline_config) a
    # job the user ended by quitting the app would otherwise restart on every
    # launch and run forever. On the server, requeueing is right -- a deploy must
    # not lose queued work -- but here quitting IS the cancel button.
    try:
        for job_id in store.requeue_stale_running():
            store.mark_failed(job_id, "Cancelled: the app was closed while this "
                                      "job was running. Submit it again to retry.")
            _log.info("offline: job %s was interrupted by a previous exit", job_id)
    except Exception:                        # never let a startup hiccup kill the loop
        _log.warning("offline: startup adoption failed", exc_info=True)
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


def _banner_lines(port: int, cfg: Config, caps: dict, open_hint: bool,
                  host: str = "127.0.0.1") -> list[str]:
    lines = ["quiverlab -- offline GUI"]
    lines.append(f"  open  http://localhost:{port}" if open_hint
                 else f"  serving on http://{host}:{port}")
    lines.append(f"  data dir:     {cfg.data_dir}")
    lines.append(f"  host:         {caps['cores']} core(s), {caps['ram_human']} RAM detected")
    wall = caps["worker_wall_seconds"]
    wall_txt = f"wall {wall}s" if wall > 0 else "no time limit"
    lines.append(f"  worker caps:  memory {caps['worker_mem_human']}, "
                 f"{wall_txt}, {caps['workers']} worker(s)")
    if caps.get("gpus"):
        lines.append(f"  {caps['gpus']} GPU(s) detected -- not used: quiverlab's "
                     "engines are exact CPU computation")
    lines.append("  big-job/email tier: disabled (all compute is local)")
    return lines


def serve_offline(port: int = 8000, data_dir=None, open_hint: bool = True,
                  host: str = "127.0.0.1") -> None:
    """Serve the offline GUI on ``host:port`` (loopback by default) with an
    embedded worker thread; block until interrupted. This is the target of the
    ``quiverlab-hpc gui`` verb -- imported lazily there so the ``run``/``render``
    paths carry no webapp dependency. Inside a container the loopback default is
    unreachable through ``docker run -p``, so the image sets
    ``QUIVERLAB_GUI_HOST=0.0.0.0`` and the verb passes it through here."""
    import threading

    import uvicorn

    app, cfg, resources = create_offline_app(data_dir)
    caps = app.state.caps
    # Thread hint from detected cores for any in-process compute (the per-job
    # child re-pins its own thread caps); overrides win via setdefault.
    threads = str(max(1, min(int(resources.get("cores") or 1), 8)))
    os.environ.setdefault("NUMBA_NUM_THREADS", threads)
    os.environ.setdefault("OMP_NUM_THREADS", threads)

    print("\n".join(_banner_lines(port, cfg, caps, open_hint, host=host)),
          flush=True)

    stop = threading.Event()
    worker = threading.Thread(target=_worker_loop, args=(cfg, stop),
                              name="ql-offline-worker", daemon=True)
    worker.start()
    try:
        # uvicorn installs its own SIGINT/SIGTERM handlers on the main thread and
        # returns after a graceful shutdown; we then drain the worker.
        uvicorn.run(app, host=host, port=port, log_level="warning")
    finally:
        stop.set()
        # A fixed grace, NOT the wall (which is unlimited offline): the loop is a
        # daemon thread and its in-flight child is killed with the process, so
        # waiting out a whole computation on Ctrl-C would be worse than leaving.
        worker.join(timeout=_SHUTDOWN_GRACE_SECONDS)
