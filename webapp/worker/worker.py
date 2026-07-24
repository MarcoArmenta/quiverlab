"""Worker: claim one pending job, run it in a resource-capped child process,
record the outcome. The child caps threads (NUMBA/OMP/BLAS = 1) and sets
RLIMIT_CPU (Linux + macOS) and RLIMIT_AS (Linux) *before* it triggers the heavy
library code, so a runaway job cannot take down the VM; it streams per-item
checkpoints to a progress file the parent drains into the store. Multiprocessing
lives at the app layer only -- the library core stays single-process.

The child uses the ``spawn`` start method: a fresh interpreter that imports
nothing heavy until ``run_spec`` runs, so the thread caps set at the top of
``_child`` are in place before numba/BLAS ever initialise their pools (quiverlab
imports neither at module load -- numba is a lazy import inside the engine)."""
from __future__ import annotations

import json
import logging
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

from webapp.server.config import Config
from webapp.server.runner import run_spec, RunError
from webapp.server.schema import ComputeRequest
from webapp.server.store import Job, JobStore

_log = logging.getLogger("quiverlab_web.worker")

# One worker child == one job. Pin every math runtime to a single thread so N
# concurrent workers occupy N cores rather than N x (core count). Set in the
# child before run_spec triggers the lazy numba/BLAS import.
_THREAD_ENV = {
    "NUMBA_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
}


def _try_setrlimit(which: int, cap: int, name: str) -> None:
    """Apply one rlimit, degrading LOUDLY (a warning log) if the platform
    refuses -- never a silent no-op."""
    import resource

    try:
        resource.setrlimit(which, (cap, cap))
    except (ValueError, OSError) as exc:  # platform refused the cap
        _log.warning("worker child could not apply %s=%d (%s: %s); the parent "
                     "wall-time kill is the remaining guard", name, cap,
                     type(exc).__name__, exc)


def _apply_caps(wall_seconds: int, mem_bytes: int) -> None:
    """Cap the child's CPU-seconds and address space.

    - RLIMIT_CPU (SIGXCPU -> SIGKILL) is enforced on both Linux and macOS.
    - RLIMIT_AS is a hard address-space ceiling only on Linux. On macOS a hard
      RLIMIT_AS is either ignored or actively breaks the child (LLVM/numba
      reserve huge virtual regions up front), so there we log LOUDLY and rely on
      the parent wall-time kill instead of silently pretending the VM is capped.

    Production runs on Linux (Docker), where both caps are hard. The test probe
    (Linux-only) asserts RLIMIT_AS lands at the requested value."""
    import resource

    _try_setrlimit(resource.RLIMIT_CPU, wall_seconds, "RLIMIT_CPU")
    if sys.platform == "linux":
        _try_setrlimit(resource.RLIMIT_AS, mem_bytes, "RLIMIT_AS")
    else:
        _log.warning("RLIMIT_AS is not enforced on %s: the %d-byte memory cap is "
                     "advisory here and the parent wall-time kill is the only "
                     "backstop (production runs on Linux where RLIMIT_AS is hard)",
                     sys.platform, mem_bytes)


def _child(spec_dict: dict, artifact_dir: str, result_max_bytes: int,
           wall_seconds: int, mem_bytes: int, q: "mp.Queue") -> None:
    # Thread caps FIRST (before run_spec triggers the lazy numba/BLAS import),
    # then the rlimit caps (per-job: anonymous vs big).
    os.environ.update(_THREAD_ENV)
    _apply_caps(wall_seconds, mem_bytes)
    # cwd = artifact dir so any stray library ./quiverlab_traces/ lands inside the
    # job dir, never at the repo root. The runner already pins verbose=False and
    # writes worked steps to an explicit out_dir, so this is belt-and-suspenders.
    os.chdir(artifact_dir)
    progress_path = Path(artifact_dir) / "progress.json"

    def _cb(d: dict) -> None:
        tmp = progress_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(d))
        tmp.replace(progress_path)          # atomic swap: parent never sees a partial

    try:
        req = ComputeRequest.model_validate(spec_dict)
        run_spec(req, Path(artifact_dir), progress_cb=_cb,
                 result_max_bytes=result_max_bytes)
        q.put(("ok", None))
    except RunError as exc:
        q.put(("fail", f"{exc.error_type}: {exc.message}"))
    except MemoryError:
        q.put(("fail", "DepthLimitError: job exceeded the memory cap"))
    except Exception as exc:  # pragma: no cover - defensive
        q.put(("fail", f"{type(exc).__name__}: {exc}"))


def _drain_progress(store: JobStore, job_id: str, progress_path: Path) -> None:
    try:
        if progress_path.exists():
            store.update_progress(job_id, json.loads(progress_path.read_text()))
    except (json.JSONDecodeError, OSError):
        pass                                # transient partial read; retry next tick


def _notify_if_big(store: JobStore, cfg: Config, job: Job, status: str,
                   mailer) -> None:
    # Completion email for big jobs only; then delete the plaintext address
    # (spec §17). Anonymous jobs are never touched. The mail module is a lazy
    # import so Task 8 carries no hard dependency on the big-jobs tier (Task 13).
    if job.tier != "big" or not job.email:
        return
    try:
        from webapp.server.mail import notify_completion
        notify_completion(cfg, job, status, mailer=mailer)
    except Exception:  # pragma: no cover - a mail failure must not lose the result
        _log.warning("completion mail failed for job %s", job.id, exc_info=True)
    finally:
        store.clear_email(job.id)


def run_one_job(store: JobStore, cfg: Config, job: Job, mailer=None) -> None:
    """Run a claimed job in a resource-capped ``spawn`` child, draining its
    progress into the store and recording done/failed. Per-job caps come from the
    job row (anonymous vs big), falling back to the anonymous config defaults for
    rows written before caps were set."""
    artifact_dir = cfg.artifacts_dir / job.id
    artifact_dir.mkdir(parents=True, exist_ok=True)
    progress_path = artifact_dir / "progress.json"
    store.mark_running(job.id)
    wall = job.wall_seconds or cfg.job_wall_seconds
    mem = job.mem_bytes or cfg.job_mem_bytes
    ctx = mp.get_context("spawn")
    q: "mp.Queue" = ctx.Queue()
    p = ctx.Process(target=_child,
                    args=(job.spec, str(artifact_dir), cfg.result_max_bytes,
                          wall, mem, q))
    p.start()
    # Parent-side wall-time backstop: RLIMIT_CPU caps CPU-seconds inside the
    # child, this caps wall-clock (with slack) so a wedged/sleeping child dies.
    deadline = time.monotonic() + wall + 5
    while p.is_alive() and time.monotonic() < deadline:
        p.join(1)
        _drain_progress(store, job.id, progress_path)
    if p.is_alive():
        p.terminate()
        p.join()
        store.mark_failed(job.id, "DepthLimitError: job exceeded the wall-time cap")
        _notify_if_big(store, cfg, job, "failed", mailer)
        return
    _drain_progress(store, job.id, progress_path)
    if q.empty():
        store.mark_failed(job.id, "RuntimeError: worker died without a result")
        _notify_if_big(store, cfg, job, "failed", mailer)
        return
    status, err = q.get()
    if status == "ok":
        store.mark_done(job.id, str(artifact_dir))
        _notify_if_big(store, cfg, job, "done", mailer)
    else:
        store.mark_failed(job.id, err)
        _notify_if_big(store, cfg, job, "failed", mailer)


def worker_tick(store: JobStore, cfg: Config, mailer=None) -> bool:
    """Claim one pending job and run it; return whether work was done."""
    job = store.claim_next()
    if job is None:
        return False
    run_one_job(store, cfg, job, mailer=mailer)
    return True
