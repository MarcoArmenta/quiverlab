"""Worker: claim one pending job, run it in a resource-capped child process,
record the outcome. The child caps threads (NUMBA/OMP/BLAS = 1) and sets
RLIMIT_CPU (Linux + macOS) and RLIMIT_AS (Linux) *before* it triggers the heavy
kernel code, so a runaway job cannot take down the VM; it streams per-item
checkpoints to a progress file the parent drains into the store. Multiprocessing
lives at the app layer only -- the library core stays single-process.

The child uses the ``spawn`` start method. Spawn re-executes the bootstrap: the
fresh interpreter re-imports the target's module chain (this module ->
``webapp.server.runner`` -> ``quiverlab`` -> numpy) *before* ``_child`` runs, so
those imports are NOT deferred. What IS deferred is numba: quiverlab imports it
lazily, inside the engine kernels, only when ``run_spec`` actually computes. The
thread-env pinning at the very top of ``_child`` therefore still precedes that
lazy numba import -- the load-bearing ordering per CLAUDE.md -- so
NUMBA_NUM_THREADS is in force before numba's threading layer initialises. (BLAS
reads OMP/OPENBLAS/MKL at first use; pinning them in ``_child`` covers that.)"""
from __future__ import annotations

import json
import logging
import multiprocessing as mp
import os
import queue
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


def _read_result(q: "mp.Queue"):
    """Non-blocking read of the child's ``(status, err)`` result. Returns None
    when the queue holds nothing -- the child exited (clean-without-a-result or
    signal-killed) without a verdict. Guarding ``get_nowait`` with ``queue.Empty``
    detects both a clean exit and a SIGKILL, unlike ``q.empty()`` which is only
    advisory in multiprocessing."""
    try:
        return q.get_nowait()
    except queue.Empty:
        return None


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
    except Exception as exc:  # a mail failure must not lose the result
        # Never log the address, even on failure: the exception TYPE name and the
        # job id only. No exc_info -- an SMTP exception's str carries the recipient.
        _log.warning("completion mail failed for job %s (type=%s)",
                     job.id, type(exc).__name__)
    finally:
        store.clear_email(job.id)


def run_one_job(store: JobStore, cfg: Config, job: Job, mailer=None) -> None:
    """Run a claimed job in a resource-capped ``spawn`` child, draining its
    progress into the store and recording done/failed. Per-job caps come from the
    job row (anonymous vs big), falling back to the anonymous config defaults for
    rows written before caps were set.

    ``claim_next`` has already flipped the row to ``running``; every normal exit
    path below reaches a terminal mark. Any *unexpected* parent-side failure
    (artifact-dir mkdir, Queue/Process construction or start, an error mid-drain)
    is caught, recorded as ``failed``, and re-raised -- the worker loop may crash
    loudly, but a claimed row must never be stranded in ``running`` forever."""
    marked = False
    try:
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
            p.join(10)
            if p.is_alive():                 # terminate ignored -- escalate to SIGKILL
                p.kill()
                p.join()
            store.mark_failed(job.id, "DepthLimitError: job exceeded the wall-time cap")
            marked = True
            _notify_if_big(store, cfg, job, "failed", mailer)
            return
        _drain_progress(store, job.id, progress_path)
        result = _read_result(q)
        if result is None:
            store.mark_failed(job.id, "RuntimeError: worker died without a result")
            marked = True
            _notify_if_big(store, cfg, job, "failed", mailer)
            return
        status, err = result
        if status == "ok":
            store.mark_done(job.id, str(artifact_dir))
            marked = True
            _notify_if_big(store, cfg, job, "done", mailer)
        else:
            store.mark_failed(job.id, err)
            marked = True
            _notify_if_big(store, cfg, job, "failed", mailer)
    except Exception as exc:
        # A parent-side failure before any terminal mark would strand the row in
        # `running`; record it as failed, then re-raise (loud, but never stuck).
        # `marked` guards against clobbering a verdict already written above.
        if not marked:
            store.mark_failed(job.id, f"worker error: {type(exc).__name__}: {exc}")
        raise


def worker_tick(store: JobStore, cfg: Config, mailer=None) -> bool:
    """Claim one pending job and run it; return whether work was done."""
    job = store.claim_next()
    if job is None:
        return False
    run_one_job(store, cfg, job, mailer=mailer)
    return True
