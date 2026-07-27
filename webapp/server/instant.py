"""Hard wall-time net for the instant tier. The cost estimator (Task 5) is a
heuristic; this backs it by running the synchronous computation in a spawned
child and killing it if it exceeds ``cfg.instant_wall_seconds`` -- the caller
then falls through to queueing. This is the ~5s net mandated by the design."""
from __future__ import annotations

import collections
import multiprocessing as mp
import os
import queue
import shutil
import tempfile
import threading
import time
from pathlib import Path

from webapp.server.config import Config
from webapp.server.runner import run_spec, RunError
from webapp.server.schema import ComputeRequest


class InstantRateLimiter:
    """Per-IP sliding-window throttle for the INSTANT tier.

    The instant path never enqueues a job, so the store-backed queue limiter
    (``check_can_queue``) never sees instant traffic and cannot throttle a flood --
    yet every instant request spawns a fresh resource-capped interpreter (~1-2s CPU,
    two pinned math threads). Unbounded per-IP instant requests are therefore a
    CPU/spawn-starvation DoS even though the memory vector is already contained
    (ops classifier + RLIMIT_AS). This is the missing rate gate: at most
    ``max_per_window`` instant computes per ``window_seconds`` per IP hash.

    In-memory and per-process (each uvicorn worker enforces its own window -- a real
    mitigation, not a distributed quota) and thread-safe for the sync threadpool the
    instant endpoint runs in. State is owned by the app instance (created in
    ``create_app``), NOT module-global, so tests don't cross-contaminate. A
    non-positive ``max_per_window`` disables the limiter (records nothing, allows
    all) -- an explicit off switch."""

    def __init__(self, max_per_window: int, window_seconds: int,
                 clock=time.monotonic) -> None:
        self._max = int(max_per_window)
        self._window = float(window_seconds)
        self._clock = clock
        self._lock = threading.Lock()
        self._hits: dict[str, "collections.deque[float]"] = {}

    def allow(self, ip_hash: str) -> bool:
        """Record an instant attempt for ``ip_hash`` and return True if it is within
        the window budget, False if it must be throttled (429). Thread-safe."""
        if self._max <= 0:
            return True                        # disabled -> never throttle
        now = self._clock()
        cutoff = now - self._window
        with self._lock:
            dq = self._hits.get(ip_hash)
            if dq is None:
                dq = collections.deque()
                self._hits[ip_hash] = dq
            while dq and dq[0] <= cutoff:       # evict stamps older than the window
                dq.popleft()
            if len(dq) >= self._max:
                # Keep the (full) deque so the window keeps sliding; do not record
                # the throttled attempt (it never runs a compute).
                return False
            dq.append(now)
            self._maybe_gc(cutoff)
            return True

    def _maybe_gc(self, cutoff: float) -> None:
        """Opportunistically drop IP buckets whose stamps have all expired, so the
        map cannot grow without bound across many one-off client IPs. Cheap: only
        sweeps when the map is large, and only under the held lock."""
        if len(self._hits) <= 4096:
            return
        stale = [ip for ip, dq in self._hits.items()
                 if not dq or dq[-1] <= cutoff]
        for ip in stale:
            self._hits.pop(ip, None)


def _child(spec_dict: dict, artifact_dir: str, result_max_bytes: int,
           cpu_seconds: int, mem_bytes: int, q: "mp.Queue") -> None:
    # Pin threads before run_spec triggers the lazy numba/BLAS import (CLAUDE.md
    # load-bearing ordering), then apply the SAME resource caps the queued-tier
    # worker enforces (RLIMIT_CPU on both OSes, RLIMIT_AS on Linux) so a runaway
    # sync compute -- e.g. an oversized module driving a giant dense matrix --
    # cannot OOM the box before the parent's wall-time kill fires. Reuse the
    # worker's cap helper (ONE source of truth, same loud-degradation-on-refusal
    # behavior); the import is here (not module-top) so it runs inside the spawn
    # child, after the thread pinning and before the lazy numba/BLAS import.
    os.environ["NUMBA_NUM_THREADS"] = "2"
    os.environ["OMP_NUM_THREADS"] = "2"
    from webapp.worker.worker import _apply_caps
    _apply_caps(cpu_seconds, mem_bytes)
    # then chdir so any stray ./quiverlab_traces/ lands inside the throwaway
    # artifact dir rather than the repo root.
    os.chdir(artifact_dir)
    try:
        result = run_spec(ComputeRequest.model_validate(spec_dict), Path(artifact_dir),
                          result_max_bytes=result_max_bytes)
        q.put(("ok", result))
    except RunError as exc:
        q.put(("fail", {"error_type": exc.error_type, "message": exc.message}))
    except Exception as exc:  # pragma: no cover - defensive; genericised upstream
        q.put(("fail", {"error_type": type(exc).__name__, "message": str(exc)}))


def run_with_timeout(req: ComputeRequest, cfg: Config) -> dict | None:
    """Return the result dict, or None if the instant run exceeded the wall net.

    Raises ``RunError`` if the computation failed loudly (a real error, not a
    timeout). The per-request temp dir is unique so concurrent instant runs never
    collide, and it is removed unconditionally -- instant results are not
    retained (only queued jobs keep artifacts)."""
    ctx = mp.get_context("spawn")
    result_q: "mp.Queue" = ctx.Queue()
    cfg.artifacts_dir.mkdir(parents=True, exist_ok=True)
    artifact_dir = Path(tempfile.mkdtemp(prefix="instant-", dir=cfg.artifacts_dir))
    # Per-child resource caps. RLIMIT_CPU gets headroom over the wall net for the
    # TWO pinned math threads (NUMBA/OMP=2) plus slack, so it is a true backstop
    # that never fires before the parent wall-time kill on a legitimate instant
    # compute; the parent join below remains the primary net. RLIMIT_AS reuses the
    # anonymous job memory ceiling -- a sync compute must not outweigh a queued one.
    cpu_seconds = cfg.instant_wall_seconds * 2 + 5
    mem_bytes = cfg.job_mem_bytes
    try:
        p = ctx.Process(target=_child,
                        args=(req.model_dump(by_alias=True), str(artifact_dir),
                              cfg.result_max_bytes, cpu_seconds, mem_bytes, result_q))
        p.start()
        p.join(cfg.instant_wall_seconds)
        if p.is_alive():
            p.terminate()                  # SIGTERM: ask the child to exit
            p.join(5)
            if p.is_alive():               # ignored SIGTERM -> escalate to SIGKILL
                p.kill()
                p.join()
            return None                    # exceeded the net -> caller queues it
        try:
            status, payload = result_q.get_nowait()
        except queue.Empty:
            return None                    # child died without a verdict -> queue it
        if status == "ok":
            return payload
        raise RunError(payload["error_type"], payload["message"])
    finally:
        shutil.rmtree(artifact_dir, ignore_errors=True)
