"""Hard wall-time net for the instant tier. The cost estimator (Task 5) is a
heuristic; this backs it by running the synchronous computation in a spawned
child and killing it if it exceeds ``cfg.instant_wall_seconds`` -- the caller
then falls through to queueing. This is the ~5s net mandated by the design."""
from __future__ import annotations

import multiprocessing as mp
import os
import queue
import shutil
import tempfile
from pathlib import Path

from webapp.server.config import Config
from webapp.server.runner import run_spec, RunError
from webapp.server.schema import ComputeRequest


def _child(spec_dict: dict, artifact_dir: str, result_max_bytes: int,
           q: "mp.Queue") -> None:
    # Pin threads before run_spec triggers the lazy numba/BLAS import (CLAUDE.md
    # load-bearing ordering), then chdir so any stray ./quiverlab_traces/ lands
    # inside the throwaway artifact dir rather than the repo root.
    os.environ["NUMBA_NUM_THREADS"] = "2"
    os.environ["OMP_NUM_THREADS"] = "2"
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
    try:
        p = ctx.Process(target=_child,
                        args=(req.model_dump(by_alias=True), str(artifact_dir),
                              cfg.result_max_bytes, result_q))
        p.start()
        p.join(cfg.instant_wall_seconds)
        if p.is_alive():
            p.terminate()
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
