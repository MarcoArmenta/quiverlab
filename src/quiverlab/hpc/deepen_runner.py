"""The ONE sanctioned bridge from the wheel CLI to an engine internal (Plan 28).

CLAUDE.md keeps ``quiverlab.engine.*`` internal; the Plan-28 design doc
(``docs/plans/2026-07-25-plan-28-hpc-container.md``, "Long-job story") explicitly
sanctions routing big Hochschild-homology jobs through
``quiverlab.engine.deepen.deepen`` -- the checkpointed, resumable minimal-A^e
driver. That reach is isolated to this single module so the rest of
``quiverlab.hpc`` stays on the public surface. The import is LAZY (inside
``run_deepen``), so merely importing ``quiverlab.hpc`` does not pull the engine."""
from __future__ import annotations

_DEFAULT_PRIME = 32003


def _resolve_prime(A, hpc) -> int:
    """The field ``deepen`` computes HH_* over. When the user pinned ``hpc.prime``
    use it; otherwise use the algebra's own prime (a GF(p) request gets HH over
    F_p, matching ``A.hochschild_homology``), falling back to a generic prime."""
    if hpc.prime is not None:
        return hpc.prime
    p = getattr(getattr(A, "domain", None), "p", None)
    return int(p) if p is not None else _DEFAULT_PRIME


def run_deepen(A, top: int, hpc, log=None):
    """Compute HH_0..HH_top of ``A`` over ``F_{hpc.prime}`` through the
    checkpointed ``deepen`` driver, resuming from ``hpc.checkpoint_dir`` if a
    prior run left a checkpoint.

    Returns ``(dims, complete, reason)``:
      * ``dims``     -- HH_0..HH_top (int list; zero-padded to ``top+1`` when the
                        resolution terminated early);
      * ``complete`` -- True if the requested degree was reached (or the
                        resolution terminated -- a genuine finite answer);
      * ``reason``   -- ``deepen``'s ``stop_reason`` (``terminated``/``term``/
                        ``max_degree``/``time``/``memory``).

    ``complete=False`` (a clean time/memory checkpoint stop) tells the CLI to exit
    75 so ``sbatch`` requeues and a rerun resumes from the checkpoint."""
    # SANCTIONED engine-internal reach (Plan 28 long-job story); lazy so the
    # import boundary of quiverlab.hpc stays clean until this path is exercised.
    # ``deepen`` drives the minimal A^e resolution over the ENGINE algebra (the
    # same int64 GF(p) structure the fast HH path consumes, via ``to_engine`` on
    # the unit-adapted algebra), so a core ``Algebra`` is bridged here.
    from quiverlab.engine.deepen import deepen
    from quiverlab.engine.adapter import to_engine

    E = to_engine(A.unit_adapted())
    prime = _resolve_prime(A, hpc)
    out = deepen(E, hpc.checkpoint_dir, prime=prime, max_degree=top,
                 max_transient_bytes=hpc.max_mem_bytes,
                 time_limit_s=hpc.time_limit_s, log=log)
    hh = [int(x) for x in out["HH"]]
    reason = out["stop_reason"]
    complete = len(hh) >= top + 1 or reason in ("terminated", "term")
    if complete:
        dims = hh[:top + 1]
        if len(dims) < top + 1:                 # terminated early: HH_n = 0 above
            dims = dims + [0] * (top + 1 - len(dims))
    else:
        dims = hh                               # partial; CLI raises CheckpointStop
    return dims, complete, reason
