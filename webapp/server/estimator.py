"""Server-side cost estimate from library-provided data only (algebra
dimension, requested degree, field kind). Heuristic — a hard wall-time net in
the request path (Task 9) backs it up. All knobs are config-overridable."""
from __future__ import annotations

from webapp.server.config import Config
from webapp.server.schema import ComputeRequest, parse_compute_item

# Exact-CC arithmetic (sympy algebraic numbers) is far slower per op than the
# GF(p) kernel stack; this multiplier reflects that, not a precise timing.
_FIELD_MULT = {"GF": 1, "CC": 50}


def estimate_ops(dim: int, max_degree: int, field_kind: str) -> int:
    base = (dim ** 3) * (max_degree + 1)     # bar-complex rank cost, order of
    return base * _FIELD_MULT.get(field_kind, _FIELD_MULT["CC"])


# Bytes per matrix cell in the peak dense differential, by field kind. GF(p) rides
# the int64 kernel (8 bytes/cell); an exact CC entry is a sympy algebraic number,
# an order-of-magnitude heavier per cell -- 64 is an honest ESTIMATE, not a
# measurement.
_BYTES_PER_CELL = {"GF": 8, "CC": 64}


def estimate_bytes(dim: int, max_degree: int, field_kind: str) -> int:
    """Order-of-magnitude PEAK memory (bytes) for the exact computation -- an
    ESTIMATE, honest and integer-only (webapp/ is float-exempt, but exact int math
    is free here). Model: the heaviest step materialises a dense differential whose
    footprint scales with the bar-complex matrix -- ``~dim**2`` cells per degree,
    ``_BYTES_PER_CELL[field]`` bytes each. It sits one power of ``dim`` below the
    op model (:func:`estimate_ops`) -- a matrix, not a matmul. It is a guide for the
    memory-visibility UX, NOT the ``RLIMIT_AS`` the worker actually enforces."""
    cells = (dim ** 2) * (max_degree + 1)
    return cells * _BYTES_PER_CELL.get(field_kind, _BYTES_PER_CELL["CC"])


_UNITS = ("B", "KiB", "MiB", "GiB", "TiB", "PiB")


def human_bytes(n: int) -> str:
    """A human-readable byte count with binary units, INTEGER math only (no
    floats): the largest unit under which ``n`` is at least 1, with one decimal
    place computed by integer division. ``human_bytes(140000) == '136.7 KiB'``."""
    n = int(n)
    if n < 0:
        n = 0
    k, base = 0, 1
    while n >= base * 1024 and k < len(_UNITS) - 1:
        base *= 1024
        k += 1
    if k == 0:
        return f"{n} B"
    whole = n // base
    tenths = (n - whole * base) * 10 // base
    return f"{whole}.{tenths} {_UNITS[k]}"


def _max_degree(req: ComputeRequest) -> int:
    hi = 0
    for raw in req.compute:
        item = parse_compute_item(raw)
        # tau_tilting's `hi` is a PAIR BUDGET, not a homological degree (Plan 45), and
        # ar_quiver's `hi` is a MODULE BUDGET (wave 2) -- neither is a degree, so a budget
        # of 512 must not drive the degree-based tier classification. Both are excluded
        # here and sized on the algebra dimension (sizing_dim), like the HH kinds.
        #
        # KNOWN LIMITATION (inherited from the tau_tilting precedent, NOT introduced
        # here -- estimator redesign is out of scope, other owners): because the budget
        # is dropped and sizing_dim keys only off the algebra dimension, a
        # representation-INFINITE algebra of small dimension paired with a large
        # ar_quiver/tau_tilting budget can be mislabelled "instant" even though the knit
        # may run long before it hits the budget cap. The wall-clock/memory caps still
        # bound it once running; a budget-aware sizing heuristic is the open backlog fix.
        if item.kind in ("tau_tilting", "ar_quiver"):
            continue
        if item.hi is not None:
            hi = max(hi, item.hi)
    return hi


def _module_dim(mspec) -> int:
    """The declared total dimension of a module spec, WITHOUT building it. A
    builtin (simple/projective/injective) is bounded by the algebra itself, so it
    contributes nothing extra; an explicit module contributes the sum of its
    dimension vector (the cost driver for its resolutions / Ext)."""
    if mspec is None or mspec.builtin is not None or mspec.dims is None:
        return 0
    return sum(int(n) for n in mspec.dims.values())


def _algebra_b_dim(req: ComputeRequest) -> int:
    """The dimension of the SECOND algebra for ``derived_compare`` (wave 2), so a big
    B over a small A still sizes the job off the instant tier (derived_compare
    fingerprints BOTH algebras). Built DEFENSIVELY: an unbuildable B returns 0 here
    (its real error surfaces later in the runner as a clean 4xx), never raising inside
    the classifier -- app.py's ``classify(sizing_dim(...))`` is not wrapped."""
    spec = req.algebra_b
    if spec is None:
        return 0
    try:
        from quiverlab.hpc.spec import build_algebra
        return int(build_algebra(spec.model_dump()).dim)
    except Exception:
        return 0


def sizing_dim(algebra_dim: int, req: ComputeRequest) -> int:
    """Effective dimension for tier classification. Module resolutions, Ext and Tor
    scale with the MODULE dimension, so a big module -- INCLUDING the Tor second
    module ``tor_target`` -- must size the job even over a small algebra, otherwise
    it would be mis-classified as instant and let an oversized Tor target drive a
    multi-GB dense-matrix allocation in the sync tier. ``derived_compare`` likewise
    fingerprints a SECOND algebra ``algebra_b``, so its dimension sizes the job too.
    Falls back to the algebra dimension when there is no explicit module / second
    algebra, so every existing family/quiver request classifies exactly as before
    (Plan 26/30 + wave 2)."""
    return max(algebra_dim, _module_dim(req.module), _module_dim(req.ext_target),
               _module_dim(req.tor_target), _algebra_b_dim(req))


# Heuristic throughput used to turn the op estimate into a human "minutes"
# figure for the warning UX (config-overridable would be trivial; a constant is
# fine for an order-of-magnitude hint).
_OPS_PER_MINUTE = 500_000_000


def _fits_big(ops: int, max_deg: int, cfg: Config) -> bool:
    return ops <= cfg.big_ops_threshold and max_deg <= cfg.big_max_degree


def classify(dim: int, req: ComputeRequest, cfg: Config) -> dict:
    """Full tier decision WITH the honest numbers the warning UX shows.
    Returns {"tier", "reason", "estimate": {"cells", "minutes", "bytes",
    "mem_human"}}. `reason` is None unless tier == "reject" -- "big_disabled"
    (fits big caps but SMTP is off) or "beyond_big_cap" (exceeds big caps) -- or a
    would-be-instant request was upgraded to "queued" because it asks for the
    worked-steps report ("report_artifacts", see below). ``bytes`` is the integer
    memory ESTIMATE (see :func:`estimate_bytes`); ``mem_human`` is its binary-unit
    rendering."""
    max_deg = _max_degree(req)
    ops = estimate_ops(dim, max_deg, req.algebra.field.kind)
    minutes = max(1, -(-ops // _OPS_PER_MINUTE))          # ceil division, ≥ 1
    mem = estimate_bytes(dim, max_deg, req.algebra.field.kind)
    est = {"cells": ops, "minutes": minutes,
           "bytes": mem, "mem_human": human_bytes(mem)}
    if ops <= cfg.instant_ops_threshold and max_deg <= cfg.instant_max_degree:
        # The worked-steps report (``artifacts.pdf``) is written into the tier's
        # artifact dir, but the INSTANT tier discards that dir unconditionally AND
        # runs ``capture_reps=False`` (see webapp/server/instant.py) -- so a report
        # request served instantly would silently return NO report and no plain-HH
        # representatives. Only the queued tier keeps a persistent artifact dir, so
        # a report request that would classify instant is upgraded to queued. This
        # only ever downgrades instant->queued (a would-be-instant request always
        # fits the wider queued caps), never bypassing the big/reject logic below.
        # TikZ is NOT a trigger: it is written to the same dir and likewise lost on
        # instant, but the canvas GUI sets ``tikz: true`` on EVERY compute, so
        # gating on it would force every GUI request to queue -- and the diagram is
        # cheap and user-drawn, unlike the report.
        if req.artifacts.pdf:
            return {"tier": "queued", "reason": "report_artifacts", "estimate": est}
        return {"tier": "instant", "reason": None, "estimate": est}
    if ops <= cfg.queued_ops_threshold and max_deg <= cfg.queued_max_degree:
        return {"tier": "queued", "reason": None, "estimate": est}
    if _fits_big(ops, max_deg, cfg):
        if cfg.big_jobs_enabled:
            return {"tier": "big", "reason": None, "estimate": est}
        return {"tier": "reject", "reason": "big_disabled", "estimate": est}
    return {"tier": "reject", "reason": "beyond_big_cap", "estimate": est}


def decide_tier(dim: int, req: ComputeRequest, cfg: Config) -> str:
    return classify(dim, req, cfg)["tier"]
