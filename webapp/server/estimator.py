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


def _max_degree(req: ComputeRequest) -> int:
    hi = 0
    for raw in req.compute:
        item = parse_compute_item(raw)
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


def sizing_dim(algebra_dim: int, req: ComputeRequest) -> int:
    """Effective dimension for tier classification. Module resolutions and Ext
    scale with the MODULE dimension, so a big module must size the job even over a
    small algebra -- otherwise it would be mis-classified as instant. Falls back to
    the algebra dimension when there is no explicit module, so every existing
    family/quiver request classifies exactly as before (Plan 26)."""
    return max(algebra_dim, _module_dim(req.module), _module_dim(req.ext_target))


# Heuristic throughput used to turn the op estimate into a human "minutes"
# figure for the warning UX (config-overridable would be trivial; a constant is
# fine for an order-of-magnitude hint).
_OPS_PER_MINUTE = 500_000_000


def _fits_big(ops: int, max_deg: int, cfg: Config) -> bool:
    return ops <= cfg.big_ops_threshold and max_deg <= cfg.big_max_degree


def classify(dim: int, req: ComputeRequest, cfg: Config) -> dict:
    """Full tier decision WITH the honest numbers the warning UX shows.
    Returns {"tier", "reason", "estimate": {"cells", "minutes"}}.
    `reason` is None unless tier == "reject", where it is "big_disabled" (fits
    big caps but SMTP is off) or "beyond_big_cap" (exceeds big caps)."""
    max_deg = _max_degree(req)
    ops = estimate_ops(dim, max_deg, req.algebra.field.kind)
    minutes = max(1, -(-ops // _OPS_PER_MINUTE))          # ceil division, ≥ 1
    est = {"cells": ops, "minutes": minutes}
    if ops <= cfg.instant_ops_threshold and max_deg <= cfg.instant_max_degree:
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
