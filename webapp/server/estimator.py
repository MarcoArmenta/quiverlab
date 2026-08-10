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
        # SAME limitation for the Plan-59 `string_homological` kind: it KNITS the AR
        # quiver identically (measured >7 min churn at the default knit budget on the
        # small rep-INFINITE 2-Kronecker), yet it is a scalar kind with no budget in
        # `hi`, so it is sized purely on the algebra dimension (below) -- a small
        # rep-infinite algebra can likewise be mislabelled "instant" until the
        # wall-clock cap bounds the knit. (`toupie` is a small HH + graph scan, not
        # knit-heavy, so it is not in this caveat.)
        # Plan 69: the `barcode` kind on a COMMUTATIVE LADDER also KNITS the AR quiver
        # (knit-heavy), but -- UNLIKE ar_quiver/string_homological -- it is NOT left in
        # this caveat: `classify` catches it via `_barcode_knit_heavy` and upgrades it
        # instant->queued (reason="knit_heavy"). A plain A_n/zigzag barcode only
        # decomposes (module-sized) and stays instant-eligible.
        # wall_chamber's `hi` is a PAIR BUDGET too (Plan 63), sized on sizing_dim like
        # tau_tilting -- not a degree.
        # Plan 60: `tilted_check` is likewise knit-heavy with a MODULE budget in `hi`.
        # Plan 61: `recognizer_ladder` carries a MODULE BUDGET too (not a degree), so it
        # joins the skip tuple beside left_right_parts.
        # Plan 67: `silting` carries a RADIUS,BUDGET pair in (lo, hi) -- enumeration
        # bounds, not homological degrees -- so it joins the skip tuple (sized on A.dim).
        # Plan 65: `exceptional_sequences` carries an ENUMERATION BUDGET (not a degree).
        # Plan 64: `congruences` carries a PAIR BUDGET (like tau_tilting), not a degree.
        # Plan 70: `hh1_lie` carries a DIM BUDGET (caps A.dim for the Der solve), not a
        # homological degree -- sized on the algebra dimension (sizing_dim), so a big
        # algebra routes off the instant tier while the dim-220 Nakayama examples get an
        # honest budget refusal (the tau_tilting/products-omission precedent).
        # Plan 72: `split_extension` / `arrow_removal` carry a TOP-DEGREE budget in hi;
        # they are sized on the algebra dim (split_extension on 2*dim via sizing_dim's
        # extension-awareness below), NOT tiered by hi-as-degree, so they skip too.
        # Plan 74: `skew_group_hh` carries a TOP-DEGREE budget in hi; the smash A|xG is
        # dim |G|*dim A, so it is sized on the algebra dim (sizing_dim), not tiered by
        # hi-as-degree -- it skips too (the split_extension/arrow_removal precedent).
        # Plan 66: `tau_cluster` carries a PAIR BUDGET too (the exchange graph + the wide
        # poset + a reduction/sub-g-fan per object) -- knit-heavier than `congruences`, sized
        # on sizing_dim (A.dim), not a degree.
        if item.kind in ("tau_tilting", "wall_chamber", "ar_quiver",
                         "left_right_parts", "tilted_check", "recognizer_ladder",
                         "silting", "exceptional_sequences", "congruences",
                         "hh1_lie", "split_extension", "arrow_removal",
                         "skew_group_hh", "tau_cluster"):
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


def _coefficient_dim(req) -> int:
    """The declared dimension of an EXPLICIT coefficient bimodule (Plan 52). The bar
    cochain basis is dim_M*(m-1)^n, so a big explicit coefficient over a SMALL algebra
    drives a cost quadratic in dim_M and must size the job off the instant tier. A
    builtin coefficient is bounded by dim A (already in the max), so it adds nothing."""
    spec = getattr(req, "coefficients", None)
    if spec is None or spec.builtin is not None or spec.dim is None:
        return 0
    return int(spec.dim)


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


def _extension_dim(req: ComputeRequest, algebra_dim: int) -> int:
    """Plan 72 (DD-G1): the ``split_extension`` LES runs the CS Hom-complex over
    ``L = T(B)`` of dimension ``2*dim B``, so a ``split_extension`` compute item
    sizes the job on ``2*algebra_dim`` (it routes off instant like an oversized
    family). ``arrow_removal`` runs HH of ``A`` and ``B`` (both <= dim A), so it
    sizes on ``algebra_dim`` -- no extra term."""
    for raw in req.compute:
        try:
            if parse_compute_item(raw).kind == "split_extension":
                return 2 * algebra_dim
        except Exception:
            continue
    return 0


def sizing_dim(algebra_dim: int, req: ComputeRequest) -> int:
    """Effective dimension for tier classification. Module resolutions, Ext and Tor
    scale with the MODULE dimension, so a big module -- INCLUDING the Tor second
    module ``tor_target`` -- must size the job even over a small algebra, otherwise
    it would be mis-classified as instant and let an oversized Tor target drive a
    multi-GB dense-matrix allocation in the sync tier. ``derived_compare`` likewise
    fingerprints a SECOND algebra ``algebra_b``, so its dimension sizes the job too.
    ``split_extension`` (Plan 72) runs over the trivial extension ``2*dim B``, so it
    sizes on ``2*algebra_dim``. Falls back to the algebra dimension when there is no
    explicit module / second algebra, so every existing family/quiver request
    classifies exactly as before (Plan 26/30 + wave 2)."""
    return max(algebra_dim, _module_dim(req.module), _module_dim(req.ext_target),
               _module_dim(req.tor_target), _coefficient_dim(req), _algebra_b_dim(req),
               _extension_dim(req, algebra_dim))


# Heuristic throughput used to turn the op estimate into a human "minutes"
# figure for the warning UX (config-overridable would be trivial; a constant is
# fine for an order-of-magnitude hint).
_OPS_PER_MINUTE = 500_000_000


def _fits_big(ops: int, max_deg: int, cfg: Config) -> bool:
    return ops <= cfg.big_ops_threshold and max_deg <= cfg.big_max_degree


def _barcode_knit_heavy(req: ComputeRequest) -> bool:
    """True iff the request asks for ``barcode`` AND its algebra is a commutative ladder
    (Plan 69 / H2). A CL barcode runs a FULL AR knit (minutes even for a small algebra),
    so it must not be served on the instant tier -- where artifacts are discarded and
    ``capture_reps=False`` -- and is upgraded instant->queued in :func:`classify`. An
    A_n/zigzag barcode only ``decompose``s (module-sized) -> ``False``, so it stays
    instant-eligible (byte-identical classification to before). Built DEFENSIVELY -- an
    unbuildable algebra returns ``False`` (its real error surfaces later as a clean 4xx),
    because ``app.py``'s ``classify(...)`` is not wrapped."""
    if not any(parse_compute_item(r).kind == "barcode" for r in req.compute):
        return False
    try:
        from quiverlab.hpc.spec import build_algebra
        A = build_algebra(req.algebra.model_dump())
        from quiverlab.families.commutative_ladder import is_commutative_ladder
        return bool(is_commutative_ladder(A)[0])
    except Exception:
        return False


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
        # Plan 69 (H2): a `barcode` request on a COMMUTATIVE LADDER runs a full AR knit
        # (knit-heavy -- the same cost class as ar_quiver/string_homological), but it is
        # a scalar module kind with no `hi` budget, so sizing_dim would size it purely on
        # the small CL algebra dim and mislabel it instant. Unlike ar_quiver/
        # string_homological (whose knit-heaviness the KNOWN LIMITATION above cannot
        # catch), the CL barcode IS caught here and upgraded instant->queued; a plain
        # A_n/zigzag barcode only decomposes (module-sized) and stays instant-eligible.
        if _barcode_knit_heavy(req):
            return {"tier": "queued", "reason": "knit_heavy", "estimate": est}
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
