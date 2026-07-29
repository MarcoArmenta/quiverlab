"""HH^*/HH_* dimensions and representative (co)cycles from the Chouhy-Solotar
resolution, over ANY Domain (spec Plan-04 Task 7).

dim HH^n = dim C^n - rank(delta^n) - rank(delta^{n-1})   (coh; delta^{-1} = 0)
dim HH_n = dim C_n - rank(b_n)     - rank(b_{n+1})        (hom; b_0     = 0)

The binding CS gates (d^2 = 0 and Theorem 4.1's order condition) run first, so a
returned dimension is certified, never approximate. Admissibility is enforced at
the boundary: CS runs only on a certified-confluent reduction system with a finite
irreducible basis."""


def _require_admissible(rs):
    from quiverlab.errors import AdmissibilityError
    if not rs.is_confluent or not rs.irreducibles:
        raise AdmissibilityError("CS runs only on a certified-admissible reduction system",
                                 hint="Groebner completion did not certify confluence / a finite basis")


def _corners(res, n):
    """The degree-n generators' (source, target) vertices, one pair per generator.

    The CS term is ``P_n = (+)_{s in S_n} A e_{o(s)} (x)_k e_{t(s)} A``, so this is
    exactly the projective-bimodule decomposition of that term -- recorded so the
    report can NAME the summands instead of quoting a generator count (Marco
    2026-07-29). Display data only: no engine value depends on it."""
    return [(ch.o, ch.t) for ch in res.ss.S(n)]


def _emit_cohomology_trace(trace, res, dom, top, mats, ranks):
    """Record the CS COHOMOLOGY worked steps: one ResolutionTerm + one RankStep per
    cochain degree 0..top (Plan 34 fix -- the old ``_emit_resolution_trace`` recorded only
    the ResolutionTerms, so the renderers' ``derive_dims`` re-derived C_n - 0 - 0 = C_n,
    the cochain-SPACE dims, and mislabelled HH^ as HH_). Emitting the differential
    ``delta^n : C^n -> C^{n+1}`` (side="cochain") with its rank makes
    ``derive_dims``/``_dims_kind`` reconstruct HH^n = C_n - rank(delta^n) - rank(delta^{n-1})
    and the correct HH^ variance, for EVERY renderer/caller. ``mats[n]`` is
    ``res.matrix(n,"coh")`` (already built for the rank), ``ranks[n]`` its rank."""
    from quiverlab.resolutions_cs.trace import ResolutionTerm
    from quiverlab.trace.recorder import rankstep
    for n in range(top + 1):
        trace.append(ResolutionTerm(degree=n, n_generators=len(res.ss.S(n)),
                                    collapsed_dim=res.dim_C(n, "coh"),
                                    corners=_corners(res, n)))
        # delta^n : C^n -> C^{n+1}, so dim_C(n+1) rows x dim_C(n) cols (matches mats[n]).
        trace.append(rankstep(n, "cochain", mats[n],
                              res.dim_C(n + 1, "coh"), res.dim_C(n, "coh"),
                              ranks[n], dom))


def _emit_homology_trace(trace, res, dom, top, bmats, ranks):
    """Record the CS HOMOLOGY worked steps: one ResolutionTerm per chain degree 0..top,
    then one RankStep per boundary ``b_d : C_d -> C_{d-1}`` for d = 1..top+1 (side="chain",
    RankStep degree=d <-> b_d, mirroring ``hochschild.bar``). ``derive_dims`` then
    reconstructs HH_n = C_n - rank(b_n) - rank(b_{n+1}) (b_0 = 0) and the HH_ variance.
    ``bmats[i]`` is ``res.matrix(i+1,"hom") = b_{i+1}`` (already built), ``ranks[i]`` its
    rank."""
    from quiverlab.resolutions_cs.trace import ResolutionTerm
    from quiverlab.trace.recorder import rankstep
    for n in range(top + 1):
        trace.append(ResolutionTerm(degree=n, n_generators=len(res.ss.S(n)),
                                    collapsed_dim=res.dim_C(n, "hom"),
                                    corners=_corners(res, n)))
    for i in range(top + 1):
        d = i + 1                                # b_d : C_d -> C_{d-1}
        trace.append(rankstep(d, "chain", bmats[i],
                              res.dim_C(d - 1, "hom"), res.dim_C(d, "hom"),
                              ranks[i], dom))


def cs_cohomology_dims(A, top, max_cells=4_000_000, trace=None):
    from quiverlab.fields.linalg import rank
    from quiverlab.hochschild.table import HHTable
    from quiverlab.resolutions_cs.build import reduction_system_of
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    rs = reduction_system_of(A); _require_admissible(rs)
    res = ChouhySolotarResolution(A, rs, max_degree=top + 1, max_cells=max_cells)
    res.assert_dd_zero(upto=top + 1, side="coh"); res.assert_order_condition(upto=top + 1)
    dom = A.domain
    # Retain the differential matrices only when tracing (so trace=None stays memory-O(1)
    # -- one live matrix at a time); the dims are computed identically either way.
    mats = [res.matrix(n, "coh") for n in range(top + 1)] if trace is not None else None
    r = [rank(mats[n] if mats is not None else res.matrix(n, "coh"), dom)
         for n in range(top + 1)]
    dims = [res.dim_C(n, "coh") - r[n] - (r[n - 1] if n else 0) for n in range(top + 1)]
    if trace is not None:
        _emit_cohomology_trace(trace, res, dom, top, mats, r)
    return HHTable(dims, "HH^", repr(A).splitlines()[0], engine="Chouhy-Solotar")


def cs_homology_dims(A, top, max_cells=4_000_000, trace=None):
    from quiverlab.fields.linalg import rank
    from quiverlab.hochschild.table import HHTable
    from quiverlab.resolutions_cs.build import reduction_system_of
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    rs = reduction_system_of(A); _require_admissible(rs)
    res = ChouhySolotarResolution(A, rs, max_degree=top + 1, max_cells=max_cells)
    res.assert_dd_zero(upto=top + 1, side="hom"); res.assert_order_condition(upto=top + 1)
    dom = A.domain
    # rk[i] = rank(b_{i+1}) = rank(matrix(i+1, "hom")); b_0 = 0 (no map out of C_0).
    # Retain b_{i+1} only when tracing (trace=None stays memory-O(1)).
    bmats = [res.matrix(n + 1, "hom") for n in range(top + 1)] if trace is not None else None
    rk = [rank(bmats[n] if bmats is not None else res.matrix(n + 1, "hom"), dom)
          for n in range(top + 1)]
    dims = [res.dim_C(n, "hom") - (rk[n - 1] if n else 0) - rk[n] for n in range(top + 1)]
    if trace is not None:
        _emit_homology_trace(trace, res, dom, top, bmats, rk)
    return HHTable(dims, "HH_", repr(A).splitlines()[0], engine="Chouhy-Solotar")


# -- representative (co)cycles ---------------------------------------------------
def _columns(M):
    return [[row[c] for row in M] for c in range(len(M[0]))] if M and M[0] else []


def _reps_mod_image(cycles, image, dom):
    """A subset of `cycles` that is linearly independent modulo span(image)."""
    from quiverlab.fields.linalg import rank
    reps, base = [], list(image)
    base_rank = rank(base, dom) if base else 0
    for v in cycles:
        rr = rank(base + [v], dom)
        if rr > base_rank:
            reps.append(v)
            base = base + [v]
            base_rank = rr
    return reps


def cs_hh_basis(A, n, side, max_cells=4_000_000):
    """Representative (co)cycles of HH^n (side="coh") / HH_n (side="hom"): a basis
    of Z modulo the relevant image, each returned as a coordinate vector in C^n / C_n
    (the CS basis order of dim_C). Admissibility-gated; the (co)cycle space comes from
    fields.linalg.nullspace of the relevant differential."""
    from quiverlab.fields.linalg import nullspace
    from quiverlab.resolutions_cs.build import reduction_system_of
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    rs = reduction_system_of(A); _require_admissible(rs)
    res = ChouhySolotarResolution(A, rs, max_degree=n + 1, max_cells=max_cells)
    dom = A.domain
    if side == "coh":
        cycles = nullspace(res.matrix(n, "coh"), dom)          # Z^n = ker delta^n
        image = _columns(res.matrix(n - 1, "coh")) if n else []  # B^n = im delta^{n-1}
    elif side == "hom":
        if n == 0:
            d0 = res.dim_C(0, "hom")                            # Z_0 = C_0 (b_0 = 0)
            cycles = [[dom.one() if i == j else dom.zero() for j in range(d0)] for i in range(d0)]
        else:
            cycles = nullspace(res.matrix(n, "hom"), dom)      # Z_n = ker b_n
        image = _columns(res.matrix(n + 1, "hom"))             # B_n = im b_{n+1}
    else:
        raise ValueError(f"side must be 'coh' or 'hom', got {side!r}")
    return _reps_mod_image(cycles, image, dom)
