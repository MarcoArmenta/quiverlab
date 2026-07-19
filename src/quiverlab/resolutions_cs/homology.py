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


def cs_cohomology_dims(A, top, max_cells=4_000_000):
    from quiverlab.fields.linalg import rank
    from quiverlab.hochschild.table import HHTable
    from quiverlab.resolutions_cs.build import reduction_system_of
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    rs = reduction_system_of(A); _require_admissible(rs)
    res = ChouhySolotarResolution(A, rs, max_degree=top + 1, max_cells=max_cells)
    res.assert_dd_zero(upto=top + 1, side="coh"); res.assert_order_condition(upto=top + 1)
    dom = A.domain
    r = [rank(res.matrix(n, "coh"), dom) for n in range(top + 1)]
    dims = [res.dim_C(n, "coh") - r[n] - (r[n - 1] if n else 0) for n in range(top + 1)]
    return HHTable(dims, "HH^", repr(A).splitlines()[0], engine="Chouhy-Solotar")


def cs_homology_dims(A, top, max_cells=4_000_000):
    from quiverlab.fields.linalg import rank
    from quiverlab.hochschild.table import HHTable
    from quiverlab.resolutions_cs.build import reduction_system_of
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    rs = reduction_system_of(A); _require_admissible(rs)
    res = ChouhySolotarResolution(A, rs, max_degree=top + 1, max_cells=max_cells)
    res.assert_dd_zero(upto=top + 1, side="hom"); res.assert_order_condition(upto=top + 1)
    dom = A.domain
    # rk[n] = rank(b_{n+1}) = rank(matrix(n+1, "hom")); b_0 = 0 (no map out of C_0).
    rk = [rank(res.matrix(n + 1, "hom"), dom) for n in range(top + 1)]
    dims = [res.dim_C(n, "hom") - (rk[n - 1] if n else 0) - rk[n] for n in range(top + 1)]
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
