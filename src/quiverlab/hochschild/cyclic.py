"""Cyclic homology over any exact Domain (Plan 19): Connes' B on the
normalized bar complex + the (b, B) bicomplex — the generic mirror of
engine/cyclic.py (which stays authoritative over GF(p) via the fast rank).

Conventions (identical to the engine; parity is gated in
tests/invariants/test_cyclic_generic.py):

  * C_n = A (x) Abar^{(x)n} on the unit-adapted basis: basis (s, J) with
    s in 0..m-1, J in {1..m-1}^n (hochschild/bar.py's shapes; basis 0 = 1_A).
  * b = the normalized bar boundary (hochschild/bar.py::boundary_matrix).
  * B(a_0 (x) a_1 (x) ... (x) a_n) =
        sum_{i=0}^{n} (-1)^{n i} 1 (x) a_i (x) ... (x) a_n (x) a_0 (x) ... (x) a_{i-1}
    — insert 1_A in the A-slot, rotate all n+1 entries; any rotation that
    puts the unit (basis index 0) in a bar slot dies in Abar = A/k.1.
  * HC_n = dim Tot_n - rank D_n - rank D_{n+1},
    Tot_n = C_n (+) C_{n-2} (+) ..., D = b + B.
"""
from quiverlab.errors import DepthLimitError
from quiverlab.fields.linalg import rank
from quiverlab.hochschild.bar import _cochain_basis, boundary_matrix
from quiverlab.hochschild.table import HHTable

_GUARD_HINT = ("the (b, B) bicomplex on the bar basis is exponential; over "
               "GF(p) use the fast engine — raise max_cells only if you know "
               "what you are doing")


def connes_B_matrix(A, n, max_cells):
    """Matrix of B : C_n -> C_{n+1} over A.domain for unit-adapted A.
    Returns (list-of-rows, ncols, nrows); rows indexed by the C_{n+1} basis."""
    dom = A.domain
    m = A.dim
    cols = _cochain_basis(m, n)
    rows = _cochain_basis(m, n + 1)
    if len(rows) * len(cols) > max_cells:
        raise DepthLimitError(
            f"Connes B_{n}: matrix would have {len(rows)} x {len(cols)} "
            f"entries (> max_cells = {max_cells})", hint=_GUARD_HINT)
    row_index = {b: i for i, b in enumerate(rows)}
    M = [[dom.zero()] * len(cols) for _ in range(len(rows))]
    one, negone = dom.one(), dom.neg(dom.one())
    for ci, (s, J) in enumerate(cols):
        entries = (s,) + J                       # (a_0, r_1..r_n), length n+1
        for i in range(n + 1):
            rotated = entries[i:] + entries[:i]
            if 0 in rotated:                     # the unit died in a bar slot
                continue
            r = row_index[(0, rotated)]          # new A-slot = the unit
            sign = one if (n * i) % 2 == 0 else negone
            M[r][ci] = dom.add(M[r][ci], sign)
    return M, len(cols), len(rows)


def _tot_degrees(n):
    """Chain degrees in Tot_n = C_n (+) C_{n-2} (+) ... (descending)."""
    return list(range(n, -1, -2))


def cyclic_homology_dims(A, top, max_cells=4_000_000, with_reps=False):
    """HHTable of dim HC_0..HC_top over A.domain via the (b, B) bicomplex.

    Works for ANY unital algebra (no quiver needed): only the unit-adapted
    basis is required. Exponential in top like the bar oracle; max_cells
    guards every assembled matrix.

    Plan 35 wave 3b -- ``with_reps=True`` returns ``(table, raw)`` where ``raw``
    exposes the explicit HC representatives from the SAME assembled total
    differentials that produce the ranks: ``raw['reps'][n]`` is a list of coordinate
    columns (Domain) over the ordered ``Tot_n = C_n (+) C_{n-2} (+) ...`` basis (a
    basis of ker D_n modulo im D_{n+1}), ``raw['totmats'][n]`` is ``D_n`` row-major
    (``[]`` for the 0-row D_0), ``raw['col_dims'][k] = dim C_k``. Labelling +
    serialization lives in ``hochschild.cyclic_reps``."""
    B0 = A.unit_adapted()
    dom = B0.domain
    m = B0.dim
    maxdeg = top + 1
    dims = {k: m * (m - 1) ** k for k in range(maxdeg + 2)}
    bmats = {k: boundary_matrix(B0, k, max_cells)[0] for k in range(1, maxdeg + 1)}
    Bmats = {k: connes_B_matrix(B0, k, max_cells)[0] for k in range(0, maxdeg)}
    ranks = {}
    Dstore = {}
    for n in range(top + 2):
        src, tgt = _tot_degrees(n), _tot_degrees(n - 1)
        if not tgt:
            ranks[n] = 0
            Dstore[n] = []                       # 0-row total differential (n = 0)
            continue
        row_off, off = {}, 0
        for d in tgt:
            row_off[d] = off
            off += dims[d]
        nrows = off
        ncols = sum(dims[d] for d in src)
        if nrows * ncols > max_cells:
            raise DepthLimitError(
                f"cyclic total differential D_{n}: {nrows} x {ncols} entries "
                f"(> max_cells = {max_cells})", hint=_GUARD_HINT)
        D = [[dom.zero()] * ncols for _ in range(nrows)]
        c0 = 0
        for d in src:
            for r0_key, blk in ((d - 1, bmats.get(d) if d >= 1 else None),
                                (d + 1, Bmats.get(d))):
                if blk is None or r0_key not in row_off:
                    continue
                r0 = row_off[r0_key]
                for r, rowvec in enumerate(blk):
                    Dr = D[r0 + r]
                    for c, val in enumerate(rowvec):
                        if not dom.is_zero(val):
                            Dr[c0 + c] = dom.add(Dr[c0 + c], val)
            c0 += dims[d]
        ranks[n] = rank(D, dom) if nrows and ncols else 0
        Dstore[n] = D
    out = []
    for n in range(top + 1):
        tot = sum(dims[d] for d in _tot_degrees(n))
        out.append(tot - ranks[n] - ranks[n + 1])
    table = HHTable(out, "HC_", repr(A).splitlines()[0],
                    engine=f"bar (b,B) mixed complex over {dom.name}")
    if not with_reps:
        return table
    return table, _capture_generic_reps(dims, Dstore, top, dom)


def _capture_generic_reps(dims, Dstore, top, dom):
    """(``raw`` payload) -- the explicit HC representatives over ``dom`` from the
    stored total differentials ``Dstore``. Each degree's classes are a basis of
    ker D_n modulo im D_{n+1}, picked by the greedy independent-modulo filter over
    ``fields.linalg`` (im D_{n+1} lies in ker D_n by D.D = 0, so the pick returns
    exactly dim Tot_n - rank D_n - rank D_{n+1} = HC_n columns)."""
    from quiverlab.fields.linalg import nullspace
    from quiverlab.fields.linalg import rank as _rank

    def cols_to_matrix(cols):
        if not cols:
            return []
        L = len(cols[0])
        return [[cols[c][r] for c in range(len(cols))] for r in range(L)]

    reps, totmats = {}, {}
    for n in range(top + 1):
        Dn = Dstore.get(n) or []
        tot_n = sum(dims[d] for d in _tot_degrees(n))
        if not Dn:                               # 0-row D_n: the whole space is cycles
            cycles = [[dom.one() if i == j else dom.zero() for i in range(tot_n)]
                      for j in range(tot_n)]
        else:
            cycles = nullspace(Dn, dom)
        Dn1 = Dstore.get(n + 1) or []
        image = ([[Dn1[r][c] for r in range(len(Dn1))] for c in range(len(Dn1[0]))]
                 if Dn1 and Dn1[0] else [])
        chosen, base = [], list(image)
        r0 = _rank(cols_to_matrix(base), dom) if base else 0
        for v in cycles:
            rr = _rank(cols_to_matrix(base + [v]), dom)
            if rr > r0:
                chosen.append(v)
                base.append(v)
                r0 = rr
        reps[n] = chosen
        totmats[n] = Dn
    return {"reps": reps, "totmats": totmats,
            "col_dims": {k: v for k, v in dims.items()}}
