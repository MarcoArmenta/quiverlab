"""Relative Hochschild (co)homology over the vertex subalgebra E = kQ_0 (Plan 52,
Task 7). The E-relative reduced bar complex (Cibils):

    C_n^E = M (x)_{E^e} \\bar A_E^{(x)_E n},   \\bar A_E = A / E,   E = kQ_0.

E = kQ_0 is a SEPARABLE k-algebra (any field, any characteristic), so each
Bar_n^E is a projective A^e-bimodule and this computes the ABSOLUTE HH_.(A, M) /
HH^.(A, M) -- a proven, free strong oracle: relative == absolute DEGREEWISE.

Single-vertex (E = k): identical to the absolute normalized bar (quiverlab's
hochschild.bar reduces modulo k.1 only). Multi-vertex: a genuinely smaller,
distinct complex (only E-composable radical tuples survive) that agrees only
degreewise, by separability.

Convention: works in A's OWN (path-type) basis -- no unit-adaptation. The bar
factors range over the RADICAL basis (A/E), each with a source/target vertex; the
coefficient M must be E-homogeneous (each basis vector in one e_v M e_w). General
subalgebra B is refused loudly by the public dispatch.
"""
from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import rank
from quiverlab.hochschild.table import HHTable


def _grading(A, M):
    """(idem, rad, src, tgt, lv, rv): path-type data plus, per M-basis index s, the
    unique left/right vertices lv[s], rv[s] with e_{lv} . m_s = m_s = m_s . e_{rv}.
    Requires E-homogeneity of M (raise loudly otherwise)."""
    from quiverlab.invariants.pathbasis import path_type_basis
    idem, rad, src, tgt = path_type_basis(A, "relative HH over kQ_0")
    lv, rv = {}, {}
    for s in range(M.dim_M):
        es = M._unit_vec(s)
        lset = [v for v in idem if M.left_apply(A._basis_vec(v), es) == es]
        rset = [w for w in idem if M.right_apply(A._basis_vec(w), es) == es]
        if len(lset) != 1 or len(rset) != 1:
            raise QuiverlabError(
                "relative HH over kQ_0 needs an E-homogeneous coefficient (each M "
                "basis vector must lie in a single corner e_v M e_w)",
                hint="the named bimodules on a path-type algebra are E-homogeneous")
        lv[s], rv[s] = lset[0], rset[0]
    return idem, rad, src, tgt, lv, rv


def _tuples(rad, src, tgt, start, end, n):
    """All composable radical tuples (r_1,...,r_n) with src(r_1)=start,
    tgt(r_i)=src(r_{i+1}), tgt(r_n)=end. n=0: the empty tuple iff start==end."""
    if n == 0:
        return [()] if start == end else []
    by_src = {}
    for r in rad:
        by_src.setdefault(src[r], []).append(r)
    out = []

    def rec(chain, cur):
        if len(chain) == n:
            if tgt[chain[-1]] == end:
                out.append(tuple(chain))
            return
        for r in by_src.get(cur, []):
            chain.append(r)
            rec(chain, tgt[r])
            chain.pop()

    rec([], start)
    return out


def _hom_basis(rad, src, tgt, lv, rv, dim_M, n):
    # M (x)_{E^e} ...: leading factor's RIGHT vertex meets r_1, its LEFT vertex meets r_n
    return [(s, J) for s in range(dim_M)
            for J in _tuples(rad, src, tgt, rv[s], lv[s], n)]


def _coh_basis(rad, src, tgt, lv, rv, dim_M, n):
    # Hom_{E^e}(..., M): the map J -> m_s needs m_s in e_{src(r_1)} M e_{tgt(r_n)}
    return [(s, J) for s in range(dim_M)
            for J in _tuples(rad, src, tgt, lv[s], rv[s], n)]


def _boundary(A, M, rad, gr, n):
    """b_n : C_n^E -> C_{n-1}^E (homology), mirroring the reduced-bar boundary on the
    E-composable radical basis with M's actions at the outer faces."""
    _idem, _rad, src, tgt, lv, rv = gr
    dom = A.domain
    cols = _hom_basis(rad, src, tgt, lv, rv, M.dim_M, n)
    rows = _hom_basis(rad, src, tgt, lv, rv, M.dim_M, n - 1)
    ridx = {b: i for i, b in enumerate(rows)}
    D = [[dom.zero()] * len(cols) for _ in range(len(rows))]

    def bump(key, ci, val):
        if key in ridx and not dom.is_zero(val):
            r = ridx[key]
            D[r][ci] = dom.add(D[r][ci], val)

    for ci, (s, J) in enumerate(cols):
        es = M._unit_vec(s)
        # term 0 (+): (m_s . r_{J0}) (x) J[1:]   [RIGHT action]
        mv = M.right_apply(A._basis_vec(J[0]), es)
        for c in range(M.dim_M):
            bump((c, J[1:]), ci, mv[c])
        # interior i=1..n-1, sign (-1)^i: merge J[i-1], J[i] (rad . rad subset rad)
        for i in range(1, n):
            prod = A.multiply(A._basis_vec(J[i - 1]), A._basis_vec(J[i]))
            for x in rad:
                cx = prod[x]
                if not dom.is_zero(cx):
                    val = cx if i % 2 == 0 else dom.neg(cx)
                    bump((s, J[: i - 1] + (x,) + J[i + 1:]), ci, val)
        # last term, sign (-1)^n: (r_{J[n-1]} . m_s) (x) J[:n-1]   [LEFT action]
        mv = M.left_apply(A._basis_vec(J[n - 1]), es)
        for c in range(M.dim_M):
            val = mv[c] if n % 2 == 0 else dom.neg(mv[c])
            bump((c, J[: n - 1]), ci, val)
    return D


def _coboundary(A, M, rad, gr, n):
    """delta^n : C^n_E -> C^{n+1}_E (cohomology), mirroring the reduced-bar
    coboundary on the E-composable radical basis."""
    _idem, _rad, src, tgt, lv, rv = gr
    dom = A.domain
    cols = _coh_basis(rad, src, tgt, lv, rv, M.dim_M, n)
    rows = _coh_basis(rad, src, tgt, lv, rv, M.dim_M, n + 1)
    ridx = {b: i for i, b in enumerate(rows)}
    D = [[dom.zero()] * len(cols) for _ in range(len(rows))]
    for ci, (s, J) in enumerate(cols):
        es = M._unit_vec(s)
        for ri, (t, K) in enumerate(rows):
            acc = dom.zero()
            # term 0: e_{K0} . f(K[1:]) with f(J)=m_s   [LEFT action, component t]
            if K[1:] == J:
                mv = M.left_apply(A._basis_vec(K[0]), es)
                acc = dom.add(acc, mv[t])
            # interior i=1..n, sign (-1)^i: f(..., merge(K[i-1],K[i]), ...)  (t == s)
            if t == s:
                for i in range(1, n + 1):
                    if K[: i - 1] == J[: i - 1] and K[i + 1:] == J[i:]:
                        x = J[i - 1]
                        prod = A.multiply(A._basis_vec(K[i - 1]), A._basis_vec(K[i]))
                        cx = prod[x]
                        if not dom.is_zero(cx):
                            acc = dom.add(acc, cx if i % 2 == 0 else dom.neg(cx))
            # last: f(K[:n]) . e_{K[n]}   [RIGHT action, sign (-1)^{n+1}, component t]
            if K[:n] == J:
                mv = M.right_apply(A._basis_vec(K[n]), es)
                v = mv[t]
                acc = dom.add(acc, v if (n + 1) % 2 == 0 else dom.neg(v))
            if not dom.is_zero(acc):
                D[ri][ci] = dom.add(D[ri][ci], acc)
    return D


def _regular_or(coefficients, A):
    from quiverlab.hochschild.coefficients import Bimodule
    return coefficients if coefficients is not None else Bimodule.regular(A)


def relative_homology_dims(A, top, coefficients=None, max_cells=4_000_000):
    M = _regular_or(coefficients, A)
    gr = _grading(A, M)
    rad = gr[1]
    dom = A.domain
    ranks = [0]
    for n in range(1, top + 2):
        D = _boundary(A, M, rad, gr, n)
        r = rank(D, dom) if D and D[0] else 0
        ranks.append(r)
    dims = []
    for n in range(top + 1):
        cn = len(_hom_basis(rad, gr[2], gr[3], gr[4], gr[5], M.dim_M, n))
        dims.append(cn - ranks[n] - ranks[n + 1])
    return HHTable(dims, "HH_", repr(A).splitlines()[0], engine="E-relative reduced bar")


def relative_cohomology_dims(A, top, coefficients=None, max_cells=4_000_000):
    M = _regular_or(coefficients, A)
    gr = _grading(A, M)
    rad = gr[1]
    dom = A.domain
    prev = 0
    dims = []
    for n in range(top + 1):
        D = _coboundary(A, M, rad, gr, n)
        r = rank(D, dom) if D and D[0] else 0
        cn = len(_coh_basis(rad, gr[2], gr[3], gr[4], gr[5], M.dim_M, n))
        dims.append(cn - r - prev)
        prev = r
    return HHTable(dims, "HH^", repr(A).splitlines()[0], engine="E-relative reduced bar")
