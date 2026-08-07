"""Normalized bar (co)chain complexes over any Domain (spec §5, component 5:
the 'bar' oracle backend — exponential, small algebras only, always exact).

Plan 52: the (co)boundary matrices take an optional ``coefficients`` bimodule
``M`` (a :class:`quiverlab.hochschild.coefficients.Bimodule`). ``None`` is the
regular bimodule ``M = A`` and the code path is BYTE-IDENTICAL to the pre-Plan-52
one (a hard acceptance gate). With ``M`` the outer two faces use ``M``'s
left/right action (P1 of the plan): the M-slot ranges over ``M``'s basis
(``dim_M``), the interior bar factors stay ``A``'s multiplication, and
``cn = dim_M · (m-1)^n``. ``M = A`` (``Lact = Ract = A.T``) recovers this code
exactly."""
import itertools

from quiverlab.errors import DepthLimitError
from quiverlab.fields.linalg import rank
from quiverlab.hochschild.table import HHTable
from quiverlab.trace.events import ResolutionTerm
from quiverlab.trace.recorder import rankstep

_GUARD_HINT = ("the bar oracle is exponential; for a quiver-presented algebra the "
               "deeper engines serve this depth: pass engine=\"cs\" (Chouhy-Solotar, "
               "any admissible presentation over any exact field), or over GF(p) the "
               "fast/minimal/Bardzell stack — raise max_cells only if you know what "
               "you are doing")


def _abar_tuples(m, n):
    return list(itertools.product(range(1, m), repeat=n))


def _cochain_basis(m, n, dim_M=None):
    if dim_M is None:
        dim_M = m
    return [(s, J) for s in range(dim_M) for J in _abar_tuples(m, n)]


def _check_cells(rows, cols, max_cells, what):
    if rows * cols > max_cells:
        raise DepthLimitError(
            f"{what}: differential matrix would have {rows} x {cols} entries "
            f"(> max_cells = {max_cells})",
            hint=_GUARD_HINT,
        )


def _coeff_in_unit_basis(A, B, coefficients):
    """Re-express a coefficient bimodule built in ``A``'s own basis into the
    unit-adapted basis ``B = A.unit_adapted()`` the bar complex works in. ``None``
    (regular) stays ``None`` (byte-identical path); an already-unit-adapted ``A``
    needs no change (``B is A``)."""
    if coefficients is None:
        return None
    if A.is_unit_adapted:
        return coefficients
    return coefficients.change_of_basis(A._unit_adapting_change(), B)


def coboundary_matrix(A, n, max_cells, coefficients=None):
    """Matrix of d: C^n -> C^{n+1} for unit-adapted A. Rows: C^{n+1} basis, cols: C^n.
    ``coefficients`` (a Bimodule, or None=regular) enters at the outer two faces."""
    dom = A.domain
    m = A.dim
    dim_M = m if coefficients is None else coefficients.dim_M
    cols = _cochain_basis(m, n, dim_M)
    rows = _cochain_basis(m, n + 1, dim_M)
    _check_cells(len(rows), len(cols), max_cells, f"bar coboundary d^{n}")
    row_index = {b: i for i, b in enumerate(rows)}
    D = [[dom.zero()] * len(cols) for _ in range(len(rows))]

    if coefficients is None:
        def Lrow(j, s):
            return A.T[j][s]          # e_j · m_s   (m-vector)

        def Rrow(j, s):
            return A.T[s][j]          # m_s · e_j
    else:
        def Lrow(j, s):
            return coefficients.Lact[j][s]

        def Rrow(j, s):
            return coefficients.Ract[j][s]

    def bump(t, K, ci, val):
        if not dom.is_zero(val):
            r = row_index[(t, K)]
            D[r][ci] = dom.add(D[r][ci], val)

    for ci, (s, J) in enumerate(cols):
        for K in _abar_tuples(m, n + 1):
            # term 0: b_{K0} * f(K[1:])  [LEFT action on M]
            if K[1:] == J:
                vec = Lrow(K[0], s)
                for t in range(dim_M):
                    bump(t, K, ci, vec[t])
            # middle terms i = 1..n: (-1)^i f(..., proj(K[i-1] K[i]), ...)  [interior]
            for i in range(1, n + 1):
                pre, post = K[: i - 1], K[i + 1:]
                if pre == J[: i - 1] and post == J[i - 1 + 1:]:
                    x = J[i - 1]
                    coef = A.T[K[i - 1]][K[i]][x]  # Abar component only (x >= 1)
                    if not dom.is_zero(coef):
                        val = coef if i % 2 == 0 else dom.neg(coef)
                        bump(s, K, ci, val)
            # last term: (-1)^{n+1} f(K[:n]) * b_{K[n]}  [RIGHT action on M]
            if K[:n] == J:
                vec = Rrow(K[n], s)
                for t in range(dim_M):
                    val = vec[t] if (n + 1) % 2 == 0 else dom.neg(vec[t])
                    bump(t, K, ci, val)
    return D, len(cols), len(rows)


def hochschild_cohomology_dims(A, top, max_cells=4_000_000, trace=None, coefficients=None):
    B = A.unit_adapted()
    Mb = _coeff_in_unit_basis(A, B, coefficients)
    dom = B.domain
    m = B.dim
    dim_M = m if Mb is None else Mb.dim_M
    # dim HH^n = dim C^n - rank(d^n) - rank(d^{n-1}) (d^{-1} = 0); single pass carries
    # the previous rank so only one differential matrix is live at a time. When
    # trace is not None, record one ResolutionTerm + one RankStep per cochain degree.
    prev = 0
    dims = []
    for n in range(top + 1):
        D, ncols, nrows = coboundary_matrix(B, n, max_cells, coefficients=Mb)
        r = rank(D, dom) if nrows and ncols else 0
        cn = dim_M * (m - 1) ** n
        dims.append(cn - r - prev)
        if trace is not None:
            trace.append(ResolutionTerm(degree=n, n_generators=cn, collapsed_dim=cn))
            trace.append(rankstep(n, "cochain", D, nrows, ncols, r, dom))
        prev = r
    return HHTable(dims, "HH^", repr(A).splitlines()[0])


def boundary_matrix(A, n, max_cells, coefficients=None):
    """Matrix of b: C_n -> C_{n-1}, n >= 1, for unit-adapted A.
    C_n basis: (s, J) with s in 0..dim_M-1 (the M slot), J in Abar^{⊗n}."""
    dom = A.domain
    m = A.dim
    dim_M = m if coefficients is None else coefficients.dim_M
    cols = _cochain_basis(m, n, dim_M)        # same index shape: (s, J)
    rows = _cochain_basis(m, n - 1, dim_M)
    _check_cells(len(rows), len(cols), max_cells, f"bar boundary b_{n}")
    row_index = {b: i for i, b in enumerate(rows)}
    D = [[dom.zero()] * len(cols) for _ in range(len(rows))]

    if coefficients is None:
        def Lrow(j, s):
            return A.T[j][s]          # e_j · m_s

        def Rrow(j, s):
            return A.T[s][j]          # m_s · e_j
    else:
        def Lrow(j, s):
            return coefficients.Lact[j][s]

        def Rrow(j, s):
            return coefficients.Ract[j][s]

    def bump(t, K, ci, val):
        if not dom.is_zero(val):
            r = row_index[(t, K)]
            D[r][ci] = dom.add(D[r][ci], val)

    for ci, (s, J) in enumerate(cols):
        # term 0 (+): (m · b_{J0}) ⊗ J[1:]   [RIGHT action on M]
        vec = Rrow(J[0], s)
        for t in range(dim_M):
            bump(t, J[1:], ci, vec[t])
        # middle terms i = 1..n-1, sign (-1)^i: merge J[i-1], J[i], project to Abar
        for i in range(1, n):
            merged = A.T[J[i - 1]][J[i]]
            for x in range(1, m):
                coef = merged[x]
                if not dom.is_zero(coef):
                    val = coef if i % 2 == 0 else dom.neg(coef)
                    bump(s, J[: i - 1] + (x,) + J[i + 1:], ci, val)
        # last term, sign (-1)^n: (b_{J[n-1]} · m) ⊗ J[:n-1]   [LEFT action on M]
        vec = Lrow(J[n - 1], s)
        for t in range(dim_M):
            val = vec[t] if n % 2 == 0 else dom.neg(vec[t])
            bump(t, J[: n - 1], ci, val)
    return D, len(cols), len(rows)


def hochschild_homology_dims(A, top, max_cells=4_000_000, trace=None, coefficients=None):
    B = A.unit_adapted()
    Mb = _coeff_in_unit_basis(A, B, coefficients)
    dom = B.domain
    m = B.dim
    dim_M = m if Mb is None else Mb.dim_M
    ranks = [0]  # rank of b_n, with b_0 = 0
    for n in range(1, top + 2):
        D, ncols, nrows = boundary_matrix(B, n, max_cells, coefficients=Mb)
        r = rank(D, dom) if nrows and ncols else 0
        # RankStep for the degree-n boundary b_n: C_n -> C_{n-1} (chain side). Emitted
        # as each matrix is built, so only one differential is live at a time.
        if trace is not None:
            trace.append(rankstep(n, "chain", D, nrows, ncols, r, dom))
        ranks.append(r)
    dims = []
    for n in range(top + 1):
        cn = dim_M * (m - 1) ** n
        dims.append(cn - ranks[n] - ranks[n + 1])
        if trace is not None:
            trace.append(ResolutionTerm(degree=n, n_generators=cn, collapsed_dim=cn))
    return HHTable(dims, "HH_", repr(A).splitlines()[0])
