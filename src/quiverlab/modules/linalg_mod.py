"""Domain-generic matrix arithmetic for modules, built on fields.linalg.

Matrices are list[list[elt]] (rows x cols); a vector is list[elt] (a column). All
coefficient work goes through the Domain (dom.add/sub/neg/mul/inv/is_zero) or through
fields.linalg, so everything is exact over QQ / QQ(alpha) / GF(p) / GF(p^n)."""
from quiverlab.fields import linalg


def zeros(r, c, dom):
    z = dom.zero()
    return [[z for _ in range(c)] for _ in range(r)]


def identity(n, dom):
    o, z = dom.one(), dom.zero()
    return [[o if i == j else z for j in range(n)] for i in range(n)]


def transpose(M):
    if not M:
        return []
    return [list(col) for col in zip(*M)]


def matmul(A, B, dom):
    if not A or not B:
        return []
    p, q, r = len(A), len(B), len(B[0])
    out = zeros(p, r, dom)
    for i in range(p):
        Ai = A[i]
        for k in range(q):
            a = Ai[k]
            if dom.is_zero(a):
                continue
            Bk = B[k]
            oi = out[i]
            for j in range(r):
                oi[j] = dom.add(oi[j], dom.mul(a, Bk[j]))
    return out


def matvec(A, v, dom):
    out = []
    for row in A:
        s = dom.zero()
        for a, x in zip(row, v):
            if not dom.is_zero(a):
                s = dom.add(s, dom.mul(a, x))
        out.append(s)
    return out


def kron(A, B, dom):
    """Kronecker product (row-block layout)."""
    ra, ca = len(A), len(A[0]) if A else 0
    rb, cb = len(B), len(B[0]) if B else 0
    out = zeros(ra * rb, ca * cb, dom)
    for i in range(ra):
        for j in range(ca):
            a = A[i][j]
            if dom.is_zero(a):
                continue
            for k in range(rb):
                for l in range(cb):
                    out[i * rb + k][j * cb + l] = dom.mul(a, B[k][l])
    return out


def vstack(mats):
    out = []
    for M in mats:
        out.extend(M)
    return out


def col(M, j):
    return [row[j] for row in M]


def cols_to_matrix(cols):
    """Assemble a matrix whose columns are the given column-vectors."""
    if not cols:
        return []
    n = len(cols[0])
    return [[c[i] for c in cols] for i in range(n)]


def mat_rank(M, dom):
    return linalg.rank(M, dom) if M else 0


def kernel_columns(M, dom):
    """Basis of {x : M x = 0} as a list of column-vectors."""
    if not M:
        return []
    return linalg.nullspace(M, dom)


def column_space_pivots(M, dom):
    """Indices of a maximal independent set of COLUMNS of M. `linalg.rref(M)` returns
    the pivot COLUMN indices of M directly, and the pivot columns of the original matrix
    are a basis of its column space -- so return those (NOT rref of the transpose, which
    would give row-space indices; that bug selected the zero column of P_1 on kA_2 and
    broke rad/top)."""
    if not M:
        return []
    _, piv = linalg.rref(M, dom)
    return sorted(piv)


def solve_columns(B, V, dom):
    """Express every column of V in the span of the columns of B. Returns the list of
    coefficient column-vectors x with B x = V[:, j], or None if any column is not in
    the span."""
    out = []
    for j in range(len(V[0]) if V else 0):
        x = linalg.solve(B, col(V, j), dom)
        if x is None:
            return None
        out.append(x)
    return out


def independent_modulo(cands, sub_cols, dom):
    """Greedy rank-growth: return the indices of the columns in `cands` that are
    linearly independent modulo the span of `sub_cols`. Used to pick minimal
    (mod-radical) generators."""
    basis = [list(c) for c in sub_cols]
    base_rank = mat_rank(cols_to_matrix(basis), dom) if basis else 0
    chosen = []
    for idx, c in enumerate(cands):
        trial = basis + [list(c)]
        r = mat_rank(cols_to_matrix(trial), dom)
        if r > base_rank:
            basis = trial
            base_rank = r
            chosen.append(idx)
    return chosen
