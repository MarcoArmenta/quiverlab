# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
"""Sparse exact F_p rank for the (very sparse) bar / CS differentials (PLAN item B3).

The dense `hh_engine.rank_mod_p` is the reference oracle: correct, simple, and the
ground truth behind every backend.  But it materialises an `(rows, cols)` dense F_p
matrix, and the bar differential `dim C_n = m·(m−1)ⁿ` blows up — a deep differential is
a huge but extremely sparse matrix (each column has O(n) nonzeros), so the dense array
is mostly zeros and can exhaust memory (this is what forced the depth caps in the QPA /
cyclic-Nakayama oracles).

`sparse_rank_mod_p` does column-wise Gaussian elimination over F_p on a
dict-of-columns representation, never allocating the dense array.  It returns exactly
the same rank as `rank_mod_p` (cross-checked in tests/test_linalg_fast.py on random
sparse matrices and on the real bar differentials), and `rank_mod_p_auto` dispatches to
it automatically when the matrix is large and sparse, falling back to the dense
reference otherwise.

A GPU / FLINT / `galois` backend (PLAN item B3 (iii)) is not testable on this host (no
GPU, base deps are numpy + sympy only); the sparse path is the portable win, and the
dense reference is always kept for validation.
"""

import numpy as np

from quiverlab.engine import _kernels


def to_sparse_cols(M):
    """Dense integer matrix -> (list of {row: value} columns, n_rows)."""
    M = np.asarray(M)
    if M.size == 0:
        return [], (M.shape[0] if M.ndim == 2 else 0)
    nrows, ncols = M.shape
    cols = []
    for j in range(ncols):
        col = M[:, j]
        nz = np.nonzero(col)[0]
        cols.append({int(r): int(col[r]) for r in nz})
    return cols, nrows


def sparse_rank_mod_p(cols, p, nrows=None):
    """Rank over F_p of a sparse matrix given as a list of {row: value} columns.

    Column-wise elimination: each column is reduced against the already-chosen pivot
    columns, then (if still nonzero) contributes a new pivot.  Returns the rank.
    """
    pivots = {}        # pivot_row -> reduced pivot column (dict, leading 1 at pivot_row)
    rank = 0
    for raw in cols:
        c = {}
        for r, v in raw.items():
            v %= p
            if v:
                c[r] = v
        # eliminate every entry of c that sits on a known pivot row
        active = [r for r in c if r in pivots]
        while active:
            r = active.pop()
            f = c.get(r)
            if not f:
                continue
            pcol = pivots[r]
            for rr, vv in pcol.items():
                nv = (c.get(rr, 0) - f * vv) % p
                if nv:
                    c[rr] = nv
                else:
                    c.pop(rr, None)
            # elimination may have created entries on other pivot rows
            active = [rr for rr in c if rr in pivots]
        if c:
            r0 = min(c)                       # deterministic pivot choice
            inv = pow(c[r0], p - 2, p)
            pivots[r0] = {rr: (vv * inv) % p for rr, vv in c.items()}
            rank += 1
    return rank


def rank_mod_p_auto(M, p, sparse_min_size=200_000, max_density_permille=50):
    """Rank over F_p, dispatching to the sparse path for large sparse matrices.

    Falls back to the dense reference `rank_mod_p` for small or dense matrices, so the
    result is identical to the reference in every case.
    """
    from quiverlab.engine.hh_engine import rank_mod_p
    M = np.asarray(M)
    if M.size == 0:
        return 0
    size = M.size
    if size < sparse_min_size:
        return _kernels.rank_mod_p_kernel(np.ascontiguousarray(M.astype(np.int64)), p) \
            if _kernels.USE_KERNELS else rank_mod_p(M, p)
    nnz = int(np.count_nonzero(M % p))
    if 1000 * nnz > max_density_permille * size:
        return _kernels.rank_mod_p_kernel(np.ascontiguousarray(M.astype(np.int64)), p) \
            if _kernels.USE_KERNELS else rank_mod_p(M, p)
    cols, nrows = to_sparse_cols(M)
    return sparse_rank_mod_p(cols, p, nrows)
