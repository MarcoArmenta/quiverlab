"""Item 2 (TODO P1): direct unit tests for the exact F_p linear-algebra primitives.

These routines (`coxeter.rref_mod_p`, `nullspace_mod_p`, `colspace_basis_mod_p`,
`solve_mod_p`, `quotient_induced`, and `hh_engine.rank_mod_p`) are the arithmetic
foundation under every dimension and every induced-action matrix in the project. A
bug here corrupts results SILENTLY -- no exception, just wrong ranks -- so they are
tested directly here on hand-built matrices, with positive, edge, and negative cases,
rather than only transitively through homology dims.
"""
import numpy as np
import pytest
# hanlab sys.path shim dropped: quiverlab uses absolute package imports.
from quiverlab.engine.hh_engine import rank_mod_p
from quiverlab.engine.coxeter import (
    rref_mod_p, nullspace_mod_p, colspace_basis_mod_p, solve_mod_p, quotient_induced,
)

P = 32003


# ----------------------------------------------------------------------
# rank_mod_p / rref_mod_p
# ----------------------------------------------------------------------
def test_rank_full_and_deficient():
    assert rank_mod_p(np.eye(4, dtype=np.int64), P) == 4
    # rank-1 matrix (all rows proportional)
    assert rank_mod_p(np.array([[1, 2, 3], [2, 4, 6], [3, 6, 9]], dtype=np.int64), P) == 1
    assert rank_mod_p(np.zeros((3, 3), dtype=np.int64), P) == 0


def test_rank_empty_matrix():
    assert rank_mod_p(np.zeros((0, 0), dtype=np.int64), P) == 0
    assert rank_mod_p(np.zeros((3, 0), dtype=np.int64), P) == 0


def test_rank_is_characteristic_dependent():
    # [[1,1],[1,1]] has rank 1 over any field; [[2,0],[0,2]] has rank 2 unless p|2.
    assert rank_mod_p(np.array([[2, 0], [0, 2]], dtype=np.int64), P) == 2
    assert rank_mod_p(np.array([[2, 0], [0, 2]], dtype=np.int64), 2) == 0   # collapses mod 2
    # entries equal to p (== 0 mod p) and entries near p both reduce correctly
    assert rank_mod_p(np.array([[P, 0], [0, P]], dtype=np.int64), P) == 0
    assert rank_mod_p(np.array([[P - 1, 0], [0, P - 1]], dtype=np.int64), P) == 2


def test_rank_plus_nullity_equals_cols():
    A = np.array([[1, 2, 3, 4], [2, 4, 6, 8], [1, 0, 1, 0]], dtype=np.int64)
    r = rank_mod_p(A, P)
    N = nullspace_mod_p(A, P)
    assert r + N.shape[1] == A.shape[1]


def test_rref_is_reduced_and_pivot_count_is_rank():
    A = np.array([[0, 2, 4], [1, 1, 1], [2, 0, -2]], dtype=np.int64)
    R, piv = rref_mod_p(A, P)
    assert len(piv) == rank_mod_p(A, P)
    # each pivot column is a unit vector in R (1 in its pivot row, 0 elsewhere)
    for i, c in enumerate(piv):
        col = R[:, c] % P
        expected = np.zeros(R.shape[0], dtype=np.int64); expected[i] = 1
        assert np.array_equal(col, expected)


# ----------------------------------------------------------------------
# nullspace_mod_p
# ----------------------------------------------------------------------
def test_nullspace_annihilates():
    A = np.array([[1, 1, 0], [0, 0, 1]], dtype=np.int64)   # rank 2, cols 3 -> nullity 1
    N = nullspace_mod_p(A, P)
    assert N.shape[1] == 1
    assert np.array_equal((A @ N) % P, np.zeros((A.shape[0], N.shape[1]), dtype=np.int64))


def test_nullspace_trivial_for_full_column_rank():
    A = np.array([[1, 0], [0, 1], [1, 1]], dtype=np.int64)  # full column rank
    N = nullspace_mod_p(A, P)
    assert N.shape == (2, 0)


def test_nullspace_full_for_zero_matrix():
    A = np.zeros((2, 3), dtype=np.int64)
    N = nullspace_mod_p(A, P)
    assert N.shape[1] == 3
    assert rank_mod_p(N, P) == 3


# ----------------------------------------------------------------------
# colspace_basis_mod_p
# ----------------------------------------------------------------------
def test_colspace_basis_independent_and_spans():
    A = np.array([[1, 2, 1], [2, 4, 0], [0, 0, 1]], dtype=np.int64)  # col2 = 2*col1
    B = colspace_basis_mod_p(A, P)
    assert B.shape[1] == rank_mod_p(A, P) == 2
    assert rank_mod_p(B, P) == B.shape[1]                       # independent
    # span unchanged: rank[B] == rank[B | A]
    assert rank_mod_p(np.concatenate([B, A], axis=1), P) == B.shape[1]


def test_colspace_basis_zero_matrix():
    B = colspace_basis_mod_p(np.zeros((3, 2), dtype=np.int64), P)
    assert B.shape == (3, 0)


# ----------------------------------------------------------------------
# solve_mod_p
# ----------------------------------------------------------------------
def test_solve_roundtrip_consistent():
    A = np.array([[1, 0], [0, 1], [1, 1]], dtype=np.int64)  # full column rank
    x_true = np.array([2, 3], dtype=np.int64)
    y = (A @ x_true) % P
    x = solve_mod_p(A, y, P)
    assert x is not None
    assert np.array_equal((A @ x) % P, y % P)


def test_solve_returns_none_when_inconsistent():
    A = np.array([[1, 0], [0, 1], [1, 1]], dtype=np.int64)  # col3 row forces y3 = y1+y2
    y = np.array([1, 0, 0], dtype=np.int64)                  # 0 != 1 -> not in colspan
    assert solve_mod_p(A, y, P) is None


# ----------------------------------------------------------------------
# quotient_induced: endomorphism on colspan(B_ker)/colspan(B_im)
# ----------------------------------------------------------------------
def _e(i, n):
    v = np.zeros(n, dtype=np.int64); v[i] = 1
    return v


def test_quotient_induced_known_matrix_and_dim():
    # ambient F_p^3; ker = span(e0, e1), im = span(e0); quotient ~ span(e1), dim 1.
    B_ker = np.stack([_e(0, 3), _e(1, 3)], axis=1)
    B_im = _e(0, 3).reshape(3, 1)
    # Sigma = diag(1, 2, 7): preserves span(e0,e1); on the quotient e1 |-> 2 e1.
    Sigma = np.diag([1, 2, 7]).astype(np.int64)
    M, d = quotient_induced(Sigma, B_im, B_ker, P)
    assert d == 1
    assert np.array_equal(M % P, np.array([[2]], dtype=np.int64))


def test_quotient_induced_action_modulo_image_is_identity():
    # Sigma: e1 |-> e1 + e0 (= e1 modulo the image span(e0)); acts as identity on quotient.
    B_ker = np.stack([_e(0, 3), _e(1, 3)], axis=1)
    B_im = _e(0, 3).reshape(3, 1)
    Sigma = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]], dtype=np.int64)
    M, d = quotient_induced(Sigma, B_im, B_ker, P)
    assert d == 1
    assert np.array_equal(M % P, np.array([[1]], dtype=np.int64))


def test_quotient_induced_raises_when_not_preserved():
    # Sigma: e1 |-> e2, which leaves colspan(B_ker)=span(e0,e1): must raise (not a chain map).
    B_ker = np.stack([_e(0, 3), _e(1, 3)], axis=1)
    B_im = _e(0, 3).reshape(3, 1)
    Sigma = np.array([[1, 0, 0], [0, 0, 0], [0, 1, 0]], dtype=np.int64)
    with pytest.raises(RuntimeError):
        quotient_induced(Sigma, B_im, B_ker, P)


def test_quotient_induced_empty_kernel():
    B_ker = np.zeros((3, 0), dtype=np.int64)
    B_im = np.zeros((3, 0), dtype=np.int64)
    M, d = quotient_induced(np.eye(3, dtype=np.int64), B_im, B_ker, P)
    assert d == 0 and M.shape == (0, 0)
