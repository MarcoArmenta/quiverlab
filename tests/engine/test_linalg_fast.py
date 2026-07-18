"""Sparse exact F_p rank, cross-checked against the dense reference (PLAN item B3).

The dense `rank_mod_p` is the ground truth; `sparse_rank_mod_p` must return exactly the
same rank on every input (random sparse matrices over several primes, the real bar
differentials, and edge cases), and `rank_mod_p_auto` must agree with the reference
whichever path it dispatches to.
"""

import numpy as np
import pytest

from quiverlab.engine.linalg_fast import sparse_rank_mod_p, to_sparse_cols, rank_mod_p_auto

PRIMES = [2, 3, 5, 32003]


def test_edge_cases():
    p = 32003
    # empty
    assert sparse_rank_mod_p([], p, 0) == 0
    assert rank_mod_p_auto(np.zeros((0, 0), dtype=np.int64), p) == 0
    # zero matrix
    Z = np.zeros((5, 4), dtype=np.int64)
    cols, nr = to_sparse_cols(Z)
    assert sparse_rank_mod_p(cols, p, nr) == 0
    # identity
    I = np.eye(7, dtype=np.int64)
    cols, nr = to_sparse_cols(I)
    assert sparse_rank_mod_p(cols, p, nr) == 7
    # entries that are multiples of p collapse to 0
    M = np.array([[p, 2 * p], [0, p]], dtype=np.int64)
    cols, nr = to_sparse_cols(M)
    assert sparse_rank_mod_p(cols, p, nr) == 0


def test_large_sparse_known_rank():
    """A large, very sparse matrix where the dense array would be wasteful."""
    p = 101
    n = 1500
    M = np.zeros((n, n), dtype=np.int64)
    for i in range(n - 1):
        M[i, i] = 1
        M[i + 1, i] = 1                    # bidiagonal -> full column rank n-1 on cols 0..n-2
    cols, nrows = to_sparse_cols(M)
    assert sparse_rank_mod_p(cols, p, nrows) == n - 1


from quiverlab.engine.hh_engine import (
    rank_mod_p,
    cn_basis,
    differential_matrix,
    truncated_polynomial,
    two_gen_local,
)


@pytest.mark.parametrize("seed", range(25))
def test_sparse_matches_dense_random(seed):
    rng = np.random.default_rng(seed)
    r = int(rng.integers(1, 40))
    c = int(rng.integers(1, 40))
    p = int(rng.choice(PRIMES))
    M = np.zeros((r, c), dtype=np.int64)
    for _ in range(int(rng.integers(0, r * c // 2 + 1))):
        M[rng.integers(0, r), rng.integers(0, c)] = int(rng.integers(0, p))
    cols, nrows = to_sparse_cols(M)
    assert sparse_rank_mod_p(cols, p, nrows) == rank_mod_p(M, p)


@pytest.mark.parametrize("alg_key", ["k[x]/(x^3)", "qci3", "cyclic_nakayama(3,2)"])
def test_sparse_matches_dense_on_bar_differentials(alg_key):
    # enabled in Task 11 (coxeter2): the cyclic_nakayama case needs coxeter2;
    # gating the whole test keeps the bar-differential parity check together and
    # avoids building cyclic_nakayama in the parametrize at collection time.
    pytest.importorskip("quiverlab.engine.coxeter2")
    from quiverlab.engine.coxeter2 import cyclic_nakayama
    builders = {
        "k[x]/(x^3)": lambda: truncated_polynomial(3),
        "qci3": lambda: two_gen_local([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, -3], "qci3"),
        "cyclic_nakayama(3,2)": lambda: cyclic_nakayama(3, 2)[0],
    }
    alg = builders[alg_key]()
    for n in range(1, 4):
        bn = cn_basis(alg, n)
        idx = {g: i for i, g in enumerate(cn_basis(alg, n - 1))}
        M = differential_matrix(alg, n, bn, idx)
        cols, nrows = to_sparse_cols(M)
        for p in (2, 3, 32003):
            assert sparse_rank_mod_p(cols, p, nrows) == rank_mod_p(M, p)


def test_auto_dispatch_agrees_with_reference():
    p = 32003
    # small/dense -> dense path
    rng = np.random.default_rng(1)
    M = rng.integers(0, p, size=(40, 40), dtype=np.int64)
    assert rank_mod_p_auto(M, p) == rank_mod_p(M, p)
    # large + sparse -> sparse path, still equal to the reference
    big = np.zeros((1000, 1000), dtype=np.int64)
    for i in range(1000):
        big[i, i] = 1                      # rank-1000 identity diagonal
    big[0, 999] = 7                        # a stray off-diagonal entry
    assert big.size >= 200_000
    assert rank_mod_p_auto(big, p) == rank_mod_p(big, p) == 1000
