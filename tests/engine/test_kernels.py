"""Unit equality tests: each numba kernel must match its pure-Python twin.

The pure-Python functions in resolutions_minimal/hh_engine are ground truth. If a
kernel disagrees, fix the kernel.
"""
import numpy as np
import pytest

from quiverlab.engine import _kernels
from quiverlab.engine._kernels import HAS_NUMBA, pow_mod, inv_mod

PRIMES = [2, 3, 5, 32003]
requires_numba = pytest.mark.skipif(not HAS_NUMBA, reason="numba not installed")


def test_has_numba_is_bool():
    assert isinstance(HAS_NUMBA, bool)
    assert isinstance(_kernels.USE_KERNELS, bool)


@pytest.mark.parametrize("p", PRIMES)
def test_inv_mod_matches_pow(p):
    for a in range(1, p if p < 50 else 50):
        inv = inv_mod(a, p)
        assert (a * inv) % p == 1
    # large-prime spot checks
    if p == 32003:
        for a in (2, 7, 12345, 32002):
            assert (a * inv_mod(a, p)) % p == 1


def test_pow_mod_matches_builtin():
    p = 32003
    for a, e in [(2, 10), (5, 0), (12345, 7), (32002, 31)]:
        assert pow_mod(a, e, p) == pow(a, e, p)


@requires_numba
@pytest.mark.parametrize("seed", range(20))
def test_rank_kernel_matches_reference(seed):
    pytest.importorskip("quiverlab.engine.hh_engine")  # ported in Task 3
    from quiverlab.engine.hh_engine import rank_mod_p
    from quiverlab.engine._kernels import rank_mod_p_kernel
    rng = np.random.default_rng(seed + 200)
    rows = int(rng.integers(1, 30))
    cols = int(rng.integers(1, 30))
    p = int(rng.choice(PRIMES))
    M = rng.integers(0, p, size=(rows, cols)).astype(np.int64)
    M[rng.random((rows, cols)) < 0.4] = 0
    assert rank_mod_p_kernel(np.ascontiguousarray(M), p) == rank_mod_p(M, p)


@requires_numba
def test_rank_kernel_edge_cases():
    pytest.importorskip("quiverlab.engine.hh_engine")  # ported in Task 3
    from quiverlab.engine.hh_engine import rank_mod_p
    from quiverlab.engine._kernels import rank_mod_p_kernel
    p = 32003
    assert rank_mod_p_kernel(np.eye(7, dtype=np.int64), p) == 7
    Z = np.zeros((5, 4), dtype=np.int64)
    assert rank_mod_p_kernel(Z, p) == 0
    M = np.array([[p, 2 * p], [0, p]], dtype=np.int64)
    assert rank_mod_p_kernel(M, p) == rank_mod_p(M, p) == 0


@requires_numba
@pytest.mark.parametrize("seed", range(20))
def test_nullspace_kernel_matches_pure(seed):
    pytest.importorskip("quiverlab.engine.resolutions_minimal")  # ported in Task 9
    from quiverlab.engine import _kernels
    from quiverlab.engine.resolutions_minimal import nullspace_mod_p
    from quiverlab.engine._kernels import nullspace_kernel
    rng = np.random.default_rng(seed)
    rows = int(rng.integers(1, 25))
    cols = int(rng.integers(1, 25))
    p = int(rng.choice(PRIMES))
    M = rng.integers(0, p, size=(rows, cols)).astype(np.int64)
    # zero out ~half to create real nullspaces
    mask = rng.random((rows, cols)) < 0.5
    M[mask] = 0
    # Oracle = the pure-Python path, forced on (USE_KERNELS is read at call time).
    # This stays a genuine kernel-vs-pure comparison even after Step 5 wires the
    # dispatcher — otherwise nullspace_mod_p would call the kernel and the test would
    # compare the kernel against itself.
    saved = _kernels.USE_KERNELS
    try:
        _kernels.USE_KERNELS = False
        py = nullspace_mod_p(M, p)                 # list of length-cols vectors (pure twin)
    finally:
        _kernels.USE_KERNELS = saved
    k = nullspace_kernel(np.ascontiguousarray(M), p)   # (nfree, cols)
    assert k.shape[0] == len(py)
    # every kernel row is in the nullspace and matches the pure-Python vector
    for i, v in enumerate(py):
        assert np.array_equal(k[i] % p, np.asarray(v) % p)
        assert np.array_equal((M @ (k[i] % p)) % p, np.zeros(rows, dtype=np.int64))


@requires_numba
@pytest.mark.parametrize("seed", range(20))
def test_independent_modulo_kernel_matches_pure(seed):
    pytest.importorskip("quiverlab.engine.resolutions_minimal")  # ported in Task 9
    from quiverlab.engine import _kernels
    from quiverlab.engine.resolutions_minimal import _independent_modulo
    from quiverlab.engine._kernels import independent_modulo_kernel
    rng = np.random.default_rng(seed + 100)
    L = int(rng.integers(2, 18))
    nspan = int(rng.integers(0, 12))
    ncand = int(rng.integers(1, 12))
    p = int(rng.choice(PRIMES))
    span = [rng.integers(0, p, size=L).astype(np.int64) for _ in range(nspan)]
    cand = [rng.integers(0, p, size=L).astype(np.int64) for _ in range(ncand)]
    # Oracle = the pure-Python path, forced on (USE_KERNELS is read at call time), so this
    # stays a genuine kernel-vs-pure comparison even after Step 5 wires the dispatcher.
    saved = _kernels.USE_KERNELS
    try:
        _kernels.USE_KERNELS = False
        chosen_py = _independent_modulo(span, cand, p)   # list of chosen vectors (pure twin)
    finally:
        _kernels.USE_KERNELS = saved
    span_arr = (np.array(span, dtype=np.int64).reshape(nspan, L) if nspan
                else np.zeros((0, L), dtype=np.int64))
    cand_arr = np.array(cand, dtype=np.int64).reshape(ncand, L)
    idx = independent_modulo_kernel(span_arr, cand_arr, p)
    # same NUMBER selected and same candidates (compare as the set of chosen rows)
    assert len(idx) == len(chosen_py)
    chosen_set = {tuple(int(x) % p for x in v) for v in chosen_py}
    kernel_set = {tuple(int(x) % p for x in cand_arr[i]) for i in idx}
    assert chosen_set == kernel_set


@requires_numba
def test_radK_kernel_matches_pure():
    pytest.importorskip("quiverlab.engine.resolutions_minimal")  # ported in Task 9
    pytest.importorskip("quiverlab.engine.scan3")                # ported in Task 4
    from quiverlab.engine.resolutions_minimal import (
        AeEngine, radical_basis, _rad_ae_columns, _rad_ab_pairs,
        nullspace_mod_p, _build_radK_py,
    )
    from quiverlab.engine.scan3 import quantum_ci
    from quiverlab.engine._kernels import radK_kernel
    A = quantum_ci(2)
    p = 32003
    eng = AeEngine(A, p)
    radAe = _rad_ae_columns(eng, radical_basis(A, p))
    rad_ab = _rad_ab_pairs(radAe, eng.m)
    # build a real degree-1 kernel to feed in
    aug = np.zeros((eng.m, eng.m2), dtype=np.int64)
    for a in range(eng.m):
        for b in range(eng.m):
            prod = eng.T[a, b, :]
            for c in np.nonzero(prod % p)[0]:
                aug[c, a * eng.m + b] = (aug[c, a * eng.m + b] + prod[c]) % p
    ker = nullspace_mod_p(aug, p)
    ker_arr = np.array(ker, dtype=np.int64).reshape(len(ker), eng.m2)
    pure = _build_radK_py(eng, ker, 1, p)                       # list of vectors
    kern = radK_kernel(eng.Lb, rad_ab, ker_arr, 1, eng.m2, p)   # (nnull*nrad, m2)
    assert kern.shape[0] == len(pure)
    for i, v in enumerate(pure):
        assert np.array_equal(kern[i] % p, v % p)


@requires_numba
def test_Dn_kernel_matches_pure():
    pytest.importorskip("quiverlab.engine.resolutions_minimal")  # ported in Task 9
    pytest.importorskip("quiverlab.engine.scan3")                # ported in Task 4
    from quiverlab.engine.resolutions_minimal import AeEngine, _build_Dn_py
    from quiverlab.engine.scan3 import quantum_ci
    from quiverlab.engine._kernels import Dn_kernel
    A = quantum_ci(2)
    p = 32003
    eng = AeEngine(A, p)
    cur_r = 1
    rng = np.random.default_rng(7)
    gens = [rng.integers(0, p, size=eng.m2 * cur_r).astype(np.int64) for _ in range(3)]
    pure = _build_Dn_py(eng, gens, cur_r, p)
    gens_arr = np.array(gens, dtype=np.int64).reshape(len(gens), eng.m2 * cur_r)
    kern = Dn_kernel(eng.Lb, gens_arr, cur_r, eng.m, eng.m2, p)
    assert np.array_equal(kern % p, pure % p)
