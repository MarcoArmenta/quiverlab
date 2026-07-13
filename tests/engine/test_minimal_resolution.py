"""Pin the minimal Aᵉ-resolution engine HH_* against the normalized-bar oracle.

Tor is resolution-independent, so minimal_homology_dims must equal
hochschild_homology_dims (the bar complex) on every algebra, in every degree the
bar complex can reach. This is the safety net for the numba rewrite: it runs on the
pure-Python engine today and must stay green after every kernel lands.
"""
import pytest

from quiverlab.engine.hh_engine import (
    truncated_polynomial,
    two_gen_local,
    hochschild_homology_dims,
)
from quiverlab.engine.scan3 import quantum_ci
from quiverlab.engine.resolutions_minimal import minimal_homology_dims
from quiverlab.engine.resolutions_minimal import minimal_homology_dims as _mhd
from quiverlab.engine import _kernels

# not ported: _open_33_0() -- k<x,y>/(x^3 - y^2, y^3, yx + xy), a genuine open-zone
# algebra -- and its two consumers (test_minimal_matches_bar_oracle_open, and the
# open-zone second half of test_kernel_path_equals_pure_path). They build the algebra
# via reduction_algebra.{algebra_from_reduction_system, make_reduction_system}, which
# is EXCLUDED from Plan 02 (deferred to Plans 04/06). The closed-form validation
# algebras below still pin minimal_homology_dims against the bar oracle across all primes.

PRIMES = (32003, 2, 3, 5)


# Validation algebras, at depths the bar complex reaches without MemoryError.
_VALIDATION = [
    (truncated_polynomial(2), 6),
    (truncated_polynomial(3), 5),
    (truncated_polynomial(4), 4),
    (quantum_ci(2), 5),
    (quantum_ci(3), 5),
    (two_gen_local([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], "k[x,y]/(x^2,y^2)"), 5),
]


@pytest.mark.parametrize("alg,N", _VALIDATION, ids=[a.name for a, _ in _VALIDATION])
def test_minimal_matches_bar_oracle(alg, N):
    mh = minimal_homology_dims(alg, N, primes=PRIMES)
    bh = hochschild_homology_dims(alg, N, primes=PRIMES)
    for p in PRIMES:
        assert mh[p] == bh[p][:len(mh[p])], f"{alg.name} p={p}"


@pytest.mark.skipif(not _kernels.HAS_NUMBA, reason="numba not installed")
def test_kernel_path_equals_pure_path():
    A = quantum_ci(2)
    save = _kernels.USE_KERNELS
    try:
        _kernels.USE_KERNELS = False
        pure = _mhd(A, 5, primes=PRIMES)
        _kernels.USE_KERNELS = True
        fast = _mhd(A, 5, primes=PRIMES)
    finally:
        _kernels.USE_KERNELS = save
    assert pure == fast
    # (open-zone second half not ported: builds _open_33_0 via reduction_algebra, excluded)


def test_truncated_poly_closed_form():
    # k[x]/x^2: HH_* char 0 proxy = [2,1,1,1,...]; char 2 = [2,2,2,...]
    A = truncated_polynomial(2)
    h = minimal_homology_dims(A, 6, primes=(32003, 2))
    assert h[32003] == [2, 1, 1, 1, 1, 1, 1]
    assert h[2] == [2, 2, 2, 2, 2, 2, 2]


def test_contracted_degree_matches_full_complex():
    """_contracted_degree(eng, cols[n], rks[n-1], n) == _contracted_complex(...)[n]."""
    import numpy as np
    from quiverlab.engine.resolutions_minimal import (
        minimal_resolution, _contracted_complex, _contracted_degree)
    A = truncated_polynomial(3)            # k[x]/x^3, a small exact oracle
    p = 32003
    rks, cols, eng, _ = minimal_resolution(A, 4, p)
    Dbar = _contracted_complex(A, rks, cols, eng, 4)
    for n in range(1, 5):
        got = _contracted_degree(eng, cols.get(n, []) or [], rks.get(n - 1, 0), n)
        assert np.array_equal(got % p, Dbar[n] % p), f"degree {n} mismatch"


def test_minimal_resolution_golden_unchanged():
    """Pin minimal_resolution output (refactor must be behavior-preserving)."""
    from quiverlab.engine.resolutions_minimal import minimal_resolution
    A = truncated_polynomial(3)
    rks, cols, eng, trunc = minimal_resolution(A, 5, 32003)
    assert rks[0] == 1 and rks[1] == 1 and rks[2] == 1          # k[x]/x^3 is periodic, r_n=1
    assert trunc is None


def test_advance_resolution_matches_loop():
    """Stepping _advance_resolution reproduces minimal_resolution's rks/cols."""
    from quiverlab.engine.resolutions_minimal import (
        minimal_resolution, _init_resolution, _advance_resolution)
    A = truncated_polynomial(3)
    p = 32003
    rks_ref, cols_ref, _, _ = minimal_resolution(A, 4, p)
    st = _init_resolution(A, p)
    for _ in range(1, 6):
        res = _advance_resolution(st, p, 20000, None)
        if res["status"] != "ok":
            break
    for n in (0, 1, 2, 3):
        assert st["rks"].get(n) == rks_ref.get(n)
