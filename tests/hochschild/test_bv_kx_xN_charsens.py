"""Plan 54 Task F -- k[x]/(x^N) BV char-sensitivity (oracle_literature).

The HH dimension pattern of k[x]/(x^N) is char-sensitive (BIKLZ arXiv:2603.04834
Sec 3.2 e=1, and the classical truncated-polynomial computation):

    dim HH^0 = N;   for n >= 1,   dim HH^n = N - 1  if  char k does not divide N,
                                             = N      if  char k divides N.

k[x]/(x^N) is SYMMETRIC (commutative Frobenius => nu = id) in EVERY characteristic,
so the BV operator is served on the symmetric/Tradler route in BOTH regimes with
Delta^2 = 0. This battery pins the char-sensitive dims (both regimes) and the
square-zero self-cert."""
import numpy as np
import pytest

import quiverlab as ql

pytestmark = pytest.mark.oracle_literature


def _dd_zero(bv, p):
    for n in range(2, bv.top + 1):
        Dn = np.array([[int(x) for x in r] for r in bv.matrices[n]], dtype=np.int64)
        Dp = np.array([[int(x) for x in r] for r in bv.matrices[n - 1]], dtype=np.int64)
        if Dn.size and Dp.size:
            assert not np.any((Dp @ Dn) % p), f"Delta^2 != 0 at {n}"


@pytest.mark.parametrize("N", [3, 4, 5])
def test_char_not_dividing_N(N):
    # char = 32003 does not divide small N -> HH^0 = N, HH^{n>=1} = N - 1.
    # top=3 (top=4 for N=5 exceeds the bar (b,B) max_cells guard on the symmetric
    # Connes B; the char-sensitive pattern is already conclusive at top=3).
    A = ql.truncated_polynomial(N, field=ql.GF(32003))
    bv = A.bv_operator(3)
    assert bv.hypothesis == "symmetric (Tradler AIF 2008)"
    assert bv.hh_dims == [N] + [N - 1] * 3
    _dd_zero(bv, 32003)


@pytest.mark.parametrize("N", [2, 3, 5])
def test_char_dividing_N(N):
    # char = N (prime) divides N -> HH^0 = N, HH^{n>=1} = N (the char-sensitive jump).
    A = ql.truncated_polynomial(N, field=ql.GF(N))
    bv = A.bv_operator(3)
    assert bv.hh_dims == [N, N, N, N]
    _dd_zero(bv, N)
