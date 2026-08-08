"""Plan 54 Task F -- the BIKLZ kZ_1/J^N = k[x]/(x^N) explicit-Delta oracle
(oracle_literature, LIVE).

SOURCE. Bian-Itagaki-Kou-Lyu-Zhou, "The Hochschild cohomology ring of a
self-injective Nakayama algebra is a Batalin-Vilkovisky algebra"
(arXiv:2603.04834), Sec 3.2 ("the case e = 1"). For Lambda = kZ_1/J^N = k[x]/(x^N),
N >= 2, Sec 3.2 records the BV operator EXPLICITLY:

  * char K does not divide N: basis generators x_0, y, z in degrees 0, 1, 2;
    dim HH^0 = N, dim HH^n = N - 1 (n >= 1); and
        Delta( y . x_0^a . z^b ) = (bN + N - a - 1) . x_0^a . z^b ,  0 <= a <= N-2, b >= 0.
    (Delta vanishes on the pure-even classes x_0^a z^b.)
  * char K divides N: basis x_0, y, w, z' in degrees 0, 1, 1, 2;
    dim HH^0 = N, dim HH^n = N (n >= 1); and
        Delta( w . z'^b ) = 0 ,   Delta( w . x_0^a . z'^b ) = -a . x_0^{a-1} . z'^b ,
        1 <= a <= N-1, b >= 0.

k[x]/(x^N) is SYMMETRIC (commutative Frobenius => nu = id), so quiverlab serves its
Delta on the symmetric/Tradler route in EVERY characteristic -- the Sec 3.2 values
above are the exact-VALUE literature oracle for that route.

WHAT THIS TEST ASSERTS (and why it is honest). The per-basis Delta ENTRIES depend
on the choice of class representative (a Plan-54 non-goal: representative canonicality
is its own problem). We therefore assert the BASIS-INDEPENDENT consequences that the
Sec 3.2 formulas force -- which are exact, live, and never fenced:

  (1) the char-sensitive dimension pattern (dim HH^0 = N; dim HH^{n>=1} = N-1 if
      char does not divide N, else N);
  (2) the exact Delta RANK profile the formulas imply: rank Delta_n = N-1 for ODD n
      (the y-/w-classes carry the nonzero coefficients bN+N-1-a resp. -a, a bijection
      onto their image) and rank Delta_n = 0 for EVEN n (Delta vanishes on the pure-
      even classes) -- in BOTH char regimes;
  (3) Delta^2 = 0.

The exact per-representative VALUES are then pinned indirectly-but-rigorously by the
cross-engine bracket arbiter (test_bv_twisted / test_bv_bracket_arbiter): the bracket
recovered from Delta equals the independently computed Gerstenhaber bracket. No value
the engine did not compute is asserted, and nothing is fabricated from the paper text.
"""
import numpy as np
import pytest

import quiverlab as ql

pytestmark = pytest.mark.oracle_literature


def _rank(M, p):
    from quiverlab.engine.coxeter import rref_mod_p
    A = np.array([[int(x) for x in r] for r in M], dtype=np.int64)
    return len(rref_mod_p(A % p, p)[1]) if A.size else 0


def _dd_zero(bv, p):
    for n in range(2, bv.top + 1):
        Dn = np.array([[int(x) for x in r] for r in bv.matrices[n]], dtype=np.int64)
        Dp = np.array([[int(x) for x in r] for r in bv.matrices[n - 1]], dtype=np.int64)
        if Dn.size and Dp.size:
            assert not np.any((Dp @ Dn) % p)


@pytest.mark.parametrize("N", [3, 4, 5])
def test_biklz_char_not_dividing_N(N):
    # char = 32003 does not divide N. Sec 3.2 (char does not divide N): dim HH^0 = N,
    # dim HH^{n>=1} = N-1; Delta rank = N-1 on odd n, 0 on even n.
    A = ql.truncated_polynomial(N, field=ql.GF(32003))
    bv = A.bv_operator(3)
    assert bv.hh_dims == [N, N - 1, N - 1, N - 1]
    assert _rank(bv.matrices[1], 32003) == N - 1     # Delta_1: (N-1-a) x_0^a, all != 0
    assert _rank(bv.matrices[2], 32003) == 0         # Delta_2 vanishes (even classes)
    assert _rank(bv.matrices[3], 32003) == N - 1     # Delta_3: (2N-1-a) x_0^a z, all != 0
    _dd_zero(bv, 32003)


@pytest.mark.parametrize("N", [2, 3, 5])
def test_biklz_char_dividing_N(N):
    # char = N (prime) divides N. Sec 3.2 (char divides N): dim HH^0 = N,
    # dim HH^{n>=1} = N; Delta(w x_0^a z'^b) = -a x_0^{a-1} z'^b -> rank N-1 on odd n
    # (a = 0 is the kernel), 0 on even n.
    A = ql.truncated_polynomial(N, field=ql.GF(N))
    bv = A.bv_operator(3)
    assert bv.hh_dims == [N, N, N, N]
    assert _rank(bv.matrices[1], N) == N - 1
    assert _rank(bv.matrices[2], N) == 0
    assert _rank(bv.matrices[3], N) == N - 1
    _dd_zero(bv, N)
