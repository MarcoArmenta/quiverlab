"""Plan 54 Task D -- GF(2) IS SERVED (oracle_selfcert).

Rationale (spec Major-3): the transpose arrangement and the BV sign eps are a
GLOBAL code-path property fixed ONCE by the odd-prime arbiter, not a per-instance
choice; over GF(2) the identical code path runs and is certified per instance by
Delta^2 = 0, the BV relation mod 2 (a genuine constraint -- the (-1) factors
collapse but Delta(aUb) - Delta(a)Ub - aUDelta(b) must still equal the independent
GF(2) bracket), and pairing invertibility. HONEST-SCOPE CAVEAT: char-2
certification is Delta^2=0 + BV-relation + pairing-invertibility, NOT
sign-determination (GF(2) cannot arbitrate the signs, which collapse)."""
import numpy as np
import pytest

import quiverlab as ql

pytestmark = pytest.mark.oracle_selfcert


def test_bv_over_gf2_served():
    # products_loop_gf2 = k[x]/(x^3) over GF(2) (symmetric, commutative Frobenius).
    A = ql.truncated_polynomial(3, field=ql.GF(2))
    bv = A.bv_operator(3)
    # Delta^2 = 0
    for n in range(2, bv.top + 1):
        Dn = np.array([[int(x) for x in r] for r in bv.matrices[n]], dtype=np.int64)
        Dprev = np.array([[int(x) for x in r] for r in bv.matrices[n - 1]], dtype=np.int64)
        if Dn.size and Dprev.size:
            assert not np.any((Dprev @ Dn) % 2)
    # BV relation mod 2 vs the independent GF(2) bracket (pairing invertibility is
    # asserted inside the transport; a returned BVOperator certifies it).
    assert bv.bracket_check is not None
    assert bv.bracket_check["agrees"] is True
    assert bv.basis == "bar/GF(2)"
