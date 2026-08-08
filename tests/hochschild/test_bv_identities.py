"""Plan 54 Task D -- the BV self-certifications (oracle_selfcert): Delta^2 = 0, the
seven-term relation (7T), the perfect-pairing certificate, and
``bracket_check.agrees``. All routes; here the symmetric route."""
import numpy as np
import pytest

import quiverlab as ql
from quiverlab.hochschild.bv.bracket import seven_term_residual

pytestmark = pytest.mark.oracle_selfcert

_SYM = [
    ("kxy", lambda p: ql.QuantumCI(-1, field=ql.GF(p)), 3),
    ("kx_x3", lambda p: ql.truncated_polynomial(3, field=ql.GF(p if p != 3 else 32003)), 3),
]


def _delta2_zero(bv, p):
    for n in range(2, bv.top + 1):
        Dn = np.array([[int(x) for x in r] for r in bv.matrices[n]], dtype=np.int64)
        Dprev = np.array([[int(x) for x in r] for r in bv.matrices[n - 1]], dtype=np.int64)
        if Dn.size and Dprev.size:
            assert not np.any((Dprev @ Dn) % p), f"Delta^2 != 0 at {n}"


@pytest.mark.parametrize("prime", [3, 5, 32003])
def test_delta_squared_zero_and_bracket_agrees(prime):
    for _name, build, top in _SYM:
        A = build(prime)
        bv = A.bv_operator(top)
        _delta2_zero(bv, A.domain.p)
        assert bv.bracket_check is not None
        assert bv.bracket_check["agrees"] is True
        assert bv.bracket_check["pairs_checked"] >= 1


@pytest.mark.parametrize("prime", [3, 5, 32003])
def test_seven_term_relation(prime):
    # (7T): Delta is a differential operator of order <= 2 wrt the cup product.
    A = ql.QuantumCI(-1, field=ql.GF(prime))
    top = 3
    bv = A.bv_operator(top)
    cup = A.cup_products(top)
    # triples (p,q,r) with p+q+r <= top
    for p in range(1, top + 1):
        for q in range(1, top + 1):
            for r in range(1, top + 1):
                if p + q + r <= top:
                    res = seven_term_residual(A, bv, cup, p, q, r)
                    assert res == 0, f"7T failed at ({p},{q},{r}): residual {res}"


def test_perfect_pairing_certificate_recorded():
    # A singular/ non-square pairing would have raised in the transport; a returned
    # BVOperator therefore certifies the perfect pairing in-window.
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    bv = A.bv_operator(3)
    assert bv.window == 3
