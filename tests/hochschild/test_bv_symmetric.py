"""Plan 54 Task B -- the symmetric-route BV transport (Tradler anchor, P52-free).

Delta via the symmetric trace-form pairing + the ORDINARY Connes B. Asserts the
BVOperator shape, the char-sensitive hh_dims, Delta^2 = 0, and (implicitly, inside
the transport) the perfect-pairing invertibility self-cert.
"""
import numpy as np
import pytest

import quiverlab as ql
from quiverlab.hochschild.bv.transport import bv_matrices_symmetric
from quiverlab.hochschild.products import BVOperator

pytestmark = pytest.mark.oracle_selfcert


def _delta_squared_zero(bv, p):
    for n in range(2, bv.top + 1):
        Dn = np.array([[int(x) for x in row] for row in bv.matrices[n]], dtype=np.int64)
        Dprev = np.array([[int(x) for x in row] for row in bv.matrices[n - 1]],
                         dtype=np.int64)
        if Dn.size == 0 or Dprev.size == 0:
            continue
        assert not np.any((Dprev @ Dn) % p), f"Delta^2 != 0 at degree {n}"


def test_kx_x3_char_not_dividing_N():
    # k[x]/(x^3) over GF(32003): char does not divide N=3 -> HH^0 = 3, HH^n = 2.
    A = ql.truncated_polynomial(3, field=ql.GF(32003))
    bv = bv_matrices_symmetric(A, 4)
    assert isinstance(bv, BVOperator)
    assert bv.hypothesis == "symmetric (Tradler AIF 2008)"
    assert bv.hh_dims == [3, 2, 2, 2, 2]
    assert bv.basis == "bar/GF(32003)"
    # Delta_n shape: dim HH^{n-1} x dim HH^n, rows = output degree.
    assert len(bv.matrices[1]) == 3 and len(bv.matrices[1][0]) == 2
    _delta_squared_zero(bv, 32003)


def test_kx_x3_char_dividing_N():
    # k[x]/(x^3) over GF(3): char | N=3 -> HH^0 = 3, HH^n = 3 (the char-sensitive
    # jump). Served on the SAME symmetric route.
    A = ql.truncated_polynomial(3, field=ql.GF(3))
    bv = bv_matrices_symmetric(A, 3)
    assert bv.hh_dims == [3, 3, 3, 3]
    assert bv.hypothesis == "symmetric (Tradler AIF 2008)"
    _delta_squared_zero(bv, 3)


def test_kxy_x2y2_tradler_dims():
    # k[x,y]/(x^2,y^2): symmetric, HH_• = [4,4,5,6] pin (existing).
    A = ql.QuantumCI(-1, field=ql.GF(5))
    bv = bv_matrices_symmetric(A, 3)
    assert bv.hh_dims == [4, 4, 5, 6]
    _delta_squared_zero(bv, 5)


def test_dual_numbers_degree0_is_2():
    # k[x]/(x^2): HH^0 = 2 (= dim A, commutative), HH^n = 1 for n >= 1 (Plan-35
    # correction). Delta_1 is 2x1.
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    bv = bv_matrices_symmetric(A, 4)
    assert bv.hh_dims == [2, 1, 1, 1, 1]
    assert len(bv.matrices[1]) == 2 and len(bv.matrices[1][0]) == 1
    _delta_squared_zero(bv, 7)


def test_blocks_shape():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    bv = bv_matrices_symmetric(A, 3)
    b = bv.blocks()
    assert b["kind"] == "bv_operator"
    assert b["hh_dims"] == [2, 1, 1, 1]
    assert set(b["matrices"]) == {"1", "2", "3"}
    assert b["nakayama"]["inner"] is True
    assert "references" in b and b["window"] == 3
