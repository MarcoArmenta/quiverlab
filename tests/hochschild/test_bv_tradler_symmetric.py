"""Plan 54 Task F -- the Tradler symmetric anchor (oracle_literature).

k[x,y]/(x^2,y^2) is SYMMETRIC (nu = id), so its BV structure is Tradler's
(Ann. Inst. Fourier 58 (2008) 2351-2379): Delta is the dual of the ORDINARY
Connes B under the symmetric trace form. The HH_• = [4,4,5,6] dimension pin is the
existing Plan-33 truncated-tensor literature value; here we anchor that the BV
operator on those groups is nontrivial (Delta != 0) and that the derived bracket
equals the independently computed Gerstenhaber bracket (the sign convention is
arbitrated, not assumed)."""
import numpy as np
import pytest

import quiverlab as ql

pytestmark = pytest.mark.oracle_literature


def test_kxy_x2y2_bv_structure():
    A = ql.QuantumCI(-1, field=ql.GF(5))          # k[x,y]/(x^2,y^2)
    top = 3
    bv = A.bv_operator(top)
    assert bv.hypothesis == "symmetric (Tradler AIF 2008)"
    assert bv.hh_dims == [4, 4, 5, 6]             # HH_• literature pin
    # Delta is nontrivial (Tradler BV is not the zero operator here)
    nonzero = any(int(x) for n in bv.matrices for row in bv.matrices[n] for x in row)
    assert nonzero, "Tradler Delta should be nonzero on k[x,y]/(x^2,y^2)"
    # derived bracket == independent (the arbiter certifies the whole structure)
    assert bv.bracket_check["agrees"] is True
    p = 5
    for n in range(2, top + 1):
        Dn = np.array([[int(x) for x in r] for r in bv.matrices[n]], dtype=np.int64)
        Dp = np.array([[int(x) for x in r] for r in bv.matrices[n - 1]], dtype=np.int64)
        if Dn.size and Dp.size:
            assert not np.any((Dp @ Dn) % p)
