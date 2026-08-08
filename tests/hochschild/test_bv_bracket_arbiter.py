"""Plan 54 Task D -- the bracket ARBITER (oracle_crossengine): the bracket derived
from Delta via the BV relation equals the INDEPENDENT Gerstenhaber bracket
(quiverlab computes the latter with NO Frobenius input) table-for-table, in-window,
over ODD primes (3, 5, 32003) -- where the (-1) factors are visible and pin the
sign/twist convention."""
import pytest

import quiverlab as ql

pytestmark = pytest.mark.oracle_crossengine

_ALGS = [
    ("kxy", lambda p: ql.QuantumCI(-1, field=ql.GF(p))),
    ("kx_x4", lambda p: ql.truncated_polynomial(4, field=ql.GF(p))),
]


@pytest.mark.parametrize("prime", [3, 5, 32003])
@pytest.mark.parametrize("name,build", _ALGS)
def test_derived_bracket_equals_independent(name, build, prime):
    A = build(prime)
    top = 3
    bv = A.bv_operator(top)
    independent = A.gerstenhaber_brackets(top)
    derived = bv.derived_bracket
    assert derived is not None
    assert bv.bracket_check["agrees"] is True
    p = A.domain.p
    checked = 0
    for key, dt in derived.tables.items():
        it = independent.tables[key]
        assert dt.dims == it.dims, key
        dl, dr, dout = dt.dims
        for k in range(dout):
            for i in range(dl):
                for j in range(dr):
                    assert int(dt.constants[k][i][j]) % p == \
                        int(it.constants[k][i][j]) % p, (key, k, i, j)
        checked += 1
    assert checked >= 1
