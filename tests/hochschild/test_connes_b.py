"""Induced Connes B on HH (Plan 35): B^2 = 0 at the induced level, rank
consistency with the (b,B) cyclic dims, and GF(p)/generic agreement."""
import pytest

import quiverlab as ql
from quiverlab.fields import QQ

pytestmark = [pytest.mark.oracle_selfcert]


@pytest.fixture(scope="module", params=[ql.GF(7), QQ], ids=["GF7", "QQ"])
def A(request):
    return ql.truncated_polynomial(2, field=request.param)


def test_induced_B_squares_to_zero(A):
    from quiverlab.hochschild.products import connes_b_tables
    cb = connes_b_tables(A, 3, max_cells=4_000_000)
    # composite HH_n -> HH_{n+1} -> HH_{n+2} must vanish entry-wise
    dom = A.domain
    for n in range(2):
        M1, M2 = cb.matrices[n], cb.matrices[n + 1]
        rows, mid, cols = len(M2), len(M1), (len(M1[0]) if M1 else 0)
        for r in range(rows):
            for c in range(cols):
                acc = dom.zero()
                for k in range(mid):
                    acc = dom.add(acc, dom.mul(dom.coerce(M2[r][k]),
                                               dom.coerce(M1[k][c])))
                assert dom.is_zero(acc), f"(B∘B)[{r}][{c}] != 0 at n={n}"


def test_ranks_and_dims_recorded(A):
    from quiverlab.hochschild.products import connes_b_tables
    cb = connes_b_tables(A, 2, max_cells=4_000_000)
    hh = list(A.hochschild_homology(2, verbose=False).dims)
    assert cb.hh_dims == hh
    assert set(cb.matrices) == {0, 1} and set(cb.ranks) == {0, 1}
    for n in (0, 1):
        assert len(cb.matrices[n]) == hh[n + 1]
        assert 0 <= cb.ranks[n] <= min(hh[n], hh[n + 1])


def test_gfp_and_generic_ranks_agree():
    from quiverlab.hochschild.products import connes_b_tables
    A7 = ql.truncated_polynomial(3, field=ql.GF(32003))
    Aq = ql.truncated_polynomial(3, field=QQ)
    r7 = connes_b_tables(A7, 3, max_cells=4_000_000).ranks
    rq = connes_b_tables(Aq, 3, max_cells=4_000_000).ranks
    assert r7 == rq          # char-0-shaped p: ranks agree (32003 is the big prime)
