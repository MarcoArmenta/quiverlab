"""The injective coresolution now RETAINS its differentials (the transposes of
the minimal projective resolution of DM over A^op, in the dual bases) so
reports can show them. Certify the shipped matrices genuinely form the exact
coresolution 0 -> M -> E^0 -> E^1 -> ...: iota is injective, consecutive maps
compose to zero, and the rank bookkeeping is exact at every term."""
import pytest

import quiverlab as ql
from quiverlab.modules import linalg_mod as lm

pytestmark = pytest.mark.oracle_selfcert


def _rank(D, dom):
    return len(lm.column_space_pivots(D, dom)) if (D and D[0]) else 0


@pytest.fixture(scope="module", params=["nakayama", "hereditary"])
def resolved(request):
    if request.param == "nakayama":
        # self-injective: the coresolution never terminates in the window
        A = ql.NakayamaAlgebra(n=4, l=3, cyclic=True, field=ql.GF(7))
    else:
        # hereditary kA3: everything resolves in one step
        A = ql.PathAlgebra("A3", field=ql.GF(5))
    M = A.projective(1).radical()       # interior uniserial: not simple/proj/inj
    assert M.dim > 0
    return A, M, M.injective_resolution(4)


def test_iota_is_injective_with_the_right_shape(resolved):
    A, M, res = resolved
    D0 = res.differential(0)            # iota: M -> E^0
    assert len(D0) == res.terms[0].dim and len(D0[0]) == M.dim
    assert _rank(D0, A.domain) == M.dim


def test_consecutive_differentials_compose_to_zero(resolved):
    A, M, res = resolved
    dom = A.domain
    for n in range(len(res.dmats) - 1):
        Dn, Dnext = res.differential(n), res.differential(n + 1)
        if not (Dn and Dn[0] and Dnext and Dnext[0]):
            continue
        prod = lm.matmul(Dnext, Dn, dom)
        z = dom.zero()
        assert all(x == z for row in prod for x in row), \
            f"d^{n + 1} o d^{n} != 0"


def test_rank_bookkeeping_is_exact_at_every_term(resolved):
    # exactness at E^n: im(d^n) = ker(d^{n+1}), i.e. rank d^n + rank d^{n+1}
    # = dim E^n (with d^0 = iota).
    A, M, res = resolved
    dom = A.domain
    for n in range(len(res.dmats) - 1):
        En = res.terms[n]
        if En is None or En.dim == 0:
            break
        r_in = _rank(res.differential(n), dom)
        r_out = _rank(res.differential(n + 1), dom)
        assert r_in + r_out == En.dim, f"not exact at E^{n}"
