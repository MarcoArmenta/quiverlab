import pytest

from quiverlab import GF, Quiver
from quiverlab.modules.morphism import direct_sum, is_direct_summand

pytestmark = pytest.mark.oracle_selfcert


def _a2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(7))


def test_direct_sum_identities():
    A = _a2()
    P1, S2, S1 = A.projective(1), A.simple(2), A.simple(1)
    D, incls, projs = direct_sum(P1, S2, S1)
    assert D.dim == P1.dim + S2.dim + S1.dim
    for i, (inc, prj) in enumerate(zip(incls, projs)):
        assert inc.then(prj).matrix == inc.src.identity_hom().matrix
    # sum of the idempotents is the identity on D
    dom = D.domain
    total = [[dom.zero()] * D.dim for _ in range(D.dim)]
    for inc, prj in zip(incls, projs):
        e = prj.then(inc)          # D -> Mi -> D
        for r in range(D.dim):
            for c in range(D.dim):
                total[r][c] = dom.add(total[r][c], e.matrix[r][c])
    assert total == D.identity_hom().matrix


def test_is_direct_summand_true_and_false():
    A = _a2()
    P1, S1, S2 = A.projective(1), A.simple(1), A.simple(2)
    D, _, _ = direct_sum(P1, S1)
    assert is_direct_summand(S1, D) is True
    assert is_direct_summand(P1, D) is True
    # S2 = rad P1 is a SUBmodule of P1 but P1 is indecomposable:
    assert is_direct_summand(S2, P1) is False
