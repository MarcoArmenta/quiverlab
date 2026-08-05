"""End(M) acts on Ext^1(M, N) by precomposition of the lifted chain map.
Self-cert: id acts as identity; the action is a representation of End(M)
(phi.psi acts as the product of the action matrices, left-to-right)."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.ar import (_ext1_action_matrix, _ext1_data,
                                  end_action_on_ext1)

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return linear_path_algebra(3, field=QQ)


def _action_of(phi, M, N, dom):
    """The same internal action matrix end_action_on_ext1 uses, for one endo phi."""
    cocycle_mats, coboundary_mats, terms, dmats = _ext1_data(M.algebra, M, N)
    return _ext1_action_matrix(phi, cocycle_mats, coboundary_mats, (terms, dmats), dom)


def test_identity_acts_as_identity():
    A = _a3()
    M, N = A.simple(2), A.simple(2).tau()          # Ext^1(M, tau M) is the AR case
    basis, action = end_action_on_ext1(M, N)
    e = len(action[0]) if action else 0
    idx = next(i for i, f in enumerate(basis) if f.is_iso() and f.matrix
               == M.identity_hom().matrix)
    assert action[idx] == lm.identity(e, A.domain)


def test_action_is_a_representation():
    A = _a3()
    M, N = A.simple(2), A.simple(2).tau()
    basis, action = end_action_on_ext1(M, N)
    dom = A.domain
    # pick two basis endos; their composite (left-to-right, f.then(g)) must act as
    # the product of action matrices (in the matching order).
    for i in range(len(basis)):
        for j in range(len(basis)):
            comp = basis[i].then(basis[j])          # M -> M
            direct = _action_of(comp, M, N, dom)    # recompute via ar
            combined = lm.matmul(action[i], action[j], dom)  # precompose order
            assert direct == combined
