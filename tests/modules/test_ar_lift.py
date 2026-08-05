"""General chain-map lift of a module endomorphism over its minimal resolution.
Self-certifying: every square d_n . phi_n == phi_{n-1} . d_n is asserted exactly,
and identity lifts to identities."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.ar import lift_endomorphism_along_resolution
from quiverlab.modules.morphism import hom_basis
from quiverlab.modules.resolution import minimal_resolution

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(32003))


def _squares_commute(phis, dmats, dom):
    # eps . phi_0 == phi . eps is the caller's cover condition; here check the
    # interior chain-map squares d_n . phi_n == phi_{n-1} . d_n.
    for n in range(1, len(phis)):
        dn = dmats[n]
        lhs = lm.matmul(dn, phis[n], dom)
        rhs = lm.matmul(phis[n - 1], dn, dom)
        if lhs != rhs:
            return False
    return True


def test_identity_lifts_to_identities():
    A = _a3()
    for v in (1, 2, 3):
        M = A.projective(v) if v == 1 else A.simple(v)
        idM = M.identity_hom()
        phis = lift_endomorphism_along_resolution(idM, degrees=3)
        for n, ph in enumerate(phis):
            assert ph == lm.identity(len(ph), M.domain)


def test_endomorphism_squares_commute():
    A = _a3()
    S2 = A.simple(2)
    terms, dmats = minimal_resolution(S2, 4)
    for phi in hom_basis(S2, S2):                 # every endo of S2
        phis = lift_endomorphism_along_resolution(phi, (terms, dmats), degrees=3)
        # cover condition eps . phi_0 == phi . eps
        eps = dmats[0]
        assert lm.matmul(eps, phis[0], A.domain) == lm.matmul(phi.matrix, eps,
                                                              A.domain)
        assert _squares_commute(phis, dmats, A.domain)


def test_lift_is_deterministic():
    A = _a3()
    P1 = A.projective(1)
    phi = hom_basis(P1, P1)[0]
    a = lift_endomorphism_along_resolution(phi, degrees=3)
    b = lift_endomorphism_along_resolution(phi, degrees=3)
    assert a == b                                 # reduce_mod_nullspace => byte-stable
