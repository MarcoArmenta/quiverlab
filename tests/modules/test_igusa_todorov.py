"""Igusa-Todorov phi/psi. Literature closed forms (Barrios-Mata 2019 for
truncated path algebras) + the phi=psi=pd identity for finite proj. dimension
(Igusa-Todorov 2005) + additivity. Over QQ/GF(32003) -- the decompose char
caveat forbids char <= dim M."""
import pytest

from quiverlab import GF, TruncatedPathAlgebra, linear_path_algebra, \
    truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.homdims import igusa_todorov_phi, igusa_todorov_psi

pytestmark = pytest.mark.oracle_literature


def test_phi_equals_pd_when_finite():
    # kA3/rad^2 (a*b=0): Nakayama, gl.dim 2. pd S1 = 2, pd S2 = 1, pd S3 = 0.
    A = TruncatedPathAlgebra("A3", 2, field=QQ)
    for v, expect in ((1, 2), (2, 1), (3, 0)):
        S = A.simple(v)
        pd = S.projective_resolution(8).pd()
        assert pd == expect
        assert igusa_todorov_phi(S) == pd              # phi = pd for finite pd
        assert igusa_todorov_psi(S) == pd              # psi = pd for finite pd


def test_self_injective_truncated_all_zero():
    # k[x]/(x^a) is self-injective: every module has phi = psi = 0
    # (Omega permutes {M_1..M_{a-1}}, so ranks are stable at n=0; no finite-pd
    # non-projective summand, so fpd = 0). Barrios-Mata 2019.
    A = truncated_polynomial(4, field=QQ)
    for i in (1, 2, 3):
        M = A.simple(1)
        for _ in range(i - 1):
            M = M.syzygy()
        assert igusa_todorov_phi(M) == 0
        assert igusa_todorov_psi(M) == 0


def test_psi_geq_phi_and_projective_additivity():
    A = linear_path_algebra(3, field=QQ)               # hereditary kA3
    from quiverlab.modules.morphism import direct_sum
    S1, P3 = A.simple(1), A.projective(3)
    assert igusa_todorov_psi(S1) >= igusa_todorov_phi(S1)
    D, _, _ = direct_sum(S1, P3)                        # add a projective summand
    assert igusa_todorov_phi(D) == igusa_todorov_phi(S1)   # phi(M (+) P) = phi(M)


@pytest.mark.oracle_selfcert
def test_char_caveat_refuses_loudly():
    # over GF(2) with dim M >= 2 the decompose trace form is unreliable ->
    # the K0 bookkeeping must inherit the loud refusal, never a silent phi.
    # NOTE (devil's-advocate round 2026-08-05): k[x]/(x^3) no longer works
    # here -- it is SELF-INJECTIVE, so the theorem guard now returns phi = 0
    # before decompose is ever consulted (a correct value, not a refusal).
    # The refusal contract is exercised on a NON-self-injective algebra.
    from quiverlab import Quiver
    from quiverlab.errors import QuiverlabError
    # local radical-square-zero on two loops: NOT self-injective (socle dim 2),
    # P_1 has dim 3 with End(P_1) = A (dim 3 > 1), so decompose must consult
    # the trace form -- unreliable over GF(2) (char 2 <= 3) -> loud refusal.
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "x*y", "y*x", "y*y"], field=GF(2))
    M = A.projective(1)
    with pytest.raises(QuiverlabError):
        igusa_todorov_phi(M)


def test_phi_psi_zero_on_complexity_two_selfinjective():
    # Devil's-advocate oracle (2026-08-05): the exterior algebra Lambda(k^2)
    # is self-injective of complexity 2 -- its syzygies GROW forever
    # (dims 3, 5, 7, ...), so without the theorem guard the rank registry
    # never stabilizes. phi = psi = 0 by the stable-equivalence theorem,
    # and the call must be instant.
    from quiverlab import GF, Quiver
    from quiverlab.modules.homdims import igusa_todorov_phi, igusa_todorov_psi

    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "x*y+y*x"], field=GF(32003))
    S = A.simple(1)
    assert igusa_todorov_phi(S) == 0
    assert igusa_todorov_psi(S) == 0
    M = S.tau()   # another non-projective over a self-injective algebra
    if M.dim:
        assert igusa_todorov_phi(M) == 0
