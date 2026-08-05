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
    A = truncated_polynomial(3, field=GF(2))
    from quiverlab.errors import QuiverlabError
    M = A.simple(1).syzygy()                            # dim 2 over GF(2): char <= dim
    with pytest.raises(QuiverlabError):
        igusa_todorov_phi(M)
