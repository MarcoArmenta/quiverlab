"""SES objects: exactness certified at construction; split test; pushout/
pullback squares certified by their universal-square identities."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules.morphism import ModuleHom, hom_basis, zero_hom
from quiverlab.modules.ses import ShortExactSequence, pullback, pushout

pytestmark = pytest.mark.oracle_selfcert


def _a2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(7))


def _rad_ses(A):
    """0 -> rad P1 -> P1 -> S1 -> 0 for kA2 (rad P1 = S2, dim 1)."""
    P1, S1 = A.projective(1), A.simple(1)
    R = P1.radical()
    iota = next(f for f in hom_basis(R, P1) if f.is_mono())
    epi = next(f for f in hom_basis(P1, S1) if f.is_epi())
    return ShortExactSequence(iota, epi)


def test_ses_certifies_exactness():
    ses = _rad_ses(_a2())
    assert ses.M.dim == ses.L.dim + ses.N.dim


def test_ses_rejects_non_exact():
    A = _a2()
    S1, S2 = A.simple(1), A.simple(2)
    with pytest.raises(QuiverlabError, match="exact"):
        ShortExactSequence(zero_hom(S1, S2), zero_hom(S2, S1))


def test_split_iff_ext_vanishes():
    A = _a2()
    ses = _rad_ses(A)
    # Ext^1(S1, S2) = k for kA2 (the arrow) => the rad sequence is NOT split
    assert A.ext(ses.N, ses.L, 1) == 1
    assert ses.is_split() is False


def test_direct_sum_ses_splits():
    A = _a2()
    S1, S2 = A.simple(1), A.simple(2)
    from quiverlab.modules.morphism import direct_sum          # Task 4 name
    D, (i1, i2), (p1, p2) = direct_sum(S2, S1)
    ses = ShortExactSequence(i1, p2)
    assert ses.is_split() is True


def test_pushout_square_commutes():
    A = _a2()
    P1 = A.projective(1)
    R = P1.radical()
    f = next(h for h in hom_basis(R, P1) if h.is_mono())
    g = R.identity_hom()
    P, inB, inC = pushout(f, g)
    assert f.then(inB).matrix == g.then(inC).matrix
    assert P.dim == P1.dim + R.dim - R.dim                  # pushout along id: P ~ P1


def test_pullback_square_commutes():
    A = _a2()
    P1, S1 = A.projective(1), A.simple(1)
    f = next(h for h in hom_basis(P1, S1) if h.is_epi())
    g = S1.identity_hom()
    P, prB, prC = pullback(f, g)
    assert prB.then(f).matrix == prC.then(g).matrix
