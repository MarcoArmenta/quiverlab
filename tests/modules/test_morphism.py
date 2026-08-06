"""First-class module homomorphisms (Plan 37 / C1). Self-certifying:
the constructor validates the intertwining relations exactly."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.modules.morphism import ModuleHom, hom_basis, zero_hom

pytestmark = pytest.mark.oracle_selfcert


def _a2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(7))


def test_hom_basis_wraps_hom_space_with_validation():
    A = _a2()
    P1, S2 = A.projective(1), A.simple(2)
    basis = hom_basis(P1, S2)
    assert len(basis) == A.hom(P1, S2)          # dims agree with the old API
    for f in basis:
        assert f.src is P1 and f.tgt is S2


def test_constructor_rejects_non_module_map():
    A = _a2()
    P1, S1 = A.projective(1), A.simple(1)       # dims 2 and 1
    # [[1, 1]]: P1 -> S1 hitting the vertex-2 coordinate -- breaks the e_1
    # idempotent intertwining relation, so the constructor must refuse.
    with pytest.raises(QuiverlabError, match="module map"):
        ModuleHom(P1, S1, [[1, 1]], check=True)


def test_then_composes_left_to_right():
    A = _a2()
    P1 = A.projective(1)
    top = P1.top()                               # = S1
    # P1 ->> top(P1): from the projective-cover direction we at least have
    # SOME epi in the hom basis; compose with an endo of top.
    fs = [f for f in hom_basis(P1, top) if f.is_epi()]
    assert fs, "an epi P1 ->> top(P1) must exist"
    f = fs[0]
    idt = top.identity_hom()
    g = f.then(idt)
    assert g.matrix == f.matrix and g.src is P1 and g.tgt is top


def test_then_refuses_mismatched_middle():
    A = _a2()
    P1, S2 = A.projective(1), A.simple(2)
    f = zero_hom(P1, S2)
    with pytest.raises(QuiverlabError, match="compose"):
        f.then(zero_hom(P1, S2))                 # tgt S2 != src P1


def test_mono_epi_iso_flags():
    A = _a2()
    S2 = A.simple(2)
    idm = S2.identity_hom()
    assert idm.is_mono() and idm.is_epi() and idm.is_iso()
    z = zero_hom(S2, S2)
    assert z.is_zero() and not z.is_mono()


def test_cross_algebra_refused():
    A, B = _a2(), _a2()
    with pytest.raises(QuiverlabError):
        hom_basis(A.simple(1), B.simple(1))
