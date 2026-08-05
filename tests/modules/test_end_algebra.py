"""End(M) as a structure-constant Algebra. Self-certified by
from_structure_constants(check=True) (associativity + unit); the
regular-module oracle pins the composition-order convention."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.modules.endomorphism import end_algebra

pytestmark = pytest.mark.oracle_selfcert


def _a2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(7))


def test_end_of_simple_is_the_field():
    A = _a2()
    E = end_algebra(A.simple(1))
    assert E.dim == 1


def test_end_dim_matches_hom_dim():
    A = _a2()
    P1 = A.projective(1)
    E = end_algebra(P1)
    assert E.dim == A.hom(P1, P1)
    E._validate()                      # associativity + unit, loud on failure


def test_end_of_regular_module_has_algebra_dimension():
    # End_A(A_A) ~ A as a k-algebra (dim check + Loewy length agree)
    A = _a2()
    from quiverlab.modules.morphism import direct_sum
    regular, _, _ = direct_sum(A.projective(1), A.projective(2))
    E = end_algebra(regular)
    assert E.dim == 3                  # dim kA2 = 3
    assert E.loewy_length() == A.loewy_length() == 2


def test_end_of_regular_is_A_not_Aop():
    # Devil's-advocate arbiter (2026-08-05): kA2 is self-opposite, so the
    # original regular-module oracle could not see the product order. On
    # the source quiver 1->2, 1->3 (A != A^op) the corner dims of
    # End((+)P_v), computed THROUGH the structure constants, must equal
    # the Cartan matrix (composition order => End(A_A) ~ A); the flipped
    # order gives the transpose.
    from quiverlab import Quiver, GF
    from quiverlab.modules.endomorphism import regular_corner_dims

    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (1, 3)}).algebra(field=GF(7))
    C = A.cartan_matrix()
    assert C != [list(r) for r in zip(*C)]        # genuinely non-symmetric
    assert regular_corner_dims(A) == C


def test_end_algebra_small_char_degrades_loudly():
    # Devil's-advocate fix (2026-08-05): char <= dim M means the trace-form
    # radical is uncertified -- end_algebra returns a label-less algebra and
    # loewy_length must REFUSE loudly, never TypeError.
    import pytest as _pytest
    from quiverlab import Quiver, GF
    from quiverlab.errors import QuiverlabError
    from quiverlab.modules.morphism import direct_sum

    A = Quiver([1, 2], {"a": (1, 2)}).algebra(field=GF(2))
    D, _, _ = direct_sum(A.projective(1), A.projective(2))    # dim 3, char 2
    E = D.end_algebra()
    assert E.basis_labels is None
    with _pytest.raises(QuiverlabError, match="label"):
        E.loewy_length()
