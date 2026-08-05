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
