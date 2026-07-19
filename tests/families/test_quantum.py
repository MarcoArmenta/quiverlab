"""QuantumCI k<x,y>/(x^2, y^2, xy + q yx). Fixtures N6, N7."""
from fractions import Fraction

from quiverlab.families import ExteriorAlgebra, QuantumCI
from quiverlab.fields import CC, GF


def test_quantum_ci_i_dim4_hh0_2():
    A = QuantumCI(q="i", field=CC)                     # Fixture N7
    assert A.dim == 4
    assert A.hochschild_cohomology(0).dims == [2]      # Z = span{1, xy}


def test_quantum_ci_rational_q():
    assert QuantumCI(q=3, field=GF(32003)).dim == 4


def test_exterior_2_equals_quantum_ci_1():
    E2 = ExteriorAlgebra(2, field=CC)
    Q1 = QuantumCI(q=1, field=CC)                      # Fixture N6
    assert E2.dim == Q1.dim == 4
    assert E2.T == Q1.T and E2.unit == Q1.unit         # byte-identical structure constants


def test_quantum_ci_minus_one_is_commutative_dim4_hh0_4():
    A = QuantumCI(q=-1, field=CC)                       # k[x,y]/(x^2,y^2), commutative
    assert A.dim == 4
    assert A.hochschild_cohomology(0).dims == [4]


def test_quantum_ci_negative_fraction_dim4_hh0_2():
    # q = -3/2 exercises BOTH _q_token deviations: str-not-repr (a Fraction must
    # stringify to '3/2', not 'Fraction(3, 2)') AND the sign-fold (the leading '-'
    # folds into the relation's subtraction, 'x*y - 3/2*y*x'). Building dim 4 proves
    # the relation parsed. Center rule: q != -1 (generic) => HH0 = [2].
    A = QuantumCI(q=Fraction(-3, 2), field=CC)
    assert A.dim == 4
    assert A.hochschild_cohomology(0).dims == [2]


def test_quantum_ci_root_of_unity_e3_dim4_hh0_2():
    # q = E(3), a primitive cube root of unity, routed through the T4 coefficient
    # grammar as a bare non-rational token. Center rule: q != -1 (generic) => HH0 = [2].
    A = QuantumCI(q="E(3)", field=CC)
    assert A.dim == 4
    assert A.hochschild_cohomology(0).dims == [2]
