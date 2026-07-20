"""ExteriorAlgebra(n) = Lambda(k^n), dim 2^n."""
from quiverlab.families import ExteriorAlgebra
from quiverlab.fields import CC, GF


def test_dims_are_powers_of_two():
    assert ExteriorAlgebra(2, field=CC).dim == 4
    assert ExteriorAlgebra(3, field=CC).dim == 8


def test_exterior_2_hh0_char_dependence():
    assert ExteriorAlgebra(2, field=CC).hochschild_cohomology(0).dims == [2]
    assert ExteriorAlgebra(2, field=GF(2)).hochschild_cohomology(0).dims == [4]  # degenerates
