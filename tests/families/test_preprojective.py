"""PreprojectiveAlgebra Pi(Q). Fixture N8."""
from quiverlab.families import PreprojectiveAlgebra
from quiverlab.fields import CC


def test_preprojective_A2_dim4():
    A = PreprojectiveAlgebra("A2", field=CC)
    assert A.dim == 4                                   # monomial mesh relations
    assert A.hochschild_cohomology(0).dims == [1]


def test_preprojective_A3_dim10():
    A = PreprojectiveAlgebra("A3", field=CC)            # tetrahedral number T_3
    assert A.dim == 10                                  # Groebner-certified oracle
