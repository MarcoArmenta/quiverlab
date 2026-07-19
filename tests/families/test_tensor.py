"""TensorProduct(A, B) = A (x) B, Kunneth. Fixture N9."""
from quiverlab.combinat import Quiver
from quiverlab.families import TensorProduct
from quiverlab.families.basic import truncated_polynomial
from quiverlab.fields import CC


def _kA2(field=CC):
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=field)


def test_kA2_tensor_kA2_is_commutative_square_dim9():
    A = TensorProduct(_kA2(), _kA2())
    assert A.dim == 9
    assert A.hochschild_cohomology(2).dims == [1, 0, 0]     # == Plan-03 square == diamond


def test_kunneth_hh0_multiplicative():
    A = TensorProduct(truncated_polynomial(2, field=CC), truncated_polynomial(2, field=CC))
    assert A.dim == 4                                       # k[x,y]/(x^2,y^2)
    assert A.hochschild_cohomology(0).dims == [4]          # HH^0 = 2 * 2


def test_mismatched_fields_are_loud():
    import pytest
    from quiverlab.fields import GF
    from quiverlab.errors import FieldError
    with pytest.raises(FieldError):
        TensorProduct(_kA2(CC), _kA2(GF(5)))
