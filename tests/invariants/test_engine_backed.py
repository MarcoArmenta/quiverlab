import pytest
from quiverlab import CC, GF, truncated_polynomial
from quiverlab.errors import FieldError


def test_cyclic_homology_gfp():
    A = truncated_polynomial(2, field=GF(3))
    t = A.cyclic_homology(4)
    assert t.kind == "HC_"
    assert len(t.dims) == 5 and all(isinstance(d, int) for d in t.dims)


def test_cyclic_homology_cc_loud():
    with pytest.raises(FieldError):
        truncated_polynomial(2, field=CC).cyclic_homology(4)


def test_symmetric_dual_numbers_gfp():
    A = truncated_polynomial(2, field=GF(5))
    assert A.is_frobenius() is True
    assert A.is_symmetric() is True


def test_frobenius_cc_loud():
    with pytest.raises(FieldError):
        truncated_polynomial(2, field=CC).is_frobenius()
