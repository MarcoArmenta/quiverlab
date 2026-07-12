import pytest
from quiverlab.errors import ExactnessError, FieldError
from quiverlab.fields.primefield import PrimeField


def test_arithmetic_mod_5():
    F5 = PrimeField(5)
    assert F5.characteristic == 5
    assert F5.add(3, 4) == 2
    assert F5.mul(3, 4) == 2
    assert F5.inv(3) == 2
    assert F5.is_zero(F5.sub(2, 2))
    assert F5.coerce(-1) == 4
    assert F5.coerce("1/3") == F5.mul(1, F5.inv(3))


def test_denominator_divisible_by_p():
    F5 = PrimeField(5)
    with pytest.raises(FieldError):
        F5.coerce("1/5")


def test_nonprime_rejected():
    with pytest.raises(FieldError):
        PrimeField(6)


def test_float_rejected():
    with pytest.raises(ExactnessError):
        PrimeField(5).coerce(0.2)
