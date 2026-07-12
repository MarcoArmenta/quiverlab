from fractions import Fraction

import pytest
from quiverlab.errors import ExactnessError
from quiverlab.fields import QQ


def test_coerce_and_arithmetic():
    a = QQ.coerce("1/3")
    b = QQ.coerce(2)
    assert QQ.eq(QQ.add(a, b), Fraction(7, 3))
    assert QQ.eq(QQ.mul(a, QQ.inv(a)), QQ.one())
    assert QQ.is_zero(QQ.sub(b, QQ.coerce(Fraction(2))))
    assert QQ.characteristic == 0
    assert QQ.to_str(a) == "1/3"


def test_floats_fail_loudly():
    with pytest.raises(ExactnessError):
        QQ.coerce(0.5)
    with pytest.raises(ExactnessError):
        QQ.coerce("0.5")


def test_scientific_notation_fails_loudly():
    with pytest.raises(ExactnessError):
        QQ.coerce("15e-1")
    with pytest.raises(ExactnessError):
        QQ.coerce("1E-3")


def test_bool_rejected():
    with pytest.raises(ExactnessError):
        QQ.coerce(True)
