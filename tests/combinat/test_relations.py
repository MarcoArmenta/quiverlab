from fractions import Fraction as F

import pytest
from quiverlab.combinat import Quiver
from quiverlab.combinat.relations import parse_relation
from quiverlab.errors import RelationError


def _q():
    return Quiver(vertices=[1, 2, 3],
                  arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3), "x": (1, 1)})


def test_monomial():
    r = parse_relation("a*b", _q())
    assert r.is_monomial and r.terms == ((F(1), ("a", "b")),)
    assert (r.source, r.target) == (1, 3)


def test_power():
    r = parse_relation("x^3", _q())
    assert r.terms == ((F(1), ("x", "x", "x")),)


def test_linear_combination_parallel():
    r = parse_relation("a*b - 2*c", _q())
    assert not r.is_monomial
    assert r.terms == ((F(1), ("a", "b")), (F(-2), ("c",)))
    assert (r.source, r.target) == (1, 3)


def test_fraction_coefficient():
    r = parse_relation("1/2*a*b + c", _q())
    assert r.terms[0][0] == F(1, 2)


def test_not_composable():
    with pytest.raises(RelationError) as ei:
        parse_relation("b*a", _q())
    assert "target" in str(ei.value)


def test_not_parallel():
    with pytest.raises(RelationError) as ei:
        parse_relation("a - c", _q())
    assert "parallel" in str(ei.value)


def test_unknown_arrow():
    with pytest.raises(RelationError):
        parse_relation("a*z", _q())


def test_decimal_coefficient_fails():
    with pytest.raises(Exception):  # ExactnessError via parse_rational
        parse_relation("0.5*a*b", _q())
