"""Sanctioned grammar extension: exact non-rational coefficients (spec §3.3;
foreseen in the Plan-03 boundary note). Rational relations are unchanged."""
from fractions import Fraction

import pytest
import sympy

from quiverlab.combinat import Quiver
from quiverlab.combinat.relations import parse_relation
from quiverlab.errors import ExactnessError, RelationError
from quiverlab.fields import CC


def _two_loops():
    return Quiver([1], {"x": (1, 1), "y": (1, 1)})


def test_rational_relations_are_byte_compatible():
    Q = _two_loops()
    r = parse_relation("x*y - 2*y*x", Q)
    assert r.terms == ((Fraction(1), ("x", "y")), (Fraction(-2), ("y", "x")))
    for c, _w in r.terms:
        assert isinstance(c, Fraction)                 # unchanged type


def test_imaginary_unit_coefficient_parses_exactly():
    Q = _two_loops()
    r = parse_relation("x*y + i*y*x", Q)               # QuantumCI(q="i")
    coeffs = {w: c for c, w in r.terms}
    assert coeffs[("x", "y")] == 1
    assert coeffs[("y", "x")] == sympy.I               # exact i, field-agnostic
    A = Q.algebra(relations=["x^2", "y^2", "x*y + i*y*x"], field=CC)   # general route
    assert A.dim == 4


def test_root_of_unity_and_radical_tokens():
    Q = _two_loops()
    r = parse_relation("x*y + E(3)*y*x", Q)
    coeffs = {w: c for c, w in r.terms}
    assert coeffs[("y", "x")] == sympy.exp(2 * sympy.pi * sympy.I / 3)


def test_float_coefficient_still_fails_loudly():
    Q = _two_loops()
    with pytest.raises((ExactnessError, RelationError)):
        parse_relation("0.5*x*y - y*x", Q)


def test_unknown_token_is_a_relation_error_not_a_silent_arrow():
    Q = _two_loops()
    with pytest.raises(RelationError):
        parse_relation("bogus(2)*x", Q)


def test_wrapped_float_coefficient_fails_loudly():
    # A float hidden inside a parenthesised token must not be silently rationalised
    # (controller adjudication: _exact_scalar uses rational=False, matching the
    # loud-exactness contract and fields.complexfield.parse_entry).
    Q = _two_loops()
    with pytest.raises(ExactnessError):
        parse_relation("(0.5)*x*y - y*x", Q)
