from fractions import Fraction

import pytest
import sympy
from quiverlab.errors import ExactnessError, FieldError
from quiverlab.fields import CC, E


def test_rational_entries_fast_path():
    dom = CC.make_domain([1, Fraction(1, 2), "2/3"])
    a = dom.coerce("2/3")
    assert dom.characteristic == 0
    assert dom.eq(dom.add(a, a), dom.coerce("4/3"))
    assert dom.is_zero(dom.sub(a, a))


def test_i_and_sqrt2():
    dom = CC.make_domain(["i", "sqrt(2)"])
    i = dom.coerce("i")
    s = dom.coerce("sqrt(2)")
    assert dom.eq(dom.mul(i, i), dom.coerce(-1))
    assert dom.eq(dom.mul(s, s), dom.coerce(2))
    assert dom.eq(dom.mul(i, dom.inv(i)), dom.one())


def test_root_of_unity():
    dom = CC.make_domain([E(3)])
    w = dom.coerce(E(3))
    w3 = dom.mul(dom.mul(w, w), w)
    assert dom.eq(w3, dom.one())
    # 1 + w + w^2 = 0
    assert dom.is_zero(dom.add(dom.add(dom.one(), w), dom.mul(w, w)))


def test_floats_fail_loudly_everywhere():
    with pytest.raises(ExactnessError):
        CC.parse_entry(0.5)
    with pytest.raises(ExactnessError):
        CC.parse_entry("0.5")
    with pytest.raises(ExactnessError):
        CC.parse_entry(sympy.Float("0.5"))
    with pytest.raises(ExactnessError):
        CC.parse_entry(1 + 2j)


def test_non_number_rejected():
    with pytest.raises(FieldError):
        CC.parse_entry("x + 1")


def test_transcendental_entries_fail_loudly():
    with pytest.raises(FieldError):
        CC.make_domain(["pi"])
    with pytest.raises(FieldError):
        CC.make_domain(["exp(2)"])


def test_bare_E_fails_loudly():
    with pytest.raises(FieldError):
        CC.parse_entry("E")


def test_cc_repr():
    assert repr(CC) == "CC"
