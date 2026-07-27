from fractions import Fraction

import pytest
import sympy
from quiverlab.errors import ExactnessError, FieldError, QuiverlabError
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


def test_parse_entry_E_string_regression():
    """Regression: CC.parse_entry("E(3)") used to raise FieldError because E()
    guarded `isinstance(n, int)`, which rejects the sympy.Integer sympify hands it.
    It now parses to the exact primitive cube root."""
    e3 = CC.parse_entry("E(3)")
    assert sympy.simplify(e3 - sympy.exp(2 * sympy.pi * sympy.I / 3)) == 0
    # 1 + w + w^2 = 0 in the field E(3) generates
    dom = CC.make_domain(["E(3)"])
    w = dom.coerce("E(3)")
    assert dom.is_zero(dom.add(dom.add(dom.one(), w), dom.mul(w, w)))


def test_E6_matrix_entry_over_CC():
    """A root-of-unity entry E(6) is a usable exact scalar over CC."""
    dom = CC.make_domain(["E(6)"])
    w = dom.coerce("E(6)")
    w6 = w
    for _ in range(5):
        w6 = dom.mul(w6, w)
    assert dom.eq(w6, dom.one())                       # E(6)^6 = 1
    w3 = dom.mul(dom.mul(w, w), w)
    assert dom.eq(w3, dom.coerce(-1))                  # primitive: E(6)^3 = -1


def test_E_accepts_python_and_sympy_integer():
    assert sympy.simplify(E(5) - sympy.exp(2 * sympy.pi * sympy.I / 5)) == 0
    assert sympy.simplify(E(sympy.Integer(5)) - sympy.exp(2 * sympy.pi * sympy.I / 5)) == 0


def test_E_rejects_non_positive_integers_and_floats():
    # direct calls: floats / non-positive integers refuse with FieldError
    for bad in (3.5, 3.0, sympy.Float("3.0")):
        with pytest.raises(FieldError):
            E(bad)
    for bad in (0, -2):
        with pytest.raises(FieldError):
            E(bad)
    # through the string parser, decimal indices still refuse loudly
    for s in ("E(3.5)", "E(3.0)"):
        with pytest.raises(QuiverlabError):
            CC.parse_entry(s)


def test_cc_repr():
    assert repr(CC) == "CC"
