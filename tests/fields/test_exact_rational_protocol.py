"""Exact-scalar readers accept exact rational objects by PROTOCOL, not name-list.

Regression (Marco, 2026-08-02): a dim-7 quiver algebra over CC died with
``FieldError: cannot read PythonMPQ MPQ(1,1) as an exact scalar`` when the
CS-native cup/cap products (and the HH representative capture) fed sympy's own
QQ-domain elements -- ``sympy.external.pythonmpq.PythonMPQ`` in a pure install,
``gmpy2`` ``mpq``/``mpz`` when the C backend is present -- back through the
exact-scalar readers. Those values are EXACT rationals; refusing them as "not an
exact scalar" is a FALSE refusal that blamed the user for floats they never typed.

Both raisers (``ComplexField.parse_entry`` for CC, ``parse_rational`` for the
generic/GF(p) Domain) now admit anything with an integer ``.numerator`` /
``.denominator``. Genuine floats/complex/bool still refuse loudly with the hint.
"""
from decimal import Decimal
from fractions import Fraction

import pytest
import sympy

from quiverlab.errors import ExactnessError, FieldError
from quiverlab.fields import CC, QQ
from quiverlab.fields.domain import exact_rational, parse_rational

try:
    from sympy.external.pythonmpq import PythonMPQ
except Exception:                                    # pragma: no cover - sympy internal move
    PythonMPQ = None


# --------------------------------------------------------------------------- #
# The shared protocol helper
# --------------------------------------------------------------------------- #
def test_exact_rational_admits_exact_types():
    assert exact_rational(3) == (3, 1)
    assert exact_rational(Fraction(2, 3)) == (2, 3)
    assert exact_rational(sympy.Rational(3, 4)) == (3, 4)
    assert exact_rational(sympy.Integer(7)) == (7, 1)
    if PythonMPQ is not None:
        assert exact_rational(PythonMPQ(1, 1)) == (1, 1)
        assert exact_rational(PythonMPQ(-6, 8)) == (-3, 4)   # normalized by the type


def test_exact_rational_misses_inexact_and_nonrational():
    # float / complex / Decimal carry no integer num/den pair -> None (miss).
    for bad in (0.5, 0.0, 1 + 2j, Decimal("1.5"), Decimal("2"), sympy.Float("0.5"),
                "1/2", None, object()):
        assert exact_rational(bad) is None
    # A truncation trap: a protocol object with a NON-integral numerator must miss,
    # never be silently int()-truncated to an exact value.
    class Bad:
        numerator = Fraction(3, 2)
        denominator = 1
    assert exact_rational(Bad()) is None


# --------------------------------------------------------------------------- #
# CC reader (ComplexField.parse_entry, reached via SympyExactDomain.coerce)
# --------------------------------------------------------------------------- #
def test_cc_reader_accepts_exact_rational_protocol():
    assert CC.parse_entry(Fraction(3, 4)) == sympy.Rational(3, 4)
    assert CC.parse_entry(5) == sympy.Integer(5)
    assert CC.parse_entry(sympy.Rational(1, 2)) == sympy.Rational(1, 2)
    if PythonMPQ is not None:
        assert CC.parse_entry(PythonMPQ(1, 1)) == sympy.Integer(1)
        assert CC.parse_entry(PythonMPQ(3, 4)) == sympy.Rational(3, 4)


def test_cc_domain_coerce_roundtrips_its_own_elements():
    """The concrete QQ domain's OWN element representation (mpq/PythonMPQ) must
    re-coerce to itself -- this is exactly the loop the products/HH-reps code took."""
    dom = CC.make_domain([1, Fraction(1, 2)])
    native = dom.coerce("3/4")                       # a raw domain element (mpq/PythonMPQ)
    assert dom.eq(dom.coerce(native), native)        # feeding it back in is a no-op
    assert type(native).__name__ in ("mpq", "PythonMPQ", "Rational", "Half")


def test_cc_reader_still_refuses_floats_and_bool_with_hint():
    for bad in (0.5, 1 + 2j, sympy.Float("0.5")):
        with pytest.raises((ExactnessError, FieldError)):
            CC.parse_entry(bad)
    # bool is a numbers.Rational but must NOT slip through as 0/1 (house style).
    with pytest.raises(ExactnessError):
        CC.parse_entry(True)
    # the genuine-float refusal keeps its exact-only hint text.
    with pytest.raises(ExactnessError) as ei:
        CC.parse_entry(0.5)
    assert "exact-only" in str(ei.value)


# --------------------------------------------------------------------------- #
# Generic / GF(p) reader (parse_rational, reached via QQ.coerce & GF(p).coerce)
# --------------------------------------------------------------------------- #
def test_generic_reader_accepts_exact_rational_protocol():
    assert parse_rational(Fraction(3, 4)) == Fraction(3, 4)
    assert parse_rational(7) == Fraction(7)
    assert parse_rational(sympy.Rational(1, 2)) == Fraction(1, 2)
    assert parse_rational(sympy.Integer(5)) == Fraction(5)
    assert QQ.coerce(sympy.Rational(2, 3)) == Fraction(2, 3)
    if PythonMPQ is not None:
        assert parse_rational(PythonMPQ(3, 4)) == Fraction(3, 4)
        assert QQ.coerce(PythonMPQ(1, 1)) == Fraction(1, 1)


def test_generic_reader_still_refuses_floats_and_bool_with_hint():
    for bad in (0.5, 1 + 2j):
        with pytest.raises(ExactnessError):
            parse_rational(bad)
    with pytest.raises(ExactnessError):
        parse_rational(True)
    with pytest.raises(ExactnessError) as ei:
        parse_rational(0.5)
    assert "exact-only" in str(ei.value)


# --------------------------------------------------------------------------- #
# gmpy2 backend (only when the C library is importable)
# --------------------------------------------------------------------------- #
def test_gmpy2_mpq_and_mpz_coerce_exactly_in_both_readers():
    gmpy2 = pytest.importorskip("gmpy2")
    q, z = gmpy2.mpq(3, 4), gmpy2.mpz(5)
    assert exact_rational(q) == (3, 4)
    assert exact_rational(z) == (5, 1)
    # CC reader
    assert CC.parse_entry(q) == sympy.Rational(3, 4)
    assert CC.parse_entry(z) == sympy.Integer(5)
    # generic reader
    assert parse_rational(q) == Fraction(3, 4)
    assert parse_rational(z) == Fraction(5)


# --------------------------------------------------------------------------- #
# End-to-end (reader-independent): the products path that actually broke
# --------------------------------------------------------------------------- #
def test_cup_and_cap_products_over_cc_no_false_refusal():
    """The Plan-35 CS-native cup/cap over CC solve their class coordinates through
    ``dom.coerce`` on sympy-QQ matrix entries (PythonMPQ/mpq). This is the exact
    site of the false refusal; it must now compute."""
    from quiverlab import Quiver

    A = Quiver([1, 2], {"a": (1, 2), "b": (2, 2)}).algebra(["b*b*b"], field=CC)
    cup = A.cup_products(2)
    cap = A.cap_products(2)
    assert cup.tables and cap.tables            # non-empty product tables, no FieldError
