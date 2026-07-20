"""QQi: first-class public exact Gaussian-rationals field Q(i) = QQ_I.

Field-unit contract (domain arithmetic pinned to Q(i)) plus an end-to-end
integration check that the qci `i` algebra built over QQi agrees with the same
algebra over CC (both compute in the exact QQ_I domain)."""
from fractions import Fraction

import pytest
import sympy

from quiverlab.errors import ExactnessError, FieldError
from quiverlab.fields import QQi, CC


# -- field-spec surface -----------------------------------------------------
def test_name_char_and_repr():
    assert QQi.name == "QQ(i)"
    assert QQi.characteristic == 0
    assert repr(QQi) == "QQ(i)"


# -- domain arithmetic: exact Q(i) ------------------------------------------
def test_i_squared_is_minus_one():
    d = QQi.make_domain(["i"])
    assert d.characteristic == 0
    i = d.coerce("i")
    assert d.eq(d.mul(i, i), d.coerce(-1))          # i^2 = -1, exactly


def test_exact_inverse():
    d = QQi.make_domain(["i"])
    x = d.coerce("1+i")
    xinv = d.inv(x)
    # inv(1+i) = (1-i)/2, and x*x^{-1} = 1, all exact (no float anywhere)
    assert d.eq(xinv, d.coerce("1/2 - 1/2*i"))
    assert d.eq(d.mul(x, xinv), d.one())


def test_add_sub_neg_eq_is_zero():
    d = QQi.make_domain(["i"])
    a = d.coerce("1+2*i")
    b = d.coerce("3-i")
    assert d.eq(d.add(a, b), d.coerce("4+i"))
    assert d.eq(d.sub(a, b), d.coerce("-2+3*i"))
    assert d.eq(d.neg(a), d.coerce("-1-2*i"))
    assert d.is_zero(d.sub(a, a))
    assert not d.is_zero(a)


def test_rational_and_int_and_fraction_entries():
    d = QQi.make_domain([1, Fraction(1, 2), "2/3"])
    half = d.coerce(Fraction(1, 2))
    assert d.eq(d.add(half, half), d.one())
    assert d.eq(d.mul(d.coerce(2), d.coerce("2/3")), d.coerce("4/3"))


def test_domain_is_pinned_to_QQ_I_not_broader_closure():
    # QQi is Q(i), never the broader algebraic closure: sqrt(2) is out of field.
    d = QQi.make_domain(["sqrt(2)", "i"])          # entries ignored for pinning
    assert str(d.sdom) == "QQ_I"
    with pytest.raises(FieldError):
        d.coerce("sqrt(2)")


# -- exactness gate ---------------------------------------------------------
def test_float_and_complex_rejected():
    with pytest.raises(ExactnessError):
        QQi.parse_entry(0.5)
    with pytest.raises(ExactnessError):
        QQi.parse_entry("0.5")
    with pytest.raises(ExactnessError):
        QQi.parse_entry(1 + 2j)
    with pytest.raises(ExactnessError):
        QQi.make_domain([0.5])


def test_parse_entry_accepts_exact_qi_literals():
    for lit in ("i", "1+2*i", -1, Fraction(1, 3)):
        e = QQi.parse_entry(lit)
        assert isinstance(e, sympy.Basic)


# -- algebra builds with structure constants in exact Q(i) ------------------
def test_algebra_builds_over_QQi():
    from quiverlab import Quiver
    pytest.importorskip("quiverlab.groebner")
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "y*x - i*x*y"], field=QQi)
    assert A.domain.characteristic == 0
    assert str(A.domain.sdom) == "QQ_I"
    # a structure constant carries the exact i coefficient
    flat = [c for row in A.T for vec in row for c in vec]
    assert any(A.domain.to_str(c) not in ("0", "1", "-1") and "I" in A.domain.to_str(c)
               for c in flat)


# -- integration: QQi end-to-end agrees with CC on the qci `i` algebra ------
def _qci(field):
    from quiverlab import Quiver
    return Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "y*x - i*x*y"], field=field)


def test_qqi_hh_dims_match_cc():
    from quiverlab.hochschild.bar import (
        hochschild_cohomology_dims, hochschild_homology_dims)
    pytest.importorskip("quiverlab.groebner")
    A_cc = _qci(CC)
    A_qqi = _qci(QQi)
    assert (hochschild_cohomology_dims(A_qqi, 3).dims
            == hochschild_cohomology_dims(A_cc, 3).dims)
    assert (hochschild_homology_dims(A_qqi, 3).dims
            == hochschild_homology_dims(A_cc, 3).dims)


def test_qqi_cs_equals_bar():
    pytest.importorskip("quiverlab.groebner")
    from quiverlab.resolutions_cs.homology import (
        cs_cohomology_dims, cs_homology_dims)
    from quiverlab.hochschild.bar import (
        hochschild_cohomology_dims, hochschild_homology_dims)
    A = _qci(QQi)
    assert cs_homology_dims(A, 3).dims == hochschild_homology_dims(A, 3).dims
    assert cs_cohomology_dims(A, 3).dims == hochschild_cohomology_dims(A, 3).dims
