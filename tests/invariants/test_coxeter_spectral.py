"""The assembled Coxeter spectral report: verdicts, order, cross-checked against
the shipped spectral_radius/mahler_measure."""
import sympy as sp
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.invariants.spectral import spectral_radius, mahler_measure


def _kron(m, field=QQ):
    return Quiver([1, 2], {f"a{i}": (1, 2) for i in range(m)}).algebra(field=field)


@pytest.mark.oracle_literature
def test_ka2_is_phi3_order_3():
    r = linear_path_algebra(2, field=QQ).coxeter_spectral()
    assert r["cyclotomic"] and r["quasi_unipotent"]
    assert [f["cyclotomic_index"] for f in r["factorization"]] == [3]
    assert r["coxeter_order"] == 3
    assert r["outside_unit_circle_count"] == 0
    assert r["spectral_radius"]["value"] == "1" and r["mahler_measure"]["value"] == "1"


@pytest.mark.oracle_literature
def test_two_kronecker_cyclotomic_type_but_infinite_order():
    r = _kron(2).coxeter_spectral()
    assert r["cyclotomic"] and r["quasi_unipotent"]               # (t-1)^2 = Φ_1^2
    assert r["coxeter_order"] is None                            # Jordan block (affine)
    assert "Jordan" in r["coxeter_order_reason"] or "infinite" in r["coxeter_order_reason"]


@pytest.mark.oracle_literature
def test_three_kronecker_wild_certified():
    r = _kron(3).coxeter_spectral()
    assert not r["cyclotomic"] and r["coxeter_order"] is None
    assert r["outside_unit_circle_count"] == 1
    assert r["spectral_radius"]["minpoly"] == [1, -7, 1]         # x^2 - (m^2-2)x + 1
    assert sp.simplify(sp.sympify(r["spectral_radius"]["value"]) - (7 + 3*sp.sqrt(5))/2) == 0


@pytest.mark.oracle_crossengine
def test_report_matches_shipped_primitives():
    for A in (linear_path_algebra(4, field=QQ), _kron(3),
              _kron(5, field=GF(32003))):
        chi = A.coxeter_polynomial().as_expr()
        r = A.coxeter_spectral()
        if "value" in r["spectral_radius"]:
            assert sp.simplify(sp.sympify(r["spectral_radius"]["value"])
                               - spectral_radius(chi)) == 0
            assert sp.simplify(sp.sympify(r["mahler_measure"]["value"])
                               - mahler_measure(chi)) == 0


@pytest.mark.oracle_selfcert
def test_singular_cartan_refuses_per_field():
    # DEVIATION from the plan's k[x]/(x^2) example: that algebra has a NON-singular
    # Cartan (C = [[2]], det 2), so Phi = [-1] and chi = t+1 -- it computes fine and
    # is NOT a refusal.  A genuine singular-Cartan witness is the 2-cycle with
    # rad^2 = 0 (C = [[1,1],[1,1]], det 0), for which coxeter_polynomial refuses.
    A = Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=GF(5))
    r = A.coxeter_spectral()
    assert "error" in r["coxeter_polynomial"]                    # singular Cartan, captured
    assert "error" in r["spectral_radius"] and "error" in r["mahler_measure"]
