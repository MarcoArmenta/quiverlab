"""Plan-04 interface freshness gate: the REAL frozen Plan-03 reduction-system API
(docs/plans/2026-07-18-plan-03-groebner.md). STOP on any drift.

Field-object contract for Task 7's `_DomainField` shim (Step-5 probe, verified
against the committed src/quiverlab/groebner/system.py):
    grep -nE 'field\\.[A-Za-z_]+' src/quiverlab/groebner/system.py  ==>
        field.parse_entry(0)      # line 64
        field.parse_entry(1)      # line 64
        field.parse_entry(c)      # line 67  (c is a relation coefficient)
        field.make_domain(raw)    # line 68
    (the only other `field.`-adjacent hit is `dom.coerce(field.parse_entry(c))`
     on line 74, i.e. still just parse_entry.)
The complete `field` surface build_reduction_system touches is EXACTLY:
    * parse_entry(x)   -- x in {0, 1, coefficient}; coefficients arrive as
                          fractions.Fraction, and Plan-03 does
                          dom.coerce(field.parse_entry(c)), so the shim's
                          parse_entry must return something A.domain.coerce
                          accepts -> pass the Fraction straight through.
    * make_domain(entries) -- returns the coefficient Domain; the shim returns
                          the existing A.domain.
No field.name / field.characteristic / any other method is called. Task 7's
`_DomainField` shim must implement exactly {parse_entry, make_domain} -- nothing
more, nothing less.
"""
from dataclasses import fields as dc_fields
import pytest

groebner = pytest.importorskip("quiverlab.groebner",
                               reason="Plan 03 (groebner) must be committed before Plan 04.")


def test_frozen_symbols():
    for n in ("ReductionSystem", "ReductionRule", "Ambiguity", "PathOrder",
              "build_reduction_system"):
        assert hasattr(groebner, n), f"groebner.{n} missing (Plan-03 drift)"


def test_reduction_system_shape():
    names = {f.name for f in dc_fields(groebner.ReductionSystem)}
    assert {"quiver", "domain", "order", "rules", "irreducibles",
            "degree_bound", "is_confluent"} <= names
    for m in ("leading_words", "reduce", "normal_form", "ambiguities"):
        assert callable(getattr(groebner.ReductionSystem, m))


def test_rule_and_ambiguity_shape(square_rs):
    rule = square_rs.rules[0]
    assert rule.lead == ("c", "d")                                   # words are ARROW-NAME tuples
    assert rule.tail == ((square_rs.domain.coerce(1), ("a", "b")),)  # tuple of (coeff, word), NOT a dict
    assert (rule.source, rule.target) == (1, 4)
    for amb in square_rs.ambiguities():
        assert amb.kind in ("overlap", "inclusion")
        assert isinstance(amb.word, tuple) and amb.left is not None and amb.right is not None


def test_normal_form_reduces_tip(square_rs):
    one = square_rs.domain.coerce(1)
    assert square_rs.normal_form(("c", "d")) == {("a", "b"): one}    # cd -> ab (a full dict)


def test_no_algebra_reduction_system_accessor():
    from quiverlab import Quiver, CC
    A = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=CC)
    assert not hasattr(A, "reduction_system")                        # CS gets the RS via build_reduction_system
