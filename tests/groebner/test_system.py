"""ReductionSystem: the object Chouhy-Solotar (Plan 04) consumes verbatim."""
import pytest
from dataclasses import fields as dc_fields

from quiverlab.combinat import Quiver
from quiverlab.fields import CC, QQ
from quiverlab.errors import AdmissibilityError, NotFiniteDimensionalError, RelationError
from quiverlab.groebner import ReductionSystem, ReductionRule, build_reduction_system
from quiverlab.groebner.reduction import ReductionRule as _RR
from quiverlab.groebner.overlap import Ambiguity


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})


def test_frozen_public_shape():
    """Plan 04 depends on these field/method names -- freeze them."""
    names = {f.name for f in dc_fields(ReductionSystem)}
    assert {"quiver", "domain", "order", "rules", "irreducibles",
            "degree_bound", "is_confluent"} <= names
    for meth in ("leading_words", "reduce", "normal_form", "ambiguities"):
        assert callable(getattr(ReductionSystem, meth))
    assert ReductionRule is _RR


def test_build_square_system():
    rs = build_reduction_system(_square(), ["a*b - c*d"], CC)
    assert rs.leading_words() == (("c", "d"),)
    assert rs.irreducibles == (("a",), ("b",), ("c",), ("d",), ("a", "b"))
    assert rs.is_confluent is True


def test_reduce_and_normal_form_use_the_rules():
    rs = build_reduction_system(_square(), ["a*b - c*d"], CC)
    one = rs.domain.coerce(1)
    # c*d reduces to a*b
    assert rs.normal_form(("c", "d")) == {("a", "b"): one}
    assert rs.reduce({("c", "d"): one}) == {("a", "b"): one}


def test_ambiguities_of_completed_system_all_resolve():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    rs = build_reduction_system(Q, ["y^2", "x*y", "y*x - x*x"], QQ)   # Fixture 3
    assert sorted(r.lead for r in rs.rules) == [
        ("x", "x", "x"), ("x", "y"), ("y", "x"), ("y", "y")]
    for amb in rs.ambiguities():
        assert isinstance(amb, Ambiguity)
        assert rs.reduce(_s(amb, rs)) == {}


def _s(amb, rs):
    from quiverlab.groebner.complete import s_polynomial
    return s_polynomial(amb, rs.domain)


def test_build_rejects_length_one_relation():
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (1, 3)})
    with pytest.raises(AdmissibilityError):
        build_reduction_system(Q, ["a*b - c"], CC)   # c has length 1


def test_build_raises_admissibility_when_bound_too_small():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    with pytest.raises(AdmissibilityError):
        build_reduction_system(Q, ["x^2", "y^2", "y*x - x*y"], CC, degree_bound=2)  # Fixture 5


def test_build_raises_not_finite_dimensional():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    with pytest.raises(NotFiniteDimensionalError):
        build_reduction_system(Q, ["x*y - y*x"], CC)   # commutative k[x,y], infinite


def test_build_rejects_relation_vanishing_in_field():
    from quiverlab.fields import GF
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    # 2*x*y is a single monomial with coefficient 2; over GF(2) it vanishes.
    with pytest.raises(RelationError):
        build_reduction_system(Q, ["2*x*y"], GF(2))
