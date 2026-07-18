"""Finiteness certificate: degree-bound check + irreducible-word enumeration
(spec §3.3, §5 component 3). Loud AdmissibilityError / NotFiniteDimensionalError."""
import pytest

from quiverlab.combinat import Quiver
from quiverlab.fields import QQ
from quiverlab.errors import AdmissibilityError, NotFiniteDimensionalError
from quiverlab.groebner.order import path_order
from quiverlab.groebner.reduction import rule_from_comb
from quiverlab.groebner.complete import complete
from quiverlab.groebner.certificate import (
    default_degree_bound, check_degree_bound, certified_irreducibles,
)


def _two_loops():
    return Quiver([1], {"x": (1, 1), "y": (1, 1)})


def _rules(combs, Q):
    o = path_order(Q)
    return [rule_from_comb(c, o, Q, QQ) for c in combs]


def test_default_bound_covers_completion_growth():
    # Fixture 3 relations have max length 2 -> default 8, and completion reaches L=3 (need 5)
    class _R:
        max_length = 2
    assert default_degree_bound([_R(), _R()]) == 8


def test_check_degree_bound_passes_when_bound_large_enough():
    Q = _two_loops()
    rules = _rules([{("x", "x"): QQ.coerce(1)}], Q)   # L=2, need 2*2-1=3
    check_degree_bound(rules, 8)                       # no raise


def test_check_degree_bound_fixture5_raises_admissibility():
    """Fixture 5: quantum CI q=1 with degree_bound=2. L=2 needs bound >= 3."""
    Q = _two_loops()
    rules = _rules([
        {("x", "x"): QQ.coerce(1)},
        {("y", "y"): QQ.coerce(1)},
        {("y", "x"): QQ.coerce(1), ("x", "y"): QQ.coerce(-1)},
    ], Q)
    with pytest.raises(AdmissibilityError) as exc:
        check_degree_bound(rules, 2)
    assert "degree_bound=2" in str(exc.value)
    assert "at least 3" in str(exc.value)


def test_certified_irreducibles_fixture3():
    Q = _two_loops()
    o = path_order(Q)
    init = _rules([
        {("y", "y"): QQ.coerce(1)},
        {("x", "y"): QQ.coerce(1)},
        {("y", "x"): QQ.coerce(1), ("x", "x"): QQ.coerce(-1)},
    ], Q)
    done = complete(init, o, Q, QQ, degree_bound=8)
    words = certified_irreducibles(Q, done)
    assert words == (("x",), ("y",), ("x", "x"))        # basis {x, y, xx}, dim 4 with vertex


def test_certified_irreducibles_fixture4_raises_not_finite():
    """Fixture 4: two loops, no relations -> infinite, names the loop cycle."""
    Q = Quiver([1], {"a": (1, 1), "b": (1, 1)})
    with pytest.raises(NotFiniteDimensionalError) as exc:
        certified_irreducibles(Q, [])
    msg = str(exc.value)
    assert "infinite" in msg.lower()
    assert ("a" in msg or "b" in msg)                   # the doubled-arrow growth is named


def test_fixture5_and_fixture4_are_distinct_error_types():
    assert AdmissibilityError is not NotFiniteDimensionalError
