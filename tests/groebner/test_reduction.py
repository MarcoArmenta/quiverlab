"""Reduction rules and normal forms over a Domain, left-to-right (spec §5 c.3)."""
from fractions import Fraction

from quiverlab.combinat import Quiver
from quiverlab.fields import QQ
from quiverlab.groebner.order import path_order
from quiverlab.groebner.reduction import (
    ReductionRule, rule_from_comb, first_factor, reduce_comb, lc_sub,
)


def _two_loops():
    return Quiver([1], {"x": (1, 1), "y": (1, 1)})


def _c(dom, n):
    return dom.coerce(n)


def test_rule_from_monomial_comb_has_empty_tail():
    Q = _two_loops()
    o = path_order(Q)
    dom = QQ
    r = rule_from_comb({("x", "x", "x"): _c(dom, 1)}, o, Q, dom)
    assert r.lead == ("x", "x", "x")
    assert r.tail == ()
    assert (r.source, r.target) == (1, 1)


def test_rule_from_binomial_picks_leading_and_normalizes():
    Q = _two_loops()
    o = path_order(Q)
    dom = QQ
    # relation y*x - x*x  ==  {yx: 1, xx: -1};  x<y so leading is yx, tail = {xx: 1}
    comb = {("y", "x"): _c(dom, 1), ("x", "x"): _c(dom, -1)}
    r = rule_from_comb(comb, o, Q, dom)
    assert r.lead == ("y", "x")
    assert r.tail == ((_c(dom, 1), ("x", "x")),)


def test_rule_normalizes_leading_coefficient_to_one():
    Q = _two_loops()
    o = path_order(Q)
    dom = QQ
    # 2*yx - 3*xx : leading yx (coeff 2), tail = (3/2) xx
    comb = {("y", "x"): _c(dom, 2), ("x", "x"): _c(dom, -3)}
    r = rule_from_comb(comb, o, Q, dom)
    assert r.lead == ("y", "x")
    assert r.tail == ((Fraction(3, 2), ("x", "x")),)


def test_first_factor_finds_leftmost_occurrence():
    Q = _two_loops()
    dom = QQ
    r = rule_from_comb({("y", "x"): _c(dom, 1), ("x", "x"): _c(dom, -1)}, path_order(Q), Q, dom)
    hit = first_factor(("y", "x", "x"), [r])
    assert hit is not None
    rule, i = hit
    assert (rule.lead, i) == (("y", "x"), 0)
    assert first_factor(("x", "x"), [r]) is None


def test_reduce_single_monomial_rule():
    Q = _two_loops()
    o = path_order(Q)
    dom = QQ
    xxx = rule_from_comb({("x", "x", "x"): _c(dom, 1)}, o, Q, dom)  # xxx -> 0
    assert reduce_comb({("x", "x", "x", "x"): _c(dom, 1)}, [xxx], o, dom) == {}
    assert reduce_comb({("x", "x"): _c(dom, 1)}, [xxx], o, dom) == {("x", "x"): _c(dom, 1)}


def test_reduce_chains_rules_to_normal_form():
    Q = _two_loops()
    o = path_order(Q)
    dom = QQ
    yx = rule_from_comb({("y", "x"): _c(dom, 1), ("x", "x"): _c(dom, -1)}, o, Q, dom)  # yx -> xx
    xxx = rule_from_comb({("x", "x", "x"): _c(dom, 1)}, o, Q, dom)                     # xxx -> 0
    # y*x*x  --yx-->  x*x*x  --xxx-->  0
    assert reduce_comb({("y", "x", "x"): _c(dom, 1)}, [yx, xxx], o, dom) == {}


def test_reduce_emits_trace_steps_when_requested():
    Q = _two_loops()
    o = path_order(Q)
    dom = QQ
    yx = rule_from_comb({("y", "x"): _c(dom, 1), ("x", "x"): _c(dom, -1)}, o, Q, dom)
    trace = []
    reduce_comb({("y", "x"): _c(dom, 1)}, [yx], o, dom, trace=trace)
    assert len(trace) == 1
    assert trace[0].rule_lead == ("y", "x")


def test_lc_sub_cancels():
    dom = QQ
    a = {("a", "b"): _c(dom, 1)}
    b = {("a", "b"): _c(dom, 1)}
    assert lc_sub(a, b, dom) == {}
