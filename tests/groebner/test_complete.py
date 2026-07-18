"""S-polynomials + Buchberger-Mora completion (spec §5, component 3), left-to-right."""
from quiverlab.combinat import Quiver
from quiverlab.fields import QQ
from quiverlab.groebner.order import path_order
from quiverlab.groebner.reduction import rule_from_comb
from quiverlab.groebner.overlap import overlaps
from quiverlab.groebner.complete import s_polynomial, complete


def _two_loops():
    return Quiver([1], {"x": (1, 1), "y": (1, 1)})


def _rules(combs, Q):
    o = path_order(Q)
    return [rule_from_comb(c, o, Q, QQ) for c in combs]


def test_s_polynomial_of_yy_and_yx():
    Q = _two_loops()
    yy, yx = _rules([{("y", "y"): QQ.coerce(1)},
                     {("y", "x"): QQ.coerce(1), ("x", "x"): QQ.coerce(-1)}], Q)
    amb = overlaps(yy, yx, degree_bound=10)[0]   # word yyx
    s = s_polynomial(amb, QQ)
    # tail(yy)=0 so front reduction is 0; back is y*tail(yx)=y*xx=yxx ; S = 0 - yxx
    assert s == {("y", "x", "x"): QQ.coerce(-1)}


def test_completion_fixture3_adds_xxx():
    """Fixture 3: {yy->0, xy->0, yx->xx} completes by ADDING xxx->0."""
    Q = _two_loops()
    o = path_order(Q)
    init = _rules([
        {("y", "y"): QQ.coerce(1)},                                   # yy -> 0
        {("x", "y"): QQ.coerce(1)},                                   # xy -> 0
        {("y", "x"): QQ.coerce(1), ("x", "x"): QQ.coerce(-1)},        # yx -> xx
    ], Q)
    done = complete(init, o, Q, QQ, degree_bound=8)
    leads = sorted(r.lead for r in done)
    assert leads == [("x", "x", "x"), ("x", "y"), ("y", "x"), ("y", "y")]
    # the added rule is a genuine monomial rule xxx -> 0
    xxx = next(r for r in done if r.lead == ("x", "x", "x"))
    assert xxx.tail == ()


def test_completion_fixture6_quantum_ci_adds_nothing():
    """Fixture 6: {xx->0, yy->0, yx->xy} is already confluent."""
    Q = _two_loops()
    o = path_order(Q)
    init = _rules([
        {("x", "x"): QQ.coerce(1)},                                   # xx -> 0
        {("y", "y"): QQ.coerce(1)},                                   # yy -> 0
        {("y", "x"): QQ.coerce(1), ("x", "y"): QQ.coerce(-1)},        # yx -> xy
    ], Q)
    done = complete(init, o, Q, QQ, degree_bound=8)
    assert sorted(r.lead for r in done) == [("x", "x"), ("y", "x"), ("y", "y")]


def test_completion_square_adds_nothing():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    o = path_order(Q)
    cd = rule_from_comb({("c", "d"): QQ.coerce(1), ("a", "b"): QQ.coerce(-1)}, o, Q, QQ)
    done = complete([cd], o, Q, QQ, degree_bound=8)
    assert [r.lead for r in done] == [("c", "d")]
    assert done[0].tail == ((QQ.coerce(1), ("a", "b")),)


def test_completion_terminates_and_is_confluent():
    """After completion, every ambiguity S-polynomial reduces to zero."""
    from quiverlab.groebner.overlap import all_ambiguities
    from quiverlab.groebner.reduction import reduce_comb
    Q = _two_loops()
    o = path_order(Q)
    init = _rules([
        {("y", "y"): QQ.coerce(1)},
        {("x", "y"): QQ.coerce(1)},
        {("y", "x"): QQ.coerce(1), ("x", "x"): QQ.coerce(-1)},
    ], Q)
    done = complete(init, o, Q, QQ, degree_bound=8)
    for amb in all_ambiguities(done, degree_bound=8):
        assert reduce_comb(s_polynomial(amb, QQ), done, o, QQ) == {}
