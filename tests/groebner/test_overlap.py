"""Overlap/inclusion ambiguities in the LEFT-TO-RIGHT convention (spec §5 c.3):
b is a SUFFIX of the first lead and a PREFIX of the second."""
from quiverlab.combinat import Quiver
from quiverlab.fields import QQ
from quiverlab.groebner.order import path_order
from quiverlab.groebner.reduction import rule_from_comb
from quiverlab.groebner.overlap import overlaps, inclusions, all_ambiguities


def _two_loops():
    return Quiver([1], {"x": (1, 1), "y": (1, 1)})


def _rule(comb):
    Q = _two_loops()
    return rule_from_comb(comb, path_order(Q), Q, QQ)


def test_self_overlap_of_xyx():
    r = _rule({("x", "y", "x"): QQ.coerce(1), ("y",): QQ.coerce(1)})  # lead xyx (len3 > y)
    amb = overlaps(r, r, degree_bound=10)
    words = sorted(a.word for a in amb)
    # xyx suffix "x" = prefix "x" -> a=xy, c=yx -> word xyxyx (only overlap)
    assert words == [("x", "y", "x", "y", "x")]
    a0 = amb[0]
    assert a0.kind == "overlap"
    assert (a0.a, a0.c) == (("x", "y"), ("y", "x"))


def test_overlap_yy_with_yx():
    yy = _rule({("y", "y"): QQ.coerce(1)})               # yy -> 0
    yx = _rule({("y", "x"): QQ.coerce(1), ("x", "x"): QQ.coerce(-1)})  # yx -> xx
    amb = overlaps(yy, yx, degree_bound=10)
    assert [a.word for a in amb] == [("y", "y", "x")]
    assert (amb[0].a, amb[0].c) == (("y",), ("x",))


def test_no_overlap_when_suffix_prefix_disagree():
    cd = _rule({("y", "x"): QQ.coerce(1)})   # lead yx: suffix "x" != prefix "y"
    assert overlaps(cd, cd, degree_bound=10) == []


def test_degree_bound_filters_long_ambiguities():
    r = _rule({("x", "y", "x"): QQ.coerce(1), ("y",): QQ.coerce(1)})
    assert overlaps(r, r, degree_bound=4) == []   # xyxyx has length 5 > 4
    assert len(overlaps(r, r, degree_bound=5)) == 1


def test_inclusion_detects_proper_factor():
    inner = _rule({("x", "x"): QQ.coerce(1)})               # xx -> 0
    outer = _rule({("x", "x", "x"): QQ.coerce(1)})          # xxx -> 0
    inc = inclusions(outer, inner, degree_bound=10)
    # xx occurs in xxx at positions 0 and 1
    assert sorted((a.a, a.c) for a in inc) == [((), ("x",)), (("x",), ())]
    assert all(a.kind == "inclusion" and a.word == ("x", "x", "x") for a in inc)


def test_all_ambiguities_covers_ordered_pairs_including_self():
    yy = _rule({("y", "y"): QQ.coerce(1)})
    yx = _rule({("y", "x"): QQ.coerce(1), ("x", "x"): QQ.coerce(-1)})
    ambs = all_ambiguities([yy, yx], degree_bound=10)
    got = {(a.left.lead, a.right.lead, a.word) for a in ambs}
    # (yy,yy)->yyy ; (yy,yx)->yyx ; (yx,yy)? suffix x != prefix y -> none
    assert (("y", "y"), ("y", "y"), ("y", "y", "y")) in got
    assert (("y", "y"), ("y", "x"), ("y", "y", "x")) in got
