"""Length-lexicographic admissible path order (spec §5, component 3)."""
import itertools

from quiverlab.combinat import Quiver
from quiverlab.groebner import PathOrder, path_order, Dispatch, ReductionStep


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})


def test_arrow_ranks_follow_insertion_order():
    o = path_order(_square())
    assert o.arrow_index == {"a": 0, "b": 1, "c": 2, "d": 3}


def test_shorter_word_is_smaller():
    o = path_order(_square())
    assert o.compare(("a",), ("a", "b")) == -1
    assert o.compare(("a", "b"), ("a",)) == 1


def test_lex_breaks_ties_by_arrow_rank():
    o = path_order(_square())
    # equal length: a<c so a*b < c*d
    assert o.compare(("a", "b"), ("c", "d")) == -1
    assert o.compare(("c", "d"), ("a", "b")) == 1
    assert o.compare(("a", "b"), ("a", "b")) == 0


def test_leading_is_the_maximum_word():
    o = path_order(Quiver([1], {"x": (1, 1), "y": (1, 1)}))
    # y*x vs x*x with x<y: leading is y*x
    assert o.leading({("y", "x"): 1, ("x", "x"): 1}) == ("y", "x")
    # longer beats shorter regardless of lex
    assert o.leading({("x", "x", "x"): 1, ("y", "x"): 1}) == ("x", "x", "x")
    assert o.leading({}) is None


def test_two_sided_multiplicative_compatibility():
    """u < v  ==>  w u x < w v x for all composable contexts (checked combinatorially)."""
    o = path_order(Quiver([1], {"x": (1, 1), "y": (1, 1)}))
    letters = ["x", "y"]
    words = [tuple(p) for n in range(1, 4) for p in itertools.product(letters, repeat=n)]
    for u in words:
        for v in words:
            if o.compare(u, v) >= 0:
                continue
            for w in [()] + [(l,) for l in letters]:
                for x in [()] + [(l,) for l in letters]:
                    assert o.compare(w + u + x, w + v + x) < 0


def test_events_are_plain_dataclasses():
    d = Dispatch(route="groebner", reason="non-monomial", n_relations=1)
    assert (d.route, d.n_relations) == ("groebner", 1)
    s = ReductionStep(word=("c", "d"), rule_lead=("c", "d"),
                      before={("c", "d"): 1}, after={("a", "b"): 1})
    assert s.rule_lead == ("c", "d")
