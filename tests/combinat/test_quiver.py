import pytest
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import RelationError


def _q():
    return Quiver(vertices=[1, 2, 3],
                  arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})


def test_basic_accessors():
    Q = _q()
    assert Q.source("a") == 1 and Q.target("a") == 2
    assert Q.word_source(("a", "b")) == 1 and Q.word_target(("a", "b")) == 3
    assert Q.compose_ok(("a", "b")) and not Q.compose_ok(("b", "a"))
    assert Q.is_acyclic()


def test_loops_and_cycles():
    L = Quiver(vertices=[1], arrows={"x": (1, 1)})
    assert not L.is_acyclic()
    assert L.compose_ok(("x", "x", "x"))


def test_validation():
    with pytest.raises(RelationError):
        Quiver(vertices=[1], arrows={"a": (1, 2)})       # endpoint not a vertex
    with pytest.raises(RelationError):
        Quiver(vertices=[1], arrows={"a*b": (1, 1)})     # bad arrow name
    with pytest.raises(RelationError):
        Quiver(vertices=[1, 1], arrows={})               # duplicate vertex


def test_repr_shows_arrows():
    s = repr(_q())
    assert "1 --a--> 2" in s
