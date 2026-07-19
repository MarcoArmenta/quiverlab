from quiverlab.resolutions_cs.terms import Chain
from quiverlab.resolutions_cs import trace


def test_chain_records_blocks_and_endpoints():
    ch = Chain(word=("x", "x", "x"), blocks=(("x",), ("x",), ("x",)), o=1, t=1, degree=3)
    assert ch.n_blocks == 3 and ch.degree == 3 and ch.word == ("x", "x", "x")


def test_chain_hashable_and_equal():
    a = Chain(("x", "x"), (("x",), ("x",)), 1, 1, 2)
    b = Chain(("x", "x"), (("x",), ("x",)), 1, 1, 2)
    assert a == b and hash(a) == hash(b) and {a, b} == {a}


def test_trace_dataclasses_are_inert():
    ev = trace.DifferentialEvent(degree=2, chain=("c", "d"), terms=[(1, (), ("c",), ("d",))])
    assert ev.degree == 2 and not hasattr(ev, "render")
