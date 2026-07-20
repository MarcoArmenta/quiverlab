"""Layered layout in exact int/Fraction coordinates (spec §3.7). We golden-test
the LAYOUT DATA (coordinates + routing), never pixels."""
from fractions import Fraction
from numbers import Integral

from quiverlab import Quiver
from quiverlab.viz.layout import layout, layer, LayoutData, EdgeRoute, LoopRoute


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})


def test_layer_depths_by_longest_path():
    assert layer(_square()) == {1: 0, 2: 1, 3: 1, 4: 2}


def test_layer_handles_loops_single_column():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    assert layer(Q) == {1: 0}


def test_square_positions_are_exact():
    L = layout(_square(), relations=["a*b - c*d"])
    assert L.positions == {
        1: (0, Fraction(0)),
        2: (1, Fraction(1, 2)),
        3: (1, Fraction(-1, 2)),
        4: (2, Fraction(0)),
    }
    assert L.relations == ("a*b - c*d",)
    assert L.loops == ()
    kinds = {e.name: e.kind for e in L.edges}
    assert kinds == {"a": "straight", "b": "straight", "c": "straight", "d": "straight"}


def test_all_coordinates_are_int_or_fraction_never_float():
    """Compensating test for the AST-gate decision: no float can leak out of viz."""
    for Q, rels in [
        (_square(), ["a*b - c*d"]),
        (Quiver([1], {"x": (1, 1), "y": (1, 1)}), []),
        (Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}), []),
    ]:
        L = layout(Q, relations=rels)
        for (x, y) in L.positions.values():
            assert isinstance(x, Integral) and isinstance(y, Fraction)
            assert not isinstance(x, float) and not isinstance(y, float)
        for e in L.edges:
            assert isinstance(e.bend, Fraction) and not isinstance(e.bend, float)
        for lp in L.loops:
            assert isinstance(lp.angle_deg, Integral)


def test_two_loops_become_loop_routes_at_integer_angles():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    L = layout(Q)
    got = sorted((lp.name, lp.at, lp.angle_deg) for lp in L.loops)
    assert got == [("x", 1, 90), ("y", 1, 30)]
    assert L.edges == ()  # loops are not edges


def test_parallel_arrows_get_symmetric_fraction_bends():
    Q = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)})
    L = layout(Q)
    bends = {e.name: e.bend for e in L.edges}
    assert bends == {"a": Fraction(1, 4), "b": Fraction(-1, 4)}
    assert all(e.kind == "parallel" for e in L.edges)
