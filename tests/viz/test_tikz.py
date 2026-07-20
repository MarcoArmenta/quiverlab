"""A.tikz(): TikZ from the SAME layout coordinates as draw() (spec §3.7)."""
import pathlib
from fractions import Fraction

from quiverlab import Quiver, CC
from quiverlab.viz.tikz import tikz_quiver
from quiverlab.viz.layout import layout

GOLDEN = pathlib.Path(__file__).parent / "golden" / "ka2.tikz"


def test_ka2_matches_golden_exactly():
    Q = Quiver([1, 2], {"a": (1, 2)})
    assert tikz_quiver(Q, relations=[]) == GOLDEN.read_text()


def test_tikz_uses_the_same_layout_as_draw():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    L = layout(Q, relations=["a*b - c*d"])
    src = tikz_quiver(Q, relations=["a*b - c*d"])
    # every vertex node uses its exact layout coordinate (half-integers via pgf {p/q})
    assert "(v1) at (0, 0)" in src
    assert "(v2) at (1, {1/2})" in src
    assert "(v3) at (1, {-1/2})" in src
    assert "(v4) at (2, 0)" in src


def test_tikz_emits_relations_node_when_present():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    src = tikz_quiver(Q, relations=["a*b - c*d"])
    assert "relations: $a*b - c*d$" in src


def test_algebra_tikz_method():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(field=CC)
    assert A.tikz() == (pathlib.Path(__file__).parent / "golden" / "ka2.tikz").read_text()
