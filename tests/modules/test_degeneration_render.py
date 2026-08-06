"""The Hasse twin renderers (Plan 49 / C8). Self-cert: the layout ranks the poset
(minimum at rank 0, cover endpoints one rank apart); tikz_hasse / hasse_svg emit
non-empty markup naming every class and drawing every cover. Float-free."""
import pytest

from fractions import Fraction

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.degeneration import degeneration_order
from quiverlab.viz.hasse_html import hasse_svg
from quiverlab.viz.layout import poset_layout
from quiverlab.viz.tikz import tikz_hasse

pytestmark = pytest.mark.oracle_selfcert


def _diamond():
    return degeneration_order(linear_path_algebra(3, field=QQ), {1: 1, 2: 1, 3: 1})


def test_layout_ranks_the_poset():
    P = _diamond()
    nodes = [v["index"] for v in P.vertices]
    pos = poset_layout(nodes, P.covers)
    assert all(isinstance(pos[i][0], (int, Fraction)) for i in nodes)   # exact x
    for lo, hi in P.covers:
        assert pos[hi][1] == pos[lo][1] + 1        # a cover spans exactly one rank
    ys = {pos[i][1] for i in nodes}
    assert min(ys) == 0 and max(ys) == 2           # diamond has 3 ranks


def test_tikz_and_svg_are_nonempty_and_name_classes():
    P = _diamond()
    tk = tikz_hasse(P)
    svg = hasse_svg(P)
    assert tk.strip().startswith(r"\begin{tikzpicture}") and r"\end{tikzpicture}" in tk
    assert svg.strip().startswith("<svg") and "</svg>" in svg
    # every cover is drawn (4 edges in the diamond)
    assert tk.count("--") >= len(P.covers)
    assert svg.count("<line") >= len(P.covers)
