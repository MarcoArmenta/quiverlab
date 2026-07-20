"""A.draw(): matplotlib rendering from the exact layout. Structure-only asserts.
No matplotlib.use() here -- draw_quiver builds its own Agg canvas, so no global
backend is needed or mutated."""
from fractions import Fraction

from quiverlab import Quiver, CC
from quiverlab.viz.draw import draw_quiver
from quiverlab.viz.layout import layout


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})


def test_draw_returns_figure_with_one_axes():
    fig = draw_quiver(_square(), relations=["a*b - c*d"])
    assert fig.__class__.__name__ == "Figure"
    assert len(fig.axes) == 1


def test_draw_places_one_label_per_vertex_at_layout_coords():
    Q = _square()
    fig = draw_quiver(Q, relations=["a*b - c*d"])
    ax = fig.axes[0]
    L = layout(Q, relations=["a*b - c*d"])
    texts = {t.get_text(): t.get_position() for t in ax.texts}
    for v, (x, y) in L.positions.items():
        assert str(v) in texts
        px, py = texts[str(v)]
        # matplotlib stores float positions; compare to the exact layout coerced.
        assert (px, py) == (float(x), float(y))


def test_draw_renders_relation_list_text():
    fig = draw_quiver(_square(), relations=["a*b - c*d"])
    ax = fig.axes[0]
    assert any("a*b - c*d" in t.get_text() for t in ax.texts)


def test_draw_one_patch_per_arrow_plus_loops():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    fig = draw_quiver(Q)
    ax = fig.axes[0]
    # two loops -> two Arc patches; no straight/parallel edges here.
    arcs = [p for p in ax.patches if p.__class__.__name__ == "Arc"]
    assert len(arcs) == 2


def test_draw_file_export_png_and_svg(tmp_path):
    Q = _square()
    png = tmp_path / "A.png"
    svg = tmp_path / "A.svg"
    draw_quiver(Q, relations=["a*b - c*d"], file=str(png))
    draw_quiver(Q, relations=["a*b - c*d"], file=str(svg))
    assert png.exists() and png.stat().st_size > 0
    assert svg.exists() and svg.read_text().lstrip().startswith("<?xml")


def test_algebra_draw_method():
    A = _square().algebra(relations=["a*b - c*d"], field=CC)
    fig = A.draw()
    assert fig.__class__.__name__ == "Figure"
