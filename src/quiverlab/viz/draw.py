"""A.draw(): render the exact layered layout with matplotlib (spec §3.7, §5 c.10).

FLOAT POLICY: no float/complex literal and no float() call in this file. Vertex
radius, curvature and offsets are int/Fraction; loop angles are integer degrees;
parallel bends use ConnectionStyle.Arc3(rad=<Fraction>) (the OBJECT form -- never
the "arc3,rad=0.3" string form, which would require formatting a float). Matplotlib
coerces the Fractions to float internally, outside src/quiverlab/.

Axis LIMITS are the one place matplotlib rejects Fraction bounds: matplotlib >= 3.11
runs np.isfinite on the limits, which raises TypeError on a Fraction (even an
integer-valued one). So the limits are snapped to ints with math.floor/math.ceil --
still float()-free (floor/ceil return int). Artist geometry (circles, arcs, arrows,
text) takes Fractions fine; only set_xlim/set_ylim need ints.

BACKEND: a Figure with its own FigureCanvasAgg is built directly -- no pyplot, no
matplotlib.use() -- so drawing never mutates the user's global backend. The returned
Figure renders inline (its Agg canvas) and fig.savefig handles PNG (Agg) and SVG
(matplotlib swaps to an SVG canvas by file extension)."""
import math
from fractions import Fraction

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.patches import Arc, Circle, ConnectionStyle, FancyArrowPatch

from quiverlab.viz.layout import layout

_VR = Fraction(1, 6)     # vertex circle radius
_LOOP_W = Fraction(1, 2)  # loop arc width/height


def draw_quiver(quiver, relations=(), file=None):
    L = layout(quiver, relations=relations)
    fig = Figure(figsize=(7, 5))
    FigureCanvasAgg(fig)                 # attach an Agg canvas WITHOUT touching pyplot
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect("equal")
    ax.axis("off")

    for v, (x, y) in L.positions.items():
        ax.add_patch(Circle((x, y), _VR, fill=False, linewidth=1, zorder=2))
        ax.text(x, y, str(v), ha="center", va="center", zorder=3)

    for e in L.edges:
        sx, sy = L.positions[e.src]
        tx, ty = L.positions[e.tgt]
        style = ConnectionStyle.Arc3(rad=e.bend)  # bend is a Fraction (0 for straight)
        arrow = FancyArrowPatch(
            (sx, sy), (tx, ty), connectionstyle=style, arrowstyle="-|>",
            mutation_scale=15, shrinkA=12, shrinkB=12, linewidth=1, zorder=1)
        ax.add_patch(arrow)
        mx, my = (sx + tx) * Fraction(1, 2), (sy + ty) * Fraction(1, 2) + e.bend
        ax.text(mx, my, e.name, ha="center", va="center", zorder=3)

    for lp in L.loops:
        cx, cy = L.positions[lp.at]
        # place the self-arc centered a little off the vertex along the base angle
        ax.add_patch(Arc((cx, cy + _VR), _LOOP_W, _LOOP_W,
                         angle=lp.angle_deg, theta1=200, theta2=520, linewidth=1, zorder=1))
        ax.text(cx, cy + _VR + _LOOP_W, lp.name, ha="center", va="bottom", zorder=3)

    xs = [x for (x, _y) in L.positions.values()]
    ys = [y for (_x, y) in L.positions.values()]
    pad = 1
    # Snap axis limits to ints: matplotlib>=3.11 rejects Fraction bounds (np.isfinite).
    x_lo, x_hi = math.floor(min(xs)) - pad, math.ceil(max(xs)) + pad
    y_lo, y_hi = math.floor(min(ys)) - pad, math.ceil(max(ys)) + pad + 1
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)

    if L.relations:
        ax.text(x_lo, y_lo, "relations:  " + ";  ".join(L.relations),
                ha="left", va="bottom", zorder=3)

    if file is not None:
        fig.savefig(str(file), bbox_inches="tight")
    return fig
