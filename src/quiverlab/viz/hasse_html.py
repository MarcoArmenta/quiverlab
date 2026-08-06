"""Inline SVG Hasse diagram (Plan 49 / C8) -- the HTML/GUI twin of
``viz.tikz.tikz_hasse``. The exact rank/Fraction layout comes from
``viz.layout.poset_layout``; the ONLY integer cast is the Fraction -> pixel scale
(``numerator // denominator`` on the scaled rational -- a provably integer floor,
never a float literal, so this file passes the no-floats gate like the rest of
``viz``). Reusable by P45's exchange lattice."""
from quiverlab.viz.layout import poset_layout


def _esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _summand_text(vertex):
    """Plain-text label of a poset class from its ``summands`` [(name, mult), ...]:
    ``P_1``, ``S_1 (+) S_2``, ``S_1^2 (+) S_2`` (name None -> ``M``)."""
    parts = []
    for name, mult in vertex.get("summands", ()):
        base = name if name else "M"
        parts.append(base if mult == 1 else "%s^%d" % (base, mult))
    return " ⊕ ".join(parts) if parts else "0"


def _px(frac, scale):
    """Fraction * integer scale -> integer pixel (exact rational floor; no float)."""
    f = frac * scale
    return f.numerator // f.denominator


def hasse_svg(poset, label=_summand_text, scale=110, pad=46):
    """An inline ``<svg>`` Hasse diagram of a poset (lower classes at the bottom):
    ``<line>`` per cover, ``<circle>`` + ``<text>`` per class (the generic/open-
    orbit class tinted). Integer pixels via ``_px``; float-free."""
    nodes = [v["index"] for v in poset.vertices]
    pos = poset_layout(nodes, poset.covers)
    xs = [pos[i][0] for i in nodes]
    ys = [pos[i][1] for i in nodes]
    minx = min(xs) if xs else 0
    maxrank = max(ys) if ys else 0
    span = (max(xs) - minx) if xs else 0
    width = pad * 2 + _px(span, scale)
    height = pad * 2 + maxrank * scale
    coord = {i: (pad + _px(pos[i][0] - minx, scale),
                 pad + (maxrank - pos[i][1]) * scale) for i in nodes}
    out = ["<svg viewBox='0 0 %d %d' role='img' xmlns='http://www.w3.org/2000/svg' "
           "style='max-width:%dpx;width:100%%;display:block'>" % (width, height, width)]
    for lo, hi in poset.covers:                 # covers behind the nodes
        x1, y1 = coord[lo]
        x2, y2 = coord[hi]
        out.append("<line x1='%d' y1='%d' x2='%d' y2='%d' stroke='#444' "
                   "stroke-width='2'/>" % (x1, y1, x2, y2))
    for v in poset.vertices:
        i = v["index"]
        x, y = coord[i]
        fill = "#e8eefc" if v.get("is_generic") else "#fff"
        out.append("<circle cx='%d' cy='%d' r='16' fill='%s' stroke='#3f51b5' "
                   "stroke-width='2'/>" % (x, y, fill))
        out.append("<text x='%d' y='%d' font-size='13' text-anchor='middle' "
                   "fill='#1c1c1c'>%s</text>" % (x, y + 30, _esc(label(v))))
    out.append("</svg>")
    return "\n".join(out) + "\n"
