"""A.tikz(): the SAME layered layout as draw(), emitted as TikZ (spec §3.7).
Coordinates are exact: an integer prints as itself; a Fraction p/q prints as the
pgfmath expression {p/q}, which pgf evaluates -- so the emitted SOURCE contains
no float literal (this file is float-free like all of viz)."""
from fractions import Fraction

from quiverlab.viz.layout import layout, poset_layout


def _coord(z):
    z = Fraction(z)
    if z.denominator == 1:
        return str(z.numerator)
    return "{%d/%d}" % (z.numerator, z.denominator)


def tikz_quiver(quiver, relations=()):
    L = layout(quiver, relations=relations)
    lines = [r"\begin{tikzpicture}[>=stealth]"]
    for v in quiver.vertices:
        x, y = L.positions[v]
        lines.append(r"  \node[draw, circle] (v%s) at (%s, %s) {$%s$};"
                     % (v, _coord(x), _coord(y), v))
    for e in L.edges:
        if e.kind == "straight":
            lines.append(r"  \draw[->] (v%s) -- (v%s) node[midway, above] {$%s$};"
                         % (e.src, e.tgt, e.name))
        else:  # parallel: bend proportionally to the Fraction offset (integer degrees)
            deg = int(e.bend * 60)
            side = "left" if deg >= 0 else "right"
            lines.append(r"  \draw[->] (v%s) to[bend %s=%d] node[midway, above] {$%s$} (v%s);"
                         % (e.src, side, abs(deg), e.name, e.tgt))
    for lp in L.loops:
        lines.append(r"  \draw[->] (v%s) to[loop, in=%d, out=%d] node {$%s$} (v%s);"
                     % (lp.at, lp.angle_deg - 20, lp.angle_deg + 20, lp.name, lp.at))
    if L.relations:
        lines.append(r"  \node[align=left, below] at (current bounding box.south) "
                     r"{relations: %s};" % ";  ".join("$%s$" % r for r in L.relations))
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines) + "\n"


def _summand_math(vertex):
    """LaTeX label of a poset class from its ``summands`` [(name, mult), ...]:
    ``P_1``, ``S_1 \\oplus S_2``, ``S_1^{2} \\oplus S_2`` (name None -> a bullet)."""
    parts = []
    for name, mult in vertex.get("summands", ()):
        base = name if name else r"\bullet"
        parts.append(base if mult == 1 else "%s^{%d}" % (base, mult))
    return r" \oplus ".join(parts) if parts else "0"


def tikz_hasse(poset, label=_summand_math):
    """A standalone ``tikzpicture`` of a Hasse diagram (Plan 49 / C8): each poset
    class is a rounded node at its exact ranked (x, y); each cover is an undrawn-
    arrow line ``lo -- hi`` (lower class at the bottom). Coordinates exact (int /
    ``{p/q}`` pgfmath), float-free like the rest of ``viz``. Reusable by P45."""
    nodes = [v["index"] for v in poset.vertices]
    pos = poset_layout(nodes, poset.covers)
    lines = [r"\begin{tikzpicture}[>=stealth]"]
    for v in poset.vertices:
        i = v["index"]
        x, y = pos[i]
        lines.append(r"  \node[draw, rounded corners] (n%s) at (%s, %s) {$%s$};"
                     % (i, _coord(x), _coord(y), label(v)))
    for lo, hi in poset.covers:
        lines.append(r"  \draw (n%s) -- (n%s);" % (lo, hi))
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines) + "\n"


def tikz_order_complex(poset, face_vector=None):
    """The Hasse diagram of a ``families.poset.Poset`` plus the face vector of its ORDER
    COMPLEX (Plan 75 / R9) -- the report's picture of the object whose cohomology IS
    ``HH^*(kP)``.

    Distinct from :func:`tikz_hasse`, which draws the Plan-49 degeneration poset (nodes
    carrying ``summands``); here the nodes are the poset's own elements, labelled by
    ``str``. ``face_vector`` is stated verbatim when given -- it is NOT recomputed here,
    so this picture can never claim a count the engine did not produce."""
    nodes = list(poset.elements)
    index = {v: i for i, v in enumerate(nodes)}
    pos = poset_layout(nodes, list(poset.covers))
    lines = [r"\begin{tikzpicture}[>=stealth]"]
    for v in nodes:
        x, y = pos[v]
        lines.append(r"  \node[draw, circle] (p%d) at (%s, %s) {$%s$};"
                     % (index[v], _coord(x), _coord(y), v))
    for lo, hi in poset.covers:
        lines.append(r"  \draw (p%d) -- (p%d);" % (index[lo], index[hi]))
    if face_vector is not None:
        lines.append(r"  \node[align=left, below] at (current bounding box.south) "
                     r"{order complex $\Delta(P)$, face vector $(%s)$};"
                     % ", ".join(str(int(f)) for f in face_vector))
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines) + "\n"


def tikz_fan(fan):
    """The wall-and-chamber fan (Plan 45) as TikZ: for n=2 the g-vector rays drawn from
    the origin (exact coordinates, integer or {p/q}); for n=3 the L1/octahedron net
    positions the server pre-projected. Float-free (pgf evaluates {p/q}). Returns an empty
    picture with an honest note when the fan is budget-capped or n not in {2,3}."""
    n = fan.get("n")
    lines = [r"\begin{tikzpicture}[>=stealth, scale=2]"]
    if not fan.get("complete") or n not in (2, 3) or not fan.get("chambers"):
        lines.append(r"  \node {fan not drawn (n not in \{2,3\}, or budget-capped)};")
        lines.append(r"\end{tikzpicture}")
        return "\n".join(lines) + "\n"
    if n == 2:
        seen = set()
        for ch in fan["chambers"]:
            for ray in ch["rays"]:
                x, y = Fraction(ray[0]), Fraction(ray[1])
                key = (x, y)
                if key in seen or (x == 0 and y == 0):
                    continue
                seen.add(key)
                lines.append(r"  \draw[->, thick] (0,0) -- (%s, %s);"
                             % (_coord(x), _coord(y)))
        for w in (fan.get("walls") or []):
            bd = w.get("brick_dimvec") or {}
            lab = ",".join(str(bd[k]) for k in sorted(bd, key=str))
            nrm = w.get("normal")
            if nrm and len(nrm) == 2:
                x, y = Fraction(nrm[0]), Fraction(nrm[1])
                lines.append(r"  \node[font=\tiny, blue] at (%s, %s) {$(%s)$};"
                             % (_coord(x), _coord(y), lab))
    else:  # n == 3: draw the pre-projected octahedron-net rays
        seen = set()
        for ch in fan["chambers"]:
            for pt in (ch.get("net2d") or []):
                if pt is None:
                    continue
                x, y = Fraction(pt[0]), Fraction(pt[1])
                key = (x, y)
                if key in seen:
                    continue
                seen.add(key)
                lines.append(r"  \draw[->, thick] (0,0) -- (%s, %s);"
                             % (_coord(x), _coord(y)))
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines) + "\n"


def tikz_wall_chamber(b):
    """The Plan-63 wall-and-chamber structure as TikZ (the report's static twin of the GUI
    SVG): the chamber g-vector rays drawn faint, plus each brick-wall D(B) overlaid as a
    labelled dark ray/line (n=2 the exact extreme rays; n=3 the pre-projected octahedron-net
    facet rays). Float-free (pgf evaluates {p/q}). An honest note when the region is
    budget-capped only for the chambers, but walls are still drawn if present; a
    render='table' (n not in {2,3}) block returns the honest no-drawing note."""
    n = b.get("n")
    lines = [r"\begin{tikzpicture}[>=stealth, scale=2]"]
    if b.get("render") not in ("fan2d", "fan3d") or not b.get("chambers"):
        lines.append(r"  \node {fan not drawn (n not in \{2,3\}, or nothing found)};")
        lines.append(r"\end{tikzpicture}")
        return "\n".join(lines) + "\n"
    use3 = (n == 3)
    seen = set()
    for ch in b["chambers"]:
        rays = (ch.get("net2d") or []) if use3 else ch.get("rays", [])
        for ray in rays:
            if ray is None:
                continue
            x, y = Fraction(ray[0]), Fraction(ray[1])
            if (x, y) in seen or (x == 0 and y == 0):
                continue
            seen.add((x, y))
            lines.append(r"  \draw[->, gray!50] (0,0) -- (%s, %s);"
                         % (_coord(x), _coord(y)))
    for w in (b.get("walls") or []):
        wrays = (w.get("net2d") or []) if use3 else (w.get("rays") or [])
        bd = w.get("brick_dimvec") or {}
        lab = w.get("brick_name") or (",".join(str(bd[k]) for k in sorted(bd, key=str)))
        for ray in wrays:
            if ray is None:
                continue
            x, y = Fraction(ray[0]), Fraction(ray[1])
            if x == 0 and y == 0:
                continue
            lines.append(r"  \draw[->, very thick] (0,0) -- (%s, %s);"
                         % (_coord(x), _coord(y)))
            lines.append(r"  \node[font=\tiny] at (%s, %s) {$%s$};"
                         % (_coord(x), _coord(y), _esc_tex(lab)))
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines) + "\n"


def _esc_tex(s):
    """Minimal TeX escaping for a wall label (dim-vector string or brick name)."""
    return str(s).replace("\\", r"\textbackslash{}").replace("_", r"\_").replace(
        "{", r"\{").replace("}", r"\}").replace("$", r"\$")
def tikz_koszul(b):
    """The Plan-77 generalized-Koszulity staircase: the internal generation degree
    l_i(n) plotted against the homological degree n, with the diagonal l = n drawn
    faint.  A Koszul strand hugs the diagonal; Berger's 2-N alternation climbs in
    steps 1, N-1, 1, N-1, ...; an almost-Koszul algebra runs ON the diagonal and then
    leaves it exactly once, at the break.

    Every number is taken VERBATIM from the block -- nothing is recomputed here, so the
    picture cannot claim a degree the engine did not produce.  The certified window is
    labelled, never hidden.  Returns an honest note when the block carries no table.
    """
    gd = (b or {}).get("generation_degrees") or {}
    lines = [r"\begin{tikzpicture}[>=stealth, scale=0.6]"]
    if not gd:
        lines.append(r"  \node {no generation-degree table recorded};")
        lines.append(r"\end{tikzpicture}")
        return "\n".join(lines) + "\n"
    keys = sorted(gd, key=lambda k: int(k))
    width = max((len(gd[k]) for k in keys), default=0)
    top = max((max(gd[k]) for k in keys if gd[k]), default=0)
    lines.append(r"  \draw[->] (0,0) -- (%d,0) node[right] {$n$};" % (width + 1))
    lines.append(r"  \draw[->] (0,0) -- (0,%d) node[above] {$\ell$};" % (top + 1))
    lines.append(r"  \draw[gray!40] (0,0) -- (%d,%d) node[right, gray] "
                 r"{$\ell = n$};" % (min(width, top), min(width, top)))
    marks = ["*", "square*", "triangle*", "diamond*"]
    for idx, k in enumerate(keys):
        row = gd[k]
        if not row:
            continue
        pts = " -- ".join("(%d,%d)" % (n, d) for n, d in enumerate(row))
        lines.append(r"  \draw[thick, mark=%s] %s;" % (marks[idx % len(marks)], pts))
        lines.append(r"  \node[right, font=\small] at (%d,%d) {$S_{%s}$};"
                     % (len(row) - 1, row[-1], k))
    cdeg = (b or {}).get("certified_through_degree")
    if cdeg is not None:
        lines.append(r"  \node[below, font=\small] at (%d,-0.6) "
                     r"{certified through degree %s};" % (max(width // 2, 1), cdeg))
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines) + "\n"
