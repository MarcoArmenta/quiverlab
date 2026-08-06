"""A.tikz(): the SAME layered layout as draw(), emitted as TikZ (spec §3.7).
Coordinates are exact: an integer prints as itself; a Fraction p/q prints as the
pgfmath expression {p/q}, which pgf evaluates -- so the emitted SOURCE contains
no float literal (this file is float-free like all of viz)."""
from fractions import Fraction

from quiverlab.viz.layout import layout


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
