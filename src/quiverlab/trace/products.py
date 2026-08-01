"""Plan 35 -- the worked-steps chapter for the HH product surface.

``products_chapter(A, kind, obj)`` turns a Task-1 product result object (an
``HHProducts`` for cup/cap/bracket, a ``ConnesB`` for connes_b) into the typed
trace-event stream the renderers turn into a homework-grade chapter:

  * a prose ``StepNote`` naming the product and citing its source;
  * a ``ResultDims`` carrying the authoritative Hochschild (co)homology dimensions
    (rendered as the chapter's HH result table) -- and, in building it, the DRIFT
    GATE: every dimension the product object recorded is compared, degree by
    degree, against a FRESH ``hochschild_(co)homology`` on ``A``; a mismatch raises
    ``QuiverlabError`` rather than shipping a chapter that misstates the
    computation (mirroring the module drift gates in ``trace.modules``);
  * one ``ProductStep`` per bidegree -- the nonzero structure-constant equations
    for cup/cap/bracket, the induced Connes differential matrix for connes_b.

The gate is genuine, not decorative: it re-derives the HH dimensions INDEPENDENTLY
of the (possibly corrupted) product object, so a tampered table dimension is caught.

Float-free: every number is an int or an exact domain-element string. The fixed
per-kind definitional formula lives in ``render_html`` (the sole owner of the
worked-steps math source); this module owns only the DATA-derived equations.
"""
from quiverlab.errors import QuiverlabError
from quiverlab.trace.events import ProductStep, ResultDims, StepNote

# (heading text, prose) per kind. Plain ASCII -- a StepNote renders as escaped text,
# not typeset math (the typeset definition formula is the renderer's job). "cup
# product" / "cap product" / "Gerstenhaber bracket" / "Connes" appear verbatim so a
# reader (and the section title) name the object.
_PROSE = {
    "cup": ("The cup product",
            "The cup product HH^p (x) HH^q -> HH^{p+q} makes the Hochschild "
            "cohomology HH^*(A) a graded-commutative ring (Gerstenhaber, The "
            "cohomology structure of an associative ring, Ann. of Math. 78 (1963)). "
            "On the recorded class basis it is the structure-constant table below."),
    "cap": ("The cap product",
            "The cap product HH^p (x) HH_n -> HH_{n-p} is the action of the "
            "cohomology ring HH^*(A) on Hochschild homology HH_*(A) (Gerstenhaber "
            "1963). Its structure constants in the recorded bases are below."),
    "bracket": ("The Gerstenhaber bracket",
            "The Gerstenhaber bracket [-,-]: HH^p (x) HH^q -> HH^{p+q-1} is the "
            "graded Lie bracket that, with the cup product, makes HH^*(A) a "
            "Gerstenhaber algebra (Gerstenhaber 1963). Bracket structure constants "
            "in the recorded bases (over the served window) are below."),
    "connes_b": ("The Connes differential B",
            "The Connes (Rinehart) operator B: HH_n -> HH_{n+1} is the boundary of "
            "the (b, B) mixed complex computing cyclic homology; B^2 = 0 at the "
            "induced level. Each induced matrix in the class bases is below."),
}

# Display symbols for the operand / output generators of a table equation. Left is
# always a cohomology class; the right and output are cohomology for cup/bracket and
# homology for cap (different letters signal the variance).
_LEFT = r"\alpha"
_RIGHT = {"cup": r"\beta", "bracket": r"\beta", "cap": r"z"}
_OUT = {"cup": r"\gamma", "bracket": r"\gamma", "cap": r"w"}


def products_chapter(A, kind, obj):
    """The worked-steps event stream for the product ``kind`` computed as ``obj``.

    Raises ``QuiverlabError`` if ``obj``'s recorded dimensions drift from a fresh
    Hochschild (co)homology on ``A`` (the drift gate)."""
    if kind == "connes_b":
        return _connes_chapter(A, obj)
    if kind in ("cup", "cap", "bracket"):
        return _table_chapter(A, kind, obj)
    raise QuiverlabError(
        "unknown product kind %r for the worked-steps chapter" % (kind,))


# --------------------------------------------------------------------------- #
# cup / cap / bracket: structure-constant tables as equation lines.
# --------------------------------------------------------------------------- #

def _table_chapter(A, kind, obj):
    top = obj.top
    # verbose=False: building the chapter must not itself spew a stray trace report
    # (hochschild_* defer to quiverlab.verbose otherwise), exactly as spec._dispatch.
    coh = list(A.hochschild_cohomology(top, verbose=False).dims)
    hom = (list(A.hochschild_homology(top, verbose=False).dims)
           if kind == "cap" else None)
    steps = []
    for key in sorted(obj.tables):
        t = obj.tables[key]
        want = _expected_dims(kind, t.degrees, coh, hom)
        if tuple(t.dims) != tuple(want):
            raise QuiverlabError(
                "products chapter drift: the %s table %s records dims %s but the "
                "Hochschild dimensions are %s -- refusing to narrate a table that "
                "misstates the computation"
                % (kind, tuple(t.degrees), tuple(t.dims), tuple(want)))
        steps.append(_table_step(kind, t))
    title, prose = _PROSE[kind]
    result_kind = "HH_" if kind == "cap" else "HH^"
    result_dims = hom if kind == "cap" else coh
    return [StepNote(title, prose, heading=True),
            ResultDims(kind=result_kind, dims=list(result_dims))] + steps


def _expected_dims(kind, degrees, coh, hom):
    """The authoritative ``(dim_left, dim_right, dim_out)`` a table of this bidegree
    MUST record, read off a fresh HH computation -- the drift-gate reference."""
    p, second = degrees
    if kind == "cup":
        return (coh[p], coh[second], coh[p + second])
    if kind == "bracket":
        return (coh[p], coh[second], coh[p + second - 1])
    return (coh[p], hom[second], hom[second - p])          # cap: (HH^p, HH_n, HH_{n-p})


def _table_step(kind, t):
    """One bidegree's ProductStep: its map heading and the nonzero-product equation
    lines (an all-vanishing bidegree carries a note instead)."""
    lines = _equation_lines(kind, t)
    note = "" if lines else "every product in this bidegree vanishes."
    return ProductStep(kind=kind, degrees=tuple(t.degrees),
                       heading=_map_heading(kind, t.degrees, t.out_degree),
                       lines=lines, matrix=None, note=note)


def _map_heading(kind, degrees, out_degree):
    p, second = degrees
    if kind == "cap":
        return r"HH^{%d} \otimes HH_{%d} \to HH_{%d}" % (p, second, out_degree)
    return r"HH^{%d} \otimes HH^{%d} \to HH^{%d}" % (p, second, out_degree)


def _equation_lines(kind, t):
    """A TeX equation per nonzero product ``left_i * right_j``: the RHS lists the
    output generators with their exact structure-constant coefficients (the unit
    coefficient ``1`` is left implicit; every non-unit coefficient is shown, so
    every nonzero constant appears verbatim)."""
    dl, dr, dout = t.dims
    out_deg = t.out_degree
    left_deg, right_deg = t.degrees
    lines = []
    for i in range(dl):
        for j in range(dr):
            terms = []
            for k in range(dout):
                c = t.constants[k][i][j]
                if c != "0":
                    terms.append(_term(kind, c, out_deg, k))
            if terms:
                lines.append("%s = %s" % (
                    _lhs(kind, left_deg, i, right_deg, j), " + ".join(terms)))
    return lines


def _lhs(kind, left_deg, i, right_deg, j):
    L = r"%s^{%d}_{%d}" % (_LEFT, left_deg, i + 1)
    R = r"%s^{%d}_{%d}" % (_RIGHT[kind], right_deg, j + 1)
    if kind == "cup":
        return r"%s \cup %s" % (L, R)
    if kind == "cap":
        return r"%s \cap %s" % (L, R)
    return r"[%s, %s]" % (L, R)                        # bracket


def _term(kind, c, out_deg, k):
    g = r"%s^{%d}_{%d}" % (_OUT[kind], out_deg, k + 1)
    return g if c == "1" else r"%s\,%s" % (c, g)


# --------------------------------------------------------------------------- #
# connes_b: the induced Connes differentials B: HH_n -> HH_{n+1} as matrices.
# --------------------------------------------------------------------------- #

def _connes_chapter(A, obj):
    top = obj.top
    hom = list(A.hochschild_homology(top, verbose=False).dims)   # no stray trace
    if list(obj.hh_dims) != hom:
        raise QuiverlabError(
            "products chapter drift: the Connes object records HH_ dimensions %s "
            "but a fresh Hochschild homology gives %s -- refusing to narrate a "
            "chapter that misstates the computation" % (list(obj.hh_dims), hom))
    title, prose = _PROSE["connes_b"]
    events = [StepNote(title, prose, heading=True),
              ResultDims(kind="HH_", dims=list(obj.hh_dims))]
    for n in sorted(obj.matrices):
        mat = obj.matrices[n]
        rows = len(mat)
        cols = len(mat[0]) if (mat and mat[0]) else 0
        if rows != obj.hh_dims[n + 1] or (rows and cols != obj.hh_dims[n]):
            raise QuiverlabError(
                "products chapter drift: the induced B_%d matrix is %dx%d but the "
                "HH_ dimensions demand %dx%d"
                % (n, rows, cols, obj.hh_dims[n + 1], obj.hh_dims[n]))
        events.append(ProductStep(
            kind="connes_b", degrees=(n,),
            heading=r"B_{%d} : HH_{%d} \to HH_{%d}" % (n, n, n + 1),
            lines=(), matrix=mat, note="induced rank %d." % obj.ranks.get(n, 0)))
    return events
