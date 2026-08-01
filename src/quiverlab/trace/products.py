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
from quiverlab.trace.events import ProductBasis, ProductStep, ResultDims, StepNote

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


def notation_legend(kind, degrees_note, basis):
    """One-sentence legend DEFINING the symbols the product-table equations use
    (Marco: "I need the definitions of alphas, betas, zetas..."). Shown immediately
    before each product family's tables in every render surface.

    ``basis`` is the concrete recorded class-basis string (e.g. ``"bar/GF(7)"`` /
    ``"cs/..."``) taken from the result block -- the structure constants are
    basis-dependent, so the reader is told which basis; ``None`` (connes_b, which has
    no single basis string) omits it. ``degrees_note`` is an optional trailing clause
    (the degree indices p/q/n are otherwise per-table); appended when non-empty.

    Plain prose with Unicode math glyphs (alpha/beta/gamma/z/w), rendered as escaped
    text by every surface -- it never needs typesetting. This is the SINGLE source of
    the legend text, consumed by ``render_html`` (via ``products_chapter``) and
    ``results_html._product_tables_html``; the GUI hardcodes the same wording."""
    on_basis = ("relative to the recorded basis %s" % basis if basis
                else "relative to the recorded class basis")
    # Plan-35 UNIT 2 (Marco): the legend no longer only NAMES the symbols -- it points
    # at the explicit per-degree listings, where every α/β/γ/z/w is printed as its
    # (co)cycle term-sum + coordinate vector with the annihilating differential. β_j /
    # γ_k / w_k are the SAME degree classes as the α_i / z_j listed there, viewed as
    # the operand / output of this product.
    _EXPLICIT = (" Each class is listed explicitly by degree above -- its term-sum, "
                 "its coordinate vector, and the differential that annihilates it.")
    if kind == "cup":
        s = ("α₁,…,α_{d_p} are the recorded basis classes of "
             "HH^p and β₁,…,β_{d_q} those of HH^q; "
             "γ₁,…,γ_{d_{p+q}} the basis of HH^{p+q}. Every table "
             "line states α_i ∪ β_j = Σ_k c·γ_k, %s; "
             "the constants c are basis-dependent.%s" % (on_basis, _EXPLICIT))
    elif kind == "bracket":
        s = ("α₁,…,α_{d_p} are the recorded basis classes of "
             "HH^p and β₁,…,β_{d_q} those of HH^q; "
             "γ₁,…,γ_{d_{p+q-1}} the basis of HH^{p+q-1}. Every "
             "table line states [α_i, β_j] = Σ_k c·γ_k in "
             "degree p+q−1, %s; the constants c are basis-dependent.%s"
             % (on_basis, _EXPLICIT))
    elif kind == "cap":
        s = ("z₁,…,z_{d_n} are the recorded basis classes of HH_n "
             "(homology) and w₁,…,w_{d_{n-p}} those of HH_{n-p}; "
             "α₁,…,α_{d_p} the basis of HH^p. Every table line "
             "states α_i ∩ z_j = Σ_k c·w_k, %s; the constants c "
             "are basis-dependent.%s" % (on_basis, _EXPLICIT))
    elif kind == "connes_b":
        s = ("each induced Connes differential B_n: HH_n → HH_{n+1} is written "
             "on the recorded homology bases -- rows index HH_{n+1}, columns index "
             "HH_n; the entries are basis-dependent. The homology cycle classes z^n_j "
             "are listed explicitly by degree above, each with its coordinate vector "
             "and the boundary b_n that annihilates it.")
    else:
        raise QuiverlabError(
            "unknown product kind %r for the notation legend" % (kind,))
    return "%s (%s)" % (s, degrees_note) if degrees_note else s


def _basis_events(kind, obj):
    """The Plan-35 UNIT-2 ``ProductBasis`` event carrying the explicit
    representatives (per (side, degree) classes + ordered enumeration + annihilating
    differential), extracted from the product object's ``blocks()``. Returns ``[]``
    for a legacy object that recorded no reps -- the chapter then shows tables only."""
    b = obj.blocks()
    if b.get("basis_classes") is None:
        return []
    return [ProductBasis(kind=kind, basis_classes=b.get("basis_classes"),
                         chain_basis=b.get("chain_basis"),
                         differentials=b.get("differentials"))]


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
    # The notation legend defining the table symbols, named the recorded basis
    # (Marco, Plan-35 follow-up) -- one shared builder for every render surface.
    return ([StepNote(title, prose, heading=True),
             StepNote("Notation.", notation_legend(kind, "", obj.basis)),
             ResultDims(kind=result_kind, dims=list(result_dims))]
            + _basis_events(kind, obj) + steps)


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


def equation_lines(kind, degrees, out_degree, dims, constants):
    """A TeX equation per nonzero product ``left_i * right_j``: the RHS lists the
    output generators with their exact structure-constant coefficients (the unit
    coefficient ``1`` is left implicit; every non-unit coefficient is shown, so
    every nonzero constant appears verbatim).

    Operates on the RAW table data (``degrees``/``out_degree``/``dims``/
    ``constants``) rather than a ``ProductTable``, so the SAME builder serves both
    the worked-steps chapter (from the frozen object, via ``_equation_lines``) and
    the report's Computed-results block (from the serialized ``blocks()`` dict, via
    ``trace.results_html``) -- one implementation, no divergent copy. Zero terms are
    skipped; a fully vanishing table yields ``[]`` (the caller emits the vanish
    note)."""
    dl, dr, dout = dims
    left_deg, right_deg = degrees
    lines = []
    for i in range(dl):
        for j in range(dr):
            terms = [_term(kind, constants[k][i][j], out_degree, k)
                     for k in range(dout) if str(constants[k][i][j]) != "0"]
            if terms:
                lines.append("%s = %s" % (
                    _lhs(kind, left_deg, i, right_deg, j), " + ".join(terms)))
    return lines


def _equation_lines(kind, t):
    """The equation lines of a ``ProductTable`` (the worked-steps chapter path)."""
    return equation_lines(kind, t.degrees, t.out_degree, t.dims, t.constants)


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
              StepNote("Notation.", notation_legend("connes_b", "", None)),
              ResultDims(kind="HH_", dims=list(obj.hh_dims))]
    events += _basis_events("connes_b", obj)
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
