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
import re

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

# Display symbols for the operand / output generators of a table equation. Marco
# 2026-07-31: notation is UNIFORM -- every cohomology class is written alpha^n_j (the
# j-th basis class of HH^n; superscript = degree, subscript = index) and every
# homology class z^n_j. So cup/bracket (all cohomology) use alpha throughout, and cap
# (cohomology acting on homology) uses alpha on the left, z on the right and output.
# beta / gamma / w are gone.
_LEFT = r"\alpha"
_RIGHT = {"cup": r"\alpha", "bracket": r"\alpha", "cap": r"z"}
_OUT = {"cup": r"\alpha", "bracket": r"\alpha", "cap": r"z"}

# The product operator shown in the Cayley table's CORNER cell (the mnemonic that
# every cell of the grid is left-class OP right-class): cup ∪, cap ∩, bracket [-,-].
_CORNER = {"cup": r"\cup", "cap": r"\cap", "bracket": r"[-,-]"}

# The GF(p) prime carried in a recorded basis string ("bar/GF(7)", "cs/GF(2)"); a
# non-GF basis (QQ, an extension-field repr) yields None -> no balancing, coefficients
# verbatim.
_GF_RE = re.compile(r"GF\((\d+)\)")


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
    # Marco 2026-07-31: UNIFORM notation, α^n_j = the j-th basis class of HH^n and
    # z^n_j = the j-th of HH_n (superscript = degree, subscript = index); β/γ/w are
    # gone. The legend names the symbols, states the product in that notation, and
    # points at the explicit per-degree listings (where each class is written as a
    # combination of the ordered basis elements). Coordinate vectors live in the JSON.
    _EXPLICIT = (" Each class is listed explicitly by degree above as a combination of "
                 "the ordered basis elements; the coordinate vectors are recorded in "
                 "the JSON.")
    if kind == "cup":
        s = ("α^n_j denotes the j-th basis class of HH^n (superscript = degree, "
             "subscript = index). Every table line states "
             "α^p_i ∪ α^q_j = Σ_k c·α^{p+q}_k, %s; "
             "the constants c are basis-dependent.%s" % (on_basis, _EXPLICIT))
    elif kind == "bracket":
        s = ("α^n_j denotes the j-th basis class of HH^n (superscript = degree, "
             "subscript = index). Every table line states "
             "[α^p_i, α^q_j] = Σ_k c·α^{p+q-1}_k in degree p+q−1, %s; "
             "the constants c are basis-dependent.%s" % (on_basis, _EXPLICIT))
    elif kind == "cap":
        s = ("α^p_j denotes the j-th basis class of HH^p (cohomology) and z^n_j the "
             "j-th of HH_n (homology). Every table line states "
             "α^p_i ∩ z^n_j = Σ_k c·z^{n-p}_k, %s; the constants c "
             "are basis-dependent.%s" % (on_basis, _EXPLICIT))
    elif kind == "connes_b":
        s = ("each induced Connes differential B_n: HH_n → HH_{n+1} is written "
             "on the recorded homology bases -- rows index HH_{n+1}, columns index "
             "HH_n; the entries are basis-dependent. The homology cycle classes z^n_j "
             "(z^n_j = the j-th basis class of HH_n) are listed explicitly by degree "
             "above as combinations of the ordered basis elements.")
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
    prime = prime_from_basis(obj.basis)
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
        steps.append(_table_step(kind, t, prime))
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


def _table_step(kind, t, prime=None):
    """One bidegree's ProductStep: its map heading, the raw structure-constant data the
    HTML surface turns into a Cayley grid, AND the nonzero-product equation lines
    (kept for the JSON prose and as the vanish sentinel -- ``lines`` empty <=> the
    whole table vanishes, which carries the one-line note instead of a grid)."""
    lines = _equation_lines(kind, t)
    note = "" if lines else "every product in this bidegree vanishes."
    return ProductStep(kind=kind, degrees=tuple(t.degrees),
                       heading=_map_heading(kind, t.degrees, t.out_degree),
                       lines=lines, matrix=None, note=note,
                       dims=list(t.dims), constants=t.constants,
                       out_degree=t.out_degree, prime=prime)


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
# Cayley multiplication tables (Marco 2026-08-01). Every product bidegree renders
# as a grid: rows = left classes, cols = right classes, each cell the product
# expressed DIRECTLY in the target basis (0 / a signed linear combination). Zeros
# are SHOWN inside a table that is nonzero somewhere -- they are information.
#
# This module owns the pure/data-derived part -- the cell TeX, the balanced-rep
# coefficient display, the structural-note derivations -- as float-free functions of
# the recorded structure constants. Each HTML surface (render_html, results_html)
# turns the returned struct into a grid via ``render_html.cayley_grid_html``; the GUI
# mirrors the same logic in JS. ``equation_lines`` stays exported (JSON prose / the
# vanish test), but the RENDERED product tables are Cayley grids.
# --------------------------------------------------------------------------- #

def prime_from_basis(basis):
    """The prime ``p`` of a recorded GF(p) basis string, else ``None`` (QQ / an
    extension-field repr -> coefficients shown verbatim, no balancing)."""
    if not basis:
        return None
    m = _GF_RE.search(str(basis))
    return int(m.group(1)) if m else None


def balanced_coeff(c, prime):
    """A GF(p) residue displayed as its BALANCED representative: ``c`` with ``c > p/2``
    (equivalently ``2c > p``) shows as ``c - p`` (so p-1 -> -1, p-2 -> -2). Display
    only -- the JSON keeps the raw residue string. A coefficient that is not a residue
    in ``[0, p)`` (a fraction, an out-of-range or already-signed value, or no prime at
    all) is returned verbatim."""
    c = str(c)
    if prime is None:
        return c
    try:
        v = int(c)
    except ValueError:
        return c
    if 0 <= v < prime and 2 * v > prime:
        return str(v - prime)
    return c


def balanced_rep_note(prime):
    """The one-sentence legend for the balanced-representative display, or ``""`` for a
    non-GF basis (nothing to balance)."""
    if prime is None:
        return ""
    return ("Coefficients are shown as balanced representatives mod %d (a residue "
            "c > %d/2 is written c-%d); the JSON record keeps the raw residues."
            % (prime, prime, prime))


def _signed_join(pieces):
    """Join ``(coeff_display, generator_tex)`` terms into one signed TeX sum, the unit
    magnitude left implicit (``\\alpha^{2}_{1}``, ``-\\alpha^{2}_{1}``,
    ``\\alpha^{2}_{1} + 2\\,\\alpha^{2}_{3}``). ``[]`` -> ``"0"``."""
    if not pieces:
        return "0"
    out = []
    for idx, (c, g) in enumerate(pieces):
        neg = c.startswith("-")
        mag = c[1:] if neg else c
        term = g if mag == "1" else r"%s\,%s" % (mag, g)
        if idx == 0:
            out.append("-" + term if neg else term)
        else:
            out.append((" - " if neg else " + ") + term)
    return "".join(out)


def cell_tex(kind, out_degree, coeffs, prime):
    """One Cayley cell: the product expressed in the TARGET basis as
    ``Σ_k balanced(c_k)·g_k`` (``g_k`` the k-th output generator), or ``"0"`` when
    every coefficient vanishes. ``coeffs`` is ``[constants[k][i][j] for k]``."""
    pieces = []
    for k, c in enumerate(coeffs):
        disp = balanced_coeff(c, prime)
        if str(disp) == "0":
            continue
        pieces.append((disp, r"%s^{%d}_{%d}" % (_OUT[kind], out_degree, k + 1)))
    return _signed_join(pieces)


def _is_int(s):
    try:
        int(str(s))
        return True
    except ValueError:
        return False


def _mirror_sign(kind, p, q):
    """The graded sign relating a product table to its transpose: cup satisfies
    ``α^p∪α^q = (-1)^{pq} α^q∪α^p``; the Gerstenhaber bracket
    ``[α^p,α^q] = -(-1)^{(p-1)(q-1)} [α^q,α^p]``. Returns +1 or -1."""
    if kind == "cup":
        return -1 if (p * q) % 2 else 1
    return -1 if ((p - 1) * (q - 1)) % 2 == 0 else 1        # bracket


def structural_notes(kind, degrees, dims, constants, prime):
    """Honest structural observations DERIVED from the constants (never asserted
    blindly):

      * ``"all squares are 0"`` when every diagonal product ``x·x`` vanishes;
      * ``"the table is graded-antisymmetric"`` / ``"...graded-commutative
        (symmetric)"`` when the table equals its sign-mirrored transpose, the sign read
        off the degrees by :func:`_mirror_sign`.

    Only meaningful for cup/bracket on a SQUARE bidegree (``p == q`` and
    ``dl == dr``) -- the two operands then live in the same space, so the diagonal and
    the transpose make sense; ``[]`` otherwise. The sign-mirror note needs residue
    arithmetic and is emitted only over a GF(p) basis."""
    if kind not in ("cup", "bracket"):
        return []
    p, q = degrees
    dl, dr, dout = dims
    if p != q or dl != dr or not dl:
        return []
    n = dl
    notes = []
    if all(str(constants[k][i][i]) == "0" for i in range(n) for k in range(dout)):
        notes.append("all squares are 0")
    if prime is not None:
        all_int = all(_is_int(constants[k][i][j])
                      for i in range(n) for j in range(n) for k in range(dout))
        if all_int:
            sign = _mirror_sign(kind, p, q)
            mirrored = all(
                (int(constants[k][i][j]) - sign * int(constants[k][j][i])) % prime == 0
                for i in range(n) for j in range(n) for k in range(dout))
            if mirrored:
                notes.append("the table is graded-antisymmetric" if sign == -1
                             else "the table is graded-commutative (symmetric)")
    return notes


def structural_note_line(kind, degrees, dims, constants, prime):
    """The structural notes of one table as a single caption sentence, or ``""``."""
    notes = structural_notes(kind, degrees, dims, constants, prime)
    if not notes:
        return ""
    s = "; ".join(notes)
    return s[0].upper() + s[1:] + "."


def cayley_table(kind, degrees, out_degree, dims, constants, prime=None):
    """The structured Cayley multiplication table of one product bidegree: row labels
    (left classes ``α^p_i``), column labels (right classes ``α^q_j`` / ``z^n_j`` for
    cap), the corner operator, and one TeX cell per entry -- the product in the target
    basis (``α^out`` / ``z^out``), balanced-rep coefficients. Presentation-agnostic:
    the HTML surfaces render it via ``render_html.cayley_grid_html``, the GUI mirrors
    it. ``structural_note_line`` is carried alongside so a caller shows it above the
    grid."""
    dl, dr, dout = dims
    left_deg, right_deg = degrees
    row_labels = [r"%s^{%d}_{%d}" % (_LEFT, left_deg, i + 1) for i in range(dl)]
    col_labels = [r"%s^{%d}_{%d}" % (_RIGHT[kind], right_deg, j + 1) for j in range(dr)]
    cells = [[cell_tex(kind, out_degree,
                       [constants[k][i][j] for k in range(dout)], prime)
              for j in range(dr)] for i in range(dl)]
    return {"corner": _CORNER[kind], "row_labels": row_labels,
            "col_labels": col_labels, "cells": cells, "dl": dl, "dr": dr,
            "note": structural_note_line(kind, degrees, dims, constants, prime)}


# --------------------------------------------------------------------------- #
# ONE BIG Cayley table per family (Marco 2026-08-01 addendum). Rows/columns run over
# ALL (co)homology basis classes, degree-major, the degree read from the class'
# superscript; a heavier rule marks each degree boundary. A cell whose target degree
# lies beyond the computed window is an em dash (NOT computed, NOT zero); a computed
# vanishing product is 0. When a family's per-axis class count exceeds the display
# cap, the surface falls back to the per-bidegree ``cayley_table`` grids.
# --------------------------------------------------------------------------- #

CAYLEY_AXIS_CAP = 50            # per-axis class count above which the big form is dropped
EM_DASH = "—"


def _combined_out_degree(kind, p, q):
    if kind == "cup":
        return p + q
    if kind == "bracket":
        return p + q - 1
    return q - p                # cap: (p, n) -> n - p


def beyond_window_note():
    """The legend for the em-dash honesty mark used by the big product table."""
    return (EM_DASH + " marks a cell whose target degree lies beyond the computed "
            "window (not computed); a computed vanishing product is shown as 0.")


def _combined_note(kind, tbl, row_meta, prime):
    """The structural caption derived over the WHOLE in-window region of a big cup/
    bracket table (row and column axes coincide): "all squares are 0" when every
    COMPUTED diagonal product vanishes, and the graded (anti)symmetry law when every
    both-computed transpose pair satisfies it with its per-block sign. ``""`` for cap
    (rows are cohomology, columns homology -- no transpose) or when nothing is
    certifiable."""
    if kind not in ("cup", "bracket"):
        return ""
    notes = []
    diag_seen = diag_zero = 0
    for gi, (p, i) in enumerate(row_meta):
        t = tbl.get((p, p))
        if t is None:
            continue
        diag_seen += 1
        dout = t["dims"][2]
        if all(str(t["constants"][k][i][i]) == "0" for k in range(dout)):
            diag_zero += 1
    if diag_seen and diag_zero == diag_seen:
        notes.append("all squares are 0")
    if prime is not None:
        seen = ok = 0
        all_int = True
        for gi, (p, i) in enumerate(row_meta):
            for gj, (q, j) in enumerate(row_meta):
                t, tT = tbl.get((p, q)), tbl.get((q, p))
                if t is None or tT is None:
                    continue
                dout = t["dims"][2]
                sign = _mirror_sign(kind, p, q)
                pairs = [(t["constants"][k][i][j], tT["constants"][k][j][i])
                         for k in range(dout)]
                if not all(_is_int(a) and _is_int(b) for a, b in pairs):
                    all_int = False           # a non-residue entry: cannot certify
                    break
                good = all((int(a) - sign * int(b)) % prime == 0 for a, b in pairs)
                seen += 1
                ok += 1 if good else 0
            if not all_int:
                break
        if all_int and seen and ok == seen:
            notes.append("the cup product is graded-commutative" if kind == "cup"
                         else "the Gerstenhaber bracket is graded-antisymmetric")
    if not notes:
        return ""
    s = "; ".join(notes)
    return s[0].upper() + s[1:] + "."


def combined_cayley(kind, tables, prime=None):
    """The ONE big Cayley table of a whole product family, from the per-bidegree
    ``tables`` (each a mapping with ``degrees`` / ``out_degree`` / ``dims`` /
    ``constants``). Rows run over every left (cohomology) class degree-major; columns
    over every right class (cohomology for cup/bracket, homology for cap) degree-major.

    Returns ``{"over_cap": True, "rows", "cols"}`` when either axis exceeds
    :data:`CAYLEY_AXIS_CAP` (the caller then renders per-bidegree grids); otherwise a
    struct ``{corner, row_labels, col_labels, row_degsep, col_degsep, cells, dl, dr,
    note, has_beyond}`` where each cell is ``"0"`` / a signed combination / the em dash
    (beyond the computed window). ``row_degsep[i]`` / ``col_degsep[j]`` flag the first
    class of a new degree block (skipping the very first) for the heavier grid rule."""
    tbl, left_dims, right_dims = {}, {}, {}
    for t in tables:
        p, q = t["degrees"]
        dl, dr, _ = t["dims"]
        tbl[(p, q)] = {"out_degree": t["out_degree"], "dims": list(t["dims"]),
                       "constants": t["constants"]}
        left_dims[p] = dl
        right_dims[q] = dr
    left_degs, right_degs = sorted(left_dims), sorted(right_dims)
    total_rows = sum(left_dims[p] for p in left_degs)
    total_cols = sum(right_dims[q] for q in right_degs)
    if total_rows >= CAYLEY_AXIS_CAP or total_cols >= CAYLEY_AXIS_CAP:
        return {"over_cap": True, "rows": total_rows, "cols": total_cols}
    top = max(left_degs + right_degs) if (left_degs and right_degs) else 0

    row_labels, row_degsep, row_meta = [], [], []
    for bi, p in enumerate(left_degs):
        for i in range(left_dims[p]):
            row_labels.append(r"%s^{%d}_{%d}" % (_LEFT, p, i + 1))
            row_degsep.append(i == 0 and bi > 0)
            row_meta.append((p, i))
    col_labels, col_degsep, col_meta = [], [], []
    for bj, q in enumerate(right_degs):
        for j in range(right_dims[q]):
            col_labels.append(r"%s^{%d}_{%d}" % (_RIGHT[kind], q, j + 1))
            col_degsep.append(j == 0 and bj > 0)
            col_meta.append((q, j))

    cells, has_beyond = [], False
    for (p, i) in row_meta:
        row = []
        for (q, j) in col_meta:
            target = _combined_out_degree(kind, p, q)
            if target < 0:                          # cap below degree 0: structural zero
                row.append("0")
                continue
            t = tbl.get((p, q))
            if t is not None:
                dout = t["dims"][2]
                row.append(cell_tex(kind, target,
                                    [t["constants"][k][i][j] for k in range(dout)],
                                    prime))
            else:                                   # target > top: not computed
                has_beyond = True
                row.append(EM_DASH)
        cells.append(row)

    return {"over_cap": False, "corner": _CORNER[kind],
            "row_labels": row_labels, "col_labels": col_labels,
            "row_degsep": row_degsep, "col_degsep": col_degsep,
            "cells": cells, "dl": total_rows, "dr": total_cols,
            "has_beyond": has_beyond,
            "note": _combined_note(kind, tbl, row_meta, prime)}


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
