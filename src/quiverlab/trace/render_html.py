"""HTML worked-steps renderer -- self-contained, NO JavaScript, NO external
resources (no CDN ``<script>``, no ``<link>``, no network fetch), so the file
renders identically with the network off. Math is TYPESET as MathML (rendered
natively by every current browser, print-to-PDF included) with the original LaTeX
kept verbatim in a ``<annotation encoding="application/x-tex">`` -- so the typeset
display and the copy/paste-able source coexist in one self-contained file. The page
is print-ready (print CSS: ``@page`` margins, page-break control, a screen-only
"print to PDF" hint): the GUI opens it in a tab and calls ``window.print()``, and a
downloaded ``report.html`` prints to PDF from any browser. Float-free: all numbers
come from event fields (ints/strings).

Matrices are displayed as INDEXED GRIDS (`matrix_grid`): an HTML table with a
header row of column indices, a header column of row indices and a light-grey rule,
so an entry can be read off by position (Marco 2026-07-29). `_pmatrix` remains for
the TeX source form.

Shared helpers: `derive_dims` + `_dims_kind` are REUSED from Task 9's render_text.
The matrix->pmatrix TeX-source helper `_pmatrix` and the small TeX-source math
builders (`oplus_tex`, `factor_stack_tex`, `_dimvec_tex`, `_tex_escape`) live HERE
-- this module is the sole owner of the worked-steps math source. The TeX source
those helpers produce is carried through into the MathML ``<annotation>`` (and
typeset by the converter below), so the display and the copy/paste-able source stay
in lock-step."""
import re

from quiverlab.trace.events import (
    Dispatch, ResolutionTerm, RankStep, ModuleTerm, ModuleDifferential,
    ExtDegree, ExtReps, StepNote, ResultDims, ProductStep, ProductBasis,
)
from quiverlab.trace.render_text import (
    derive_dims, _dims_kind, compute_algebra_objects, ext_result_runs,
    ELISION_PREAMBLE,
)

# Inline-only styling (no external stylesheet); keeps the file fully self-contained
# AND print-ready. Light-forced (a report is printed on white); the screen-only
# ".ql-hint" banner is hidden by the print stylesheet.
_STYLE = (
    "<style>"
    "html{background:#fff;color:#111}"
    "body{font-family:Georgia,'Times New Roman',serif;line-height:1.5;"
    "max-width:46em;margin:2em auto;padding:0 1.5em;background:#fff;color:#111}"
    "h1{font-size:1.6em;border-bottom:2px solid #333;padding-bottom:.25em}"
    "h2{font-size:1.15em;margin-top:1.5em;border-bottom:1px solid #ccc;"
    "padding-bottom:.12em}"
    "p{margin:.5em 0}"
    "math{font-size:1.05em}"
    # Equations are shown COMPLETE -- no clip, no scrollbar (Marco 2026-07-29).
    # A wide matrix is instead typeset a size down (see _fit_pct), which is
    # deterministic and hides nothing.
    ".ql-eq{margin:.7em 0}"
    ".ql-hint{background:#eef3ff;border:1px solid #9db8e8;border-radius:6px;"
    "padding:.6em .85em;margin:0 0 1.4em;font-family:sans-serif;font-size:.92em;"
    "color:#123}"
    ".ql-note{color:#444;font-size:.95em}"
    # .ql-dims is the DEGREE table (n | value_n); .ql-table is every other
    # structural table (rad/top/soc, summands, resolution terms).
    "table.ql-dims,table.ql-table{border-collapse:collapse;margin:.7em 0;"
    "font-size:.95em}"
    "table.ql-dims td,table.ql-dims th,table.ql-table td,table.ql-table th"
    "{border:1px solid #bbb;padding:2px 10px;text-align:center}"
    "table.ql-dims th,table.ql-table th{background:#f2f2f2}"
    # An indexed matrix grid: white cells on a light-grey rule, with the row and
    # column indices in grey headers so an entry can be read off by position
    # (Marco 2026-07-29). DOUBLE ZEBRA STRIPING (Marco 2026-08-01): alternate the
    # DATA rows (even rows light grey) AND the DATA columns (even columns a
    # translucent overlay), so intersections give four subtly distinct levels and
    # an entry is easy to locate by its (row, column) index. Striping is pure CSS
    # nth-child on the existing <td> cells -- no attribute or structure change, so
    # the tests that read entries out of the grid are unaffected. Header row/column
    # (<th>) keep their distinct grey and are never striped. print-color-adjust:exact
    # keeps the shading when the report is printed (else browsers drop backgrounds).
    "table.ql-matrix{border-collapse:collapse;margin:.5em 0 .9em;"
    "font-variant-numeric:tabular-nums}"
    "table.ql-matrix td,table.ql-matrix th{border:1px solid #d0d0d0;"
    "padding:2px 9px;text-align:right;background:#fff;"
    "-webkit-print-color-adjust:exact;print-color-adjust:exact}"
    # data rows: even (data) rows light grey -- the header is tr:nth-child(1), so
    # tr:nth-child(even) selects the 1st, 3rd, ... data row.
    "table.ql-matrix tr:nth-child(even) td{background-color:#f2f2f2}"
    # data columns: even (data) columns a translucent overlay; background-image
    # layers OVER the row background-color, so the two stripings stack.
    "table.ql-matrix td:nth-child(even)"
    "{background-image:linear-gradient(rgba(0,0,0,.045),rgba(0,0,0,.045))}"
    "table.ql-matrix th{background:#f0f0f0;color:#555;font-weight:normal;"
    "font-size:.85em;text-align:center}"
    "table.ql-matrix th.ql-corner{background:#e4e4e4;border-color:#c4c4c4}"
    # Big-family Cayley table (Marco 2026-08-01): a heavier rule at each DEGREE
    # boundary so the degree blocks read apart (CSS-only; the structural read-back is
    # unaffected). ql-degrow on the first row of a new left-degree block draws its top
    # edge; ql-degcol on the first column of a new right-degree block draws its left
    # edge. Applied to both the header (<th>) and the body (<td>) cells of that
    # boundary so the rule is continuous.
    "table.ql-cayley td.ql-degrow,table.ql-cayley th.ql-degrow"
    "{border-top:2px solid #8a8a8a}"
    "table.ql-cayley td.ql-degcol,table.ql-cayley th.ql-degcol"
    "{border-left:2px solid #8a8a8a}"
    ".ql-mlabel{margin:.6em 0 .1em}"
    "ol,ul{padding-left:1.4em}"
    # Fine-grained table of contents (Marco 2026-07-31): top-level sections with
    # per-degree / per-theory / per-product subsections nested beneath.
    "ol.ql-toc,ol.ql-toc-sub{list-style:none}"
    "ol.ql-toc>li{margin:.25em 0}"
    "ol.ql-toc-sub{margin:.1em 0 .3em;font-size:.95em;color:#333}"
    "ol.ql-toc a{text-decoration:none}"
    # The JSON-guide appendix table.
    "table.ql-guide{border-collapse:collapse;margin:.7em 0;font-size:.92em}"
    "table.ql-guide td,table.ql-guide th{border:1px solid #bbb;padding:3px 9px;"
    "text-align:left;vertical-align:top}"
    "table.ql-guide th{background:#f2f2f2}"
    "table.ql-guide code{font-size:.95em}"
    # Plan 35 UNIT 2: the ordered (co)chain enumeration and the explicit basis
    # classes -- compact lists so the per-degree reps read as one block.
    ".ql-enum,.ql-classes{margin:.2em 0 .6em}"
    ".ql-enum li,.ql-classes li{margin:.05em 0}"
    "h4{font-size:1.02em;margin:1em 0 .2em;color:#222}"
    "@media print{.ql-hint{display:none}"
    "body{max-width:none;margin:0;font-size:11pt}"
    "h2,.ql-eq,math,table.ql-dims,table.ql-table,table.ql-matrix"
    "{break-inside:avoid;page-break-inside:avoid}"
    "h1,h2{break-after:avoid}}"
    "@page{margin:2cm}"
    "</style>")

# A matrix wider than this many columns is typeset one size down so it still fits
# the text column WITHOUT a scrollbar; the shrink is proportional and floored, so a
# very wide matrix is small but never clipped and never hidden. Integer arithmetic
# only (the src/ no-floats gate).
_EQ_FIT_COLS = 16
_EQ_MIN_PCT = 50


def _fit_pct(ncols):
    """Font-size percentage for an equation whose widest matrix has ``ncols``
    columns, or None when it fits at full size (shrink-only, never magnifies)."""
    if not ncols or ncols <= _EQ_FIT_COLS:
        return None
    return max(_EQ_MIN_PCT, _EQ_FIT_COLS * 100 // ncols)

# --------------------------------------------------------------------------- #
# LaTeX -> presentation MathML converter (self-contained, no JS).
#
# The trace renderers emit a small, CLOSED grammar (pmatrix/matrix environments,
# sub/superscripts, a fixed handful of commands, spacing, atoms). This converts
# exactly THAT grammar to presentation MathML. Anything OUTSIDE the covered grammar
# -- an unknown command or environment -- raises `_UnknownTeX`, and `_math` catches
# it and routes the WHOLE math run to a `<mtext>` of the escaped LaTeX source (the
# real, advertised fallback: the source is shown verbatim, WITH its backslashes,
# never a backslash-stripped command name masquerading as typeset math). So the
# document is always well-formed and never shows a garbled glyph, and any command
# the LaTeX renderer starts emitting that we have not taught the converter degrades
# HONESTLY (readable source) instead of leaking (e.g. `\frac`->"frac 1 2"). The
# fidelity guard test (tests/trace/test_render_html_fidelity_p34.py) pins this:
# every real module+HH sample converts with no leaked command name.
# Deterministic: a pure function of the input string.
# --------------------------------------------------------------------------- #


class _UnknownTeX(Exception):
    """A command/environment outside the covered grammar. Caught by `_math`, which
    then renders the whole run as `<mtext>` of the escaped source (the real fallback)."""


_ROWBREAK = "\\\\"          # the TeX row separator: two backslash characters

_TOK = re.compile(
    r"\\begin\{[A-Za-z*]+\}"    # \begin{env}
    r"|\\end\{[A-Za-z*]+\}"     # \end{env}
    r"|\\[A-Za-z]+"             # \command
    r"|\\."                     # \\  \  \; \, and other escaped symbols
    r"|[A-Za-z]+"              # identifier run
    r"|[0-9]+"                # number run
    r"|\s+"                  # whitespace
    r"|.",                   # any single character
    re.DOTALL)

# \command -> a MathML <mi> identifier (upright when >1 char, per spec).
_CMD_MI = {
    r"\dim": "dim", r"\Hom": "Hom", r"\ker": "ker", r"\coker": "coker",
    # lowercase Greek (incl. the var- forms the module symbols use, e.g. \varphi
    # from the decompose splitting endomorphism, \varepsilon from the augmentation)
    r"\alpha": "&#945;", r"\beta": "&#946;", r"\gamma": "&#947;",
    r"\delta": "&#948;", r"\epsilon": "&#949;", r"\varepsilon": "&#949;",
    r"\zeta": "&#950;", r"\eta": "&#951;", r"\theta": "&#952;",
    r"\vartheta": "&#977;", r"\iota": "&#953;", r"\kappa": "&#954;",
    r"\lambda": "&#955;", r"\mu": "&#956;", r"\nu": "&#957;", r"\xi": "&#958;",
    r"\pi": "&#960;", r"\varpi": "&#982;", r"\rho": "&#961;",
    r"\varrho": "&#961;", r"\sigma": "&#963;", r"\varsigma": "&#962;",
    r"\tau": "&#964;", r"\upsilon": "&#965;", r"\phi": "&#966;",
    r"\varphi": "&#966;", r"\chi": "&#967;", r"\psi": "&#968;", r"\omega": "&#969;",
    # uppercase Greek
    r"\Gamma": "&#915;", r"\Delta": "&#916;", r"\Theta": "&#920;",
    r"\Lambda": "&#923;", r"\Xi": "&#926;", r"\Pi": "&#928;",
    r"\Sigma": "&#931;", r"\Phi": "&#934;", r"\Psi": "&#936;", r"\Omega": "&#937;",
    r"\partial": "&#8706;", r"\infty": "&#8734;",
}
# Matrix-like environments the parser understands, and their fence delimiters. An
# environment OUTSIDE this set raises `_UnknownTeX` (honest fallback). `array` carries
# a ``{colspec}`` argument (consumed); ``matrix``/``smallmatrix`` render undelimited.
_MATRIX_ENVS = frozenset({
    "matrix", "smallmatrix", "pmatrix", "bmatrix", "vmatrix", "Vmatrix",
    "array", "cases",
})
_MATRIX_DELIMS = {
    "pmatrix": ("(", ")"), "bmatrix": ("[", "]"),
    "vmatrix": ("|", "|"), "Vmatrix": ("&#8214;", "&#8214;"),
    "cases": ("{", ""),
}
# \command -> a MathML <mo> operator (entity glyphs).
_CMD_MO = {
    r"\oplus": "&#8853;", r"\otimes": "&#8855;", r"\to": "&#8594;",
    r"\rightarrow": "&#8594;", r"\mapsto": "&#8614;", r"\cdot": "&#8901;",
    r"\times": "&#215;", r"\cong": "&#8773;", r"\leq": "&#8804;",
    r"\geq": "&#8805;",
    # Plan 35 product-chapter operators: cup, cap, the Gerstenhaber circle, the
    # cyclic-sum ellipsis and sigma.
    r"\cup": "&#8746;", r"\cap": "&#8745;", r"\circ": "&#8728;",
    r"\cdots": "&#8943;", r"\sum": "&#8721;", r"\smile": "&#8994;",
    r"\frown": "&#8995;",
}
# \command -> an <mspace> of the given width.
_CMD_SPACE = {
    r"\qquad": "2em", r"\quad": "1em", r"\;": "0.278em", r"\:": "0.222em",
    r"\,": "0.167em", "\\ ": "0.25em",
}
# Size/layout commands with no MathML analogue -- consumed, emit nothing.
_CMD_IGNORE = frozenset({
    r"\big", r"\bigl", r"\bigr", r"\Big", r"\Bigl", r"\Bigr",
    r"\left", r"\right", r"\displaystyle", r"\noindent", r"\!",
})


def _skip_ws(tokens, i):
    while i < len(tokens) and tokens[i].isspace():
        i += 1
    return i


def _read_raw_group(tokens, i):
    """Return (verbatim-text, next-i) for the ``{...}`` at ``i`` (brace-matched),
    or the single following token when there is no group. Used for the literal
    arguments of ``\\text{}`` / ``\\operatorname{}``."""
    i = _skip_ws(tokens, i)
    if i >= len(tokens) or tokens[i] != "{":
        return (tokens[i], i + 1) if i < len(tokens) else ("", i)
    i += 1
    depth, buf = 1, []
    while i < len(tokens) and depth > 0:
        t = tokens[i]
        if t == "{":
            depth += 1
            buf.append(t)
        elif t == "}":
            depth -= 1
            if depth == 0:
                return "".join(buf), i + 1
            buf.append(t)
        else:
            buf.append(t)
        i += 1
    return "".join(buf), i


def _symbol_atom(ch):
    if ch.isspace():
        return ""
    return "<mo>%s</mo>" % _esc(ch)


def _parse_matrix(tokens, i, endtok):
    """Parse a matrix body (cells split on ``&``, rows on ``\\\\``) into an
    <mtable>, consuming the closing ``\\end{...}``."""
    rows, cur = [], []
    while i < len(tokens):
        cell, i = _parse_seq(tokens, i, {"&", _ROWBREAK, endtok})
        cur.append(cell)
        if i >= len(tokens):
            break
        sep = tokens[i]
        i += 1
        if sep == "&":
            continue
        if sep == _ROWBREAK:
            rows.append(cur)
            cur = []
        else:                                    # endtok
            break
    if cur:
        rows.append(cur)
    # A trailing "\\" leaves an empty final row -- drop it.
    if rows and len(rows[-1]) == 1 and rows[-1][0] == "":
        rows.pop()
    body = "".join(
        "<mtr>%s</mtr>" % "".join("<mtd>%s</mtd>" % c for c in r) for r in rows)
    return "<mtable>%s</mtable>" % body, i


def _parse_command(tokens, i, stops):
    cmd = tokens[i]
    i += 1
    if cmd == r"\operatorname":
        content, i = _read_raw_group(tokens, i)
        return "<mi>%s</mi>" % _esc(content), i
    if cmd in (r"\text", r"\mathrm", r"\textit", r"\mathit", r"\mathbf"):
        content, i = _read_raw_group(tokens, i)
        return "<mtext>%s</mtext>" % _esc(content), i
    if cmd == r"\underline":
        inner, i = _parse_atom(tokens, i, stops)
        return ('<munder accentunder="true"><mrow>%s</mrow><mo>&#8213;</mo>'
                '</munder>' % inner), i
    if cmd in _CMD_MO:
        return "<mo>%s</mo>" % _CMD_MO[cmd], i
    if cmd in _CMD_MI:
        return "<mi>%s</mi>" % _CMD_MI[cmd], i
    if cmd in _CMD_SPACE:
        return '<mspace width="%s"></mspace>' % _CMD_SPACE[cmd], i
    if cmd in _CMD_IGNORE:
        return "", i
    if cmd == r"\resizebox":
        # \resizebox{width}{height}{content}: a graphicx sizing wrapper with no MathML
        # analogue -- discard the two size args and typeset only the content (a wide
        # matrix wrapped in \resizebox in the embedded source is unwrapped here).
        _, i = _read_raw_group(tokens, i)
        _, i = _read_raw_group(tokens, i)
        return _parse_atom(tokens, i, stops)
    # Outside the covered grammar: raise so `_math` routes the WHOLE run to the
    # <mtext>-of-escaped-source fallback -- NEVER a backslash-stripped command name
    # typeset as if it were math (the old `\frac`->"frac" leak).
    raise _UnknownTeX(cmd)


def _parse_atom(tokens, i, stops):
    i = _skip_ws(tokens, i)
    if i >= len(tokens):
        return "", i
    t = tokens[i]
    if t == "{":
        inner, i = _parse_seq(tokens, i + 1, {"}"})
        if i < len(tokens) and tokens[i] == "}":
            i += 1
        return "<mrow>%s</mrow>" % inner, i
    if t.startswith(r"\begin{"):
        env = t[len(r"\begin{"):-1]
        if env not in _MATRIX_ENVS:
            raise _UnknownTeX(t)                  # honest fallback for an unknown env
        i += 1
        if env == "array":
            _, i = _read_raw_group(tokens, i)     # consume the {colspec} argument
        table, i = _parse_matrix(tokens, i, r"\end{%s}" % env)
        delim = _MATRIX_DELIMS.get(env)
        if delim:
            close = "<mo>%s</mo>" % delim[1] if delim[1] else ""
            return "<mrow><mo>%s</mo>%s%s</mrow>" % (delim[0], table, close), i
        return table, i                          # matrix/smallmatrix: no delimiters
    if t.startswith("\\"):
        return _parse_command(tokens, i, stops)
    if t[:1].isalpha():
        return "<mi>%s</mi>" % _esc(t), i + 1
    if t[:1].isdigit():
        return "<mn>%s</mn>" % _esc(t), i + 1
    return _symbol_atom(t), i + 1


def _parse_element(tokens, i, stops):
    """One atom plus any trailing ``_``/``^`` scripts, as msub/msup/msubsup."""
    base, i = _parse_atom(tokens, i, stops)
    sub = sup = None
    while True:
        i = _skip_ws(tokens, i)
        if i < len(tokens) and tokens[i] == "_":
            sub, i = _parse_atom(tokens, i + 1, stops)
        elif i < len(tokens) and tokens[i] == "^":
            sup, i = _parse_atom(tokens, i + 1, stops)
        else:
            break
    if sub is not None and sup is not None:
        return "<msubsup>%s%s%s</msubsup>" % (base, sub, sup), i
    if sub is not None:
        return "<msub>%s%s</msub>" % (base, sub), i
    if sup is not None:
        return "<msup>%s%s</msup>" % (base, sup), i
    return base, i


def _parse_seq(tokens, i, stops):
    parts = []
    while i < len(tokens):
        i = _skip_ws(tokens, i)
        if i >= len(tokens) or tokens[i] in stops:
            break
        node, i = _parse_element(tokens, i, stops)
        parts.append(node)
    return "".join(parts), i


def _tex_to_mathml_body(expr):
    tokens = _TOK.findall(expr)
    body, _ = _parse_seq(tokens, 0, frozenset())
    return body


def _math(expr, ncols=None):
    """A display equation: typeset presentation MathML with the LaTeX source kept
    verbatim in an x-tex ``<annotation>`` (self-contained, no JavaScript). On any
    conversion surprise, fall back to ``<mtext>`` of the escaped source so the
    document stays well-formed and the source is still readable.

    ``ncols`` is the widest embedded matrix's column count, when the caller knows
    it: a wide matrix is typeset a size down so the equation fits the text column
    complete, with no scrollbar and nothing clipped (Marco 2026-07-29)."""
    try:
        pres = _tex_to_mathml_body(expr) or ("<mtext>%s</mtext>" % _esc(expr))
    except Exception:
        pres = "<mtext>%s</mtext>" % _esc(expr)
    pct = _fit_pct(ncols)
    style = ' style="font-size:%d%%"' % pct if pct is not None else ""
    return ('<div class="ql-eq"><math display="block"%s><semantics><mrow>%s</mrow>'
            '<annotation encoding="application/x-tex">%s</annotation>'
            '</semantics></math></div>' % (style, pres, _esc(expr)))


def _math_inline(expr):
    """The same MathML conversion as :func:`_math`, inline and WITHOUT the
    ``div.ql-eq`` wrapper -- for math inside a table cell."""
    try:
        pres = _tex_to_mathml_body(expr) or ("<mtext>%s</mtext>" % _esc(expr))
    except Exception:
        pres = "<mtext>%s</mtext>" % _esc(expr)
    return ('<math display="inline"><semantics><mrow>%s</mrow>'
            '<annotation encoding="application/x-tex">%s</annotation>'
            '</semantics></math>' % (pres, _esc(expr)))


def matrix_grid(matrix, label=None):
    """A matrix as an INDEXED GRID (Marco 2026-07-29): a header row of column
    indices, a header column of row indices, and a light rule between every cell,
    so an entry can be read off by its position. 1-based, the mathematician's
    convention. ``label`` is optional TeX shown as ``label =`` above the grid.

    Returns the HTML for the whole block (caption + grid). A zero-DIMENSIONAL
    matrix (no rows or no columns) renders as the symbol ``0`` -- an empty grid
    would be a stray box. Entries are copied verbatim (exact ints / fraction
    strings), never reformatted."""
    rows = matrix or []
    ncols = len(rows[0]) if rows and rows[0] is not None else 0
    if not rows or not ncols:
        return _math("%s = 0" % label) if label else _math("0")
    # Marco 2026-08-02: only matrices with fewer than DISPLAY_CAP (=20) rows AND columns
    # are shown; a larger one states its size and points at the JSON (the complete matrix
    # is always in result.json / trace.json). This is the single chokepoint for every
    # ORDINARY matrix grid. Product Cayley tables do NOT route through here -- they render
    # via cayley_grid_html and are deliberately uncapped.
    if len(rows) >= DISPLAY_CAP or ncols >= DISPLAY_CAP:
        pre = "%s = " % _esc(label) if label else ""
        return ("<p class='ql-note'>%s%d×%d matrix (%d rows and %d columns exceed the "
                "%d-line display cap); the complete matrix is in the accompanying "
                "JSON record.</p>" % (pre, len(rows), ncols, len(rows), ncols,
                                      DISPLAY_CAP))
    out = []
    if label:
        out.append('<p class="ql-mlabel">%s =</p>' % _math_inline(label))
    head = ['<th class="ql-corner"></th>']
    head += ["<th>%d</th>" % (j + 1) for j in range(ncols)]
    body = []
    for i, row in enumerate(rows):
        cells = "".join("<td>%s</td>" % _esc(str(x)) for x in row)
        body.append("<tr><th>%d</th>%s</tr>" % (i + 1, cells))
    out.append('<table class="ql-matrix"><tr>%s</tr>%s</table>'
               % ("".join(head), "".join(body)))
    return "".join(out)


def _cayley_cell_html(c, classes):
    """One Cayley body cell: ``0`` / em dash rendered as plain text, a combination
    typeset; the degree-boundary classes (CSS-only heavier rule) are attached to the
    ``<td>`` -- the read-back helper tolerates the attribute."""
    attr = ' class="%s"' % " ".join(classes) if classes else ""
    inner = c if c in ("0", "—") else _math_inline(c)
    return "<td%s>%s</td>" % (attr, inner)


def cayley_grid_html(table):
    """Render a structured Cayley multiplication table
    (:func:`quiverlab.trace.products.cayley_table` per-bidegree, or
    :func:`quiverlab.trace.products.combined_cayley` for the whole family) as an INDEXED
    GRID, reusing the ``matrix_grid`` conventions: the corner cell holds the product
    operator, the header ROW the right-factor classes, the header COLUMN the left-factor
    classes, and each body cell the product in the target basis (``0`` / a signed
    combination / an em dash beyond the computed window). Double zebra striping + print
    colours come free from the shared ``ql-matrix`` CSS; ``ql-cayley`` marks it for the
    read-back helper, and the optional ``row_degsep`` / ``col_degsep`` flags add a
    heavier CSS-only rule (``ql-degrow`` / ``ql-degcol``) at each degree boundary.

    Marco 2026-08-02: the product Cayley table is UNCAPPED -- it always renders in full,
    however many classes the axes carry (the deliberate exception to
    :data:`DISPLAY_CAP`, which bounds every ordinary matrix grid via
    :func:`matrix_grid`). A ``ql-cayley`` grid never routes through that cap."""
    dl, dr = table["dl"], table["dr"]
    if not dl or not dr:
        return _math("0")
    rds = table.get("row_degsep") or [False] * dl
    cds = table.get("col_degsep") or [False] * dr
    head = ['<th class="ql-corner">%s</th>' % _math_inline(table["corner"])]
    for j, lbl in enumerate(table["col_labels"]):
        cls = ' class="ql-degcol"' if cds[j] else ""
        head.append("<th%s>%s</th>" % (cls, _math_inline(lbl)))
    body = []
    for i, row in enumerate(table["cells"]):
        rcls = ' class="ql-degrow"' if rds[i] else ""
        cells = ["<th%s>%s</th>" % (rcls, _math_inline(table["row_labels"][i]))]
        for j, c in enumerate(row):
            classes = []
            if rds[i]:
                classes.append("ql-degrow")
            if cds[j]:
                classes.append("ql-degcol")
            cells.append(_cayley_cell_html(c, classes))
        body.append("<tr>%s</tr>" % "".join(cells))
    return ('<table class="ql-matrix ql-cayley"><tr>%s</tr>%s</table>'
            % ("".join(head), "".join(body)))


def _event_grid(ev, label=None):
    """:func:`matrix_grid` for any event carrying ``(elided, matrix, nrows, ncols,
    note)`` -- RankStep, ModuleDifferential, ExtDegree. An elided body is stated as
    its shape note instead of a fabricated grid."""
    if ev.elided or ev.matrix is None:
        note = "%s: %s" % (label, ev.note) if label else ev.note
        return "<p class='ql-note'>%s</p>" % _esc(note)
    if ev.nrows == 0 or ev.ncols == 0:
        return _math("%s = 0" % label) if label else _math("0")
    return matrix_grid([[ev.matrix[i][j] for j in range(ev.ncols)]
                        for i in range(ev.nrows)], label=label)


def _dims_table(row_label, dims, col_label="n"):
    """A degree table ``n | value_n`` -- the shape the GUI uses for HH / Ext / Tor
    (Marco 2026-07-29: the Result used to be one long equation in a scroll box)."""
    head = ["<th>%s</th>" % _esc(col_label)]
    body = ["<th>%s</th>" % _esc(row_label)]
    for i, d in enumerate(dims):
        head.append("<td>%d</td>" % i)
        body.append("<td>%s</td>" % _esc(str(d)))
    return ('<table class="ql-dims"><tr>%s</tr><tr>%s</tr></table>'
            % ("".join(head), "".join(body)))


class _MatrixEcho:
    """Remembers matrices already printed so an IDENTICAL later differential is
    referenced instead of repeated (Marco 2026-07-29). Decisive for a periodic
    resolution, where every second differential is the same matrix.

    ``label_for(event, symbol)`` returns the earlier symbol when this exact matrix
    has been shown, else None (and records it). An elided matrix is never matched:
    its body was dropped, so equality is unknowable."""

    def __init__(self):
        self._seen = {}

    def label_for(self, ev, symbol):
        if getattr(ev, "elided", False) or getattr(ev, "matrix", None) is None:
            return None
        key = (ev.nrows, ev.ncols,
               tuple(tuple(row) for row in ev.matrix))
        if key in self._seen:
            return self._seen[key]
        self._seen[key] = symbol
        return None


def _hh_row_label(kind):
    """Row label for the (co)homology table: ``HH^`` -> ``dim HH^n``, ``HH_`` -> ``dim HH_n``."""
    return "dim %sn" % kind


def _hh_heading(kind):
    """Section heading naming WHAT the table is. "Result" said nothing -- every
    section of the page is a result (Marco 2026-07-29)."""
    return "Hochschild cohomology" if kind == "HH^" else "Hochschild homology"


def hh_typing_html(kind, route):
    """The typing paragraph for a Hochschild (co)homology section (Marco 2026-07-31):
    exactly what the engine computes and what the bar-bracket/tensor notation means.
    ``kind`` is HH^ / HH_ (or hh_cohomology / hh_homology); ``route`` is "bar"/"cs"."""
    from quiverlab.trace.interpretations import hh_space_typing
    return "<p class='ql-note'>%s</p>" % _esc(hh_space_typing(kind, route))


def route_of_engine(engine):
    """"cs" when the recorded engine/basis string names Chouhy-Solotar, else "bar"
    (the bar / fast / minimal cochain bases all share the bar cochain typing)."""
    r = (str(engine) or "").lower()
    return "cs" if ("cs" in r or "chouhy" in r or "solotar" in r) else "bar"


def _route_of_dispatch(used):
    for e in used or ():
        if route_of_engine(getattr(e, "route", "")) == "cs":
            return "cs"
    return "bar"


def _results_carry_hh(results, kind):
    """True when the Computed results section already shows this exact HH table, so
    printing it again under its own heading would just repeat the same numbers."""
    if not results or not kind:
        return False
    from quiverlab.trace.results_html import normalize
    want = "hh_cohomology" if kind == "HH^" else "hh_homology"
    return any(k == want for k, _ in normalize(results))


def _term_summands_html(term, degree):
    """The resolution term NAMED as a direct sum of projective bimodules, when the
    engine recorded the generators' corners (Marco 2026-07-29 -- the report used to
    give only a generator count). ``P(v,w)`` is ``A e_v (x) e_w A``; the standard
    Chouhy-Solotar term is ``C_n = (+)_{s in S_n} A e_{o(s)} (x) e_{t(s)} A``.
    Returns '' when the corners were not recorded (e.g. the bar resolution over a
    structure-constants algebra, which is not vertex-decomposed)."""
    corners = getattr(term, "corners", None)
    if not corners:
        return ""
    counts = {}
    for pair in corners:
        key = (str(pair[0]), str(pair[1]))
        counts[key] = counts.get(key, 0) + 1
    parts = []
    for (v, w) in sorted(counts):
        base = "P(%s,%s)" % (v, w)
        c = counts[(v, w)]
        parts.append(base if c == 1 else "%s^{%d}" % (base, c))
    tex = r" \oplus ".join(parts)
    return _math(r"C_{%d} = %s \qquad P(v,w) = A e_{v} \otimes e_{w} A"
                 % (degree, tex))


def _pi_section_html(objects, note):
    if not objects:
        if note:
            return ["<h2>The projectives and injectives of A</h2>",
                    "<p>(unavailable: %s)</p>" % _esc(note)]
        return []
    out = ["<h2>The projectives and injectives of A</h2>",
           "<p><i>Loewy layers stacked top to bottom; the simples "
           "S<sub>v</sub> are omitted.</i></p>"]
    for row in objects:
        v = row["vertex"]
        for sym in ("P", "I"):
            d = row[sym]
            layers = r" \;\big|\; ".join(factor_stack_tex(L) for L in d["layers"]) \
                or "0"
            out.append(_math(r"%s_{%s} = %s \qquad \dim = %d,\ \underline{\dim} = %s"
                             % (sym, v, layers, d["dim"], _dimvec_tex(d["dimvec"]))))
    return out


def _describe_modules(modules):
    """``[(name, Module), ...]`` -> description dicts, each guarded independently:
    a module that cannot be described is skipped, never fatal (the report is a
    record of a computation that already succeeded)."""
    out = []
    for name, mod in (modules or ()):
        if mod is None:
            continue
        try:
            from quiverlab.trace.modules import module_description
            out.append(module_description(mod, name))
        except Exception:                    # descriptive only -- never sink a report
            continue
    return out


def _modules_section_html(descriptions):
    """"The modules" section: for each module the computation was about, its
    dimension vector, its Loewy series, and the exact per-arrow action matrices
    (Marco 2026-07-29 -- a dimension vector alone does not say what the module IS).
    An arrow acting as the exact zero map is named, not printed."""
    if not descriptions:
        return []
    out = ["<p><i>Each module as it was given to the engine: its Loewy layers "
           "(stacked top to bottom) and the exact matrix of every arrow.</i></p>"]
    for d in descriptions:
        name = d["name"]
        layers = r" \;\big|\; ".join(factor_stack_tex(L) for L in d["layers"]) or "0"
        out.append("<h3>%s</h3>" % _esc(name))
        out.append(_math(r"%s = %s \qquad \dim = %d,\ \underline{\dim} = %s"
                         % (name, layers, d["dim"], _dimvec_tex(d["dimvec"]))))
        out.append(_math(r"\operatorname{top} %s = %s,\qquad \operatorname{soc} %s = %s"
                         % (name, factor_stack_tex(d["top"]),
                            name, factor_stack_tex(d["socle"]))))
        if d.get("side") == "left":
            out.append("<p class='ql-note'>a LEFT module (a right module over "
                       "A<sup>op</sup>).</p>")
        if d.get("display_only"):
            out.append("<p class='ql-note'>display only — entries lie outside the "
                       "integer/fraction input grammar (e.g. GF(p^n) elements).</p>")
        out.extend(_arrow_maps_html(name, d.get("maps") or {}))
    return out


def _arrow_maps_html(label, maps):
    """One matrix per arrow acting NON-trivially; the zero arrows named in one line
    (Marco 2026-07-29). Shared by the modules section and the summand tables."""
    from quiverlab.trace.results_html import _is_zero
    live = sorted(a for a in maps if not _is_zero(maps[a]))
    zero = sorted(a for a in maps if _is_zero(maps[a]))
    if not live:
        return ["<p class='ql-note'>%s: every arrow acts as zero.</p>" % _esc(label)]
    out = []
    for a in live:
        out.append("<p>%s, arrow %s:</p>" % (_esc(label), _esc(str(a))))
        out.append(matrix_grid(maps[a]))
    if zero:
        out.append("<p class='ql-note'>%s: arrow%s %s act%s as zero.</p>"
                   % (_esc(label), "s" if len(zero) > 1 else "",
                      _esc(", ".join(zero)), "" if len(zero) > 1 else "s"))
    return out


def _ext_degree_html(e):
    """One Ext/Tor degree with the full rank arithmetic spelled out, op-aware
    (Plan 34 MAJOR-2). For Ext:
    ``\\delta^n``, ``dim Hom``, ``Ext^n`` (superscript). For Tor: ``d_{n+1}``, tensor
    dimension, ``Tor_n`` (subscript). The dimension count ``space - rank - rank =
    result`` is rendered too, so the no-toolchain HTML surface matches the PDF."""
    if e.op == "Tor":
        # homology: Tor_n = ker d_n / im d_{n+1}; the SHOWN matrix is d_{n+1}
        # (the map INTO degree n, rank_here); the outgoing d_n has rank_prev.
        shown, other = "d_{%d}" % (e.degree + 1), "d_{%d}" % e.degree
        space_word = r"\dim(Q_{%d}\otimes_A N)" % e.degree
        quotient = r"\ker d_{%d}/\operatorname{im}d_{%d}" % (e.degree, e.degree + 1)
        arith = "%d - %d - %d" % (e.space_dim, e.rank_prev, e.rank_here)
        label = r"\operatorname{Tor}_{%d}(M,N)" % e.degree
    else:  # Ext (cohomology): Ext^n = ker delta^n / im delta^{n-1}
        shown, other = r"\delta^{%d}" % e.degree, r"\delta^{%d}" % (e.degree - 1)
        space_word = r"\dim\operatorname{Hom}(Q_{%d},N)" % e.degree
        quotient = r"\ker\delta^{%d}/\operatorname{im}\delta^{%d}" % (e.degree, e.degree - 1)
        arith = "%d - %d - %d" % (e.space_dim, e.rank_here, e.rank_prev)
        label = r"\operatorname{Ext}^{%d}(M,N)" % e.degree
    return [
        _math(r"%s: %s = %d,\quad \operatorname{rank}%s = %d,\quad "
              r"\operatorname{rank}%s = %d"
              % (label, space_word, e.space_dim, shown, e.rank_here, other, e.rank_prev)),
        _event_grid(e, label=shown),
        _math(r"%s = %s,\qquad \dim = %s = %d"
              % (label, quotient, arith, e.result_dim)),
    ]


def _module_steps_html(events):
    # A products chapter owns its StepNotes (the product definition); it is rendered
    # by _products_html, not here, so a product stream never spawns a spurious
    # "Worked module steps" section from its lone definitional StepNote.
    if any(isinstance(e, ProductStep) for e in events):
        return []
    mods = [e for e in events
            if isinstance(e, (ModuleTerm, ModuleDifferential, ExtDegree, ExtReps,
                              StepNote))]
    if not mods:
        return []
    out = ["<h2>Worked module steps</h2>", "<p><i>%s</i></p>" % _esc(ELISION_PREAMBLE)]
    step_no = 0
    echo = _MatrixEcho()
    for e in mods:
        if isinstance(e, StepNote):
            if getattr(e, "heading", False):
                # A numbered worked step (homework style): the HTML mirror of the LaTeX
                # renderer's run-in ``\paragraph{Step N. ...}`` (Plan 34 MAJOR-2).
                step_no += 1
                out.append("<p><b>Step %d. %s</b>%s</p>" % (
                    step_no, _esc(e.text),
                    " " + _esc(e.detail) if e.detail else ""))
            else:
                out.append("<p>%s%s</p>" % (
                    _esc(e.text),
                    "<br><i>%s</i>" % _esc(e.detail) if e.detail else ""))
        elif isinstance(e, ModuleTerm):
            what = "Q_{%d}" % e.degree if e.sym == "P" else "E^{%d}" % e.degree
            out.append(_math(r"%s = %s \qquad \dim = %d"
                             % (what, oplus_tex(e.summands, e.sym), e.dim)))
        elif isinstance(e, ModuleDifferential):
            name = getattr(e, "mod_name", "M") or "M"
            dom = name if getattr(e, "dom_is_module", False) \
                else oplus_tex(e.dom_summands, e.sym)
            cod = name if e.cod_is_module else oplus_tex(e.cod_summands, e.sym)
            prior = echo.label_for(e, e.symbol)
            out.append(_math(r"%s : %s \to %s" % (e.symbol, dom, cod)))
            if prior is not None:                     # identical to an earlier map
                out.append(_math(r"%s = %s" % (e.symbol, prior)))
                out.append("<p class='ql-note'>(the same matrix as above; not "
                           "repeated)</p>")
            else:
                out.append(_event_grid(e, label=e.symbol))
        elif isinstance(e, ExtReps):
            secs = module_reps_sections(e.basis_classes, e.chain_basis,
                                        e.differentials, e.op, anchor_prefix="ws")
            if secs:
                out.append("<h3 id='ws-module-reps'>Explicit representatives by "
                           "degree</h3>")
                out.append("<p class='ql-note'>Each class is written over the ordered "
                           "basis, with the differential that annihilates it; the "
                           "coordinate vectors are recorded in the JSON.</p>")
                out.extend(secs)
        elif isinstance(e, ExtDegree):
            out.extend(_ext_degree_html(e))
    # ALL Ext/Tor runs (not only the first): an Ext run then a Tor run both appear,
    # each with the correct subscript/superscript (Plan 34 MINOR).
    # Named for WHAT it is (Ext / Tor), not "Result" -- the whole page is results.
    runs = ext_result_runs(events)
    for op, dims in runs:
        sep = "_" if op == "Tor" else "^"
        out.append("<h3 id='ws-run-%s'>%s</h3>" % (op.lower(), _esc(op)))
        out.append(_dims_table("dim %s%sn" % (op, sep), dims))
    return out


# --------------------------------------------------------------------------- #
# Plan 35: the HH-product worked-steps chapter. The ProductStep stream (built by
# trace.products) carries the DATA (per-bidegree equation lines / induced-B
# matrices); the fixed per-kind DEFINITIONAL FORMULA is the renderer's, since this
# module is the sole owner of the worked-steps math source (see the module
# docstring). Anything the MathML converter cannot typeset falls back to escaped
# source (readable), so the definitions never garble.
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Plan 35 UNIT 2: the per-degree EXPLICIT-REPRESENTATIVES layout, shared by the
# products worked-steps chapter (render_html, driven by a ProductBasis event) and by
# the Computed-results product block (results_html, driven by the block dict). Both
# call `product_degree_sections` with a surface-distinct anchor prefix, so a report
# carrying both surfaces gets unique, referenceable anchors. Rendering-only: the data
# is UNIT 1's, captured at table-build time; nothing is recomputed here.
#
# Each (side, degree) becomes one sub-section: (a) the ordered basis enumeration of
# the ambient (co)chain space (HTML-visible, capped with a machine-record pointer),
# (b) the explicit classes as term-sum + inline coordinate vector over that
# enumeration, (c) the annihilating differential (indexed grid, or a stated note when
# UNIT 1 elided the body) + a one-line self-cert verification sentence.
# --------------------------------------------------------------------------- #

# Marco 2026-08-02: the report shows at most the first DISPLAY_CAP basis elements of a
# LISTED space, and never a matrix grid with DISPLAY_CAP or more rows/columns -- beyond
# that it states the size and points at the accompanying JSON, which always carries the
# complete data. One constant for every ordinary display cap. TWO deliberate exceptions
# (both Marco 2026-08-02, and only these): (a) the product Cayley tables are UNCAPPED
# (products.combined_cayley / cayley_grid_html always render -- product tables can be
# big); (b) the flat (co)homology class lists (alpha^n_i / z^n_i) are UNCAPPED wherever
# they render (products, Connes, and the HH degree-section class blocks). The
# chain/cochain-space ENUMERATIONS in the HH sections keep this cap.
DISPLAY_CAP = 20

# ordered-enumeration entries shown inline before a machine-record pointer (display
# only; the full ordered enumeration always lives in result.json / trace.json).
_REPS_ENUM_DISPLAY = DISPLAY_CAP

# per side: (long name, HH symbol, differential TeX symbol, chain kind, class-symbol
# letter, cocycle/cycle word).
_REPS_SIDE = {
    "coh": ("cohomology", "HH^", r"\delta", "cochain", r"\alpha", "cocycle"),
    "hom": ("homology", "HH_", "b", "chain", "z", "cycle"),
}


def _signed_join(pieces):
    """``[(is_negative, magnitude_str), ...]`` -> a linear-combination string with the
    correct leading sign and interior ``+``/``-`` separators (``0`` when empty)."""
    out = []
    for i, (neg, mag) in enumerate(pieces):
        if i == 0:
            out.append(("-" + mag) if neg else mag)
        else:
            out.append((" - " if neg else " + ") + mag)
    return "".join(out) if out else "0"


def _coeff_split(coeff):
    """(is_negative, magnitude) of an exact coefficient string. GF(p) coeffs are
    always the non-negative residue, so only QQ/int coeffs are ever negative."""
    c = str(coeff)
    neg = c.startswith("-")
    return neg, (c[1:] if neg else c)


def _term_sum_text(terms, kind):
    """The class's labeled term-sum as readable text, e.g. ``[x → x]`` (cochain) /
    ``e_1 ⊗ x`` (chain); a non-unit coefficient is shown. Reuses UNIT 1's
    ``element_label`` (the single owner of the term-sum labelling) -- never re-derived."""
    from quiverlab.hochschild.basis_reps import element_label
    pieces = []
    for coeff, word, value in terms:
        neg, mag = _coeff_split(coeff)
        lab = element_label(tuple(word), value, kind).replace("->", "→").replace("(x)", "⊗")
        pieces.append((neg, lab if mag == "1" else "%s %s" % (mag, lab)))
    return _signed_join(pieces)


def _enumeration_html(enum, long_name, hh, n):
    """The ordered basis enumeration of the ambient (co)chain space as a numbered list
    (1-based: entry k is the symbol ``e_k`` the coordinate vectors use). Capped for
    display -- an over-long OR record-elided enumeration states its size and points at
    the machine record."""
    ambient = "C^{%d}" % n if hh == "HH^" else "C_{%d}" % n
    intro = ("<p>Ordered basis of the degree-%d %s space %s "
             "(entry <i>k</i> is the <i>k</i>-th basis element; the JSON coordinate "
             "vectors index into it):</p>"
             % (n, long_name, _math_inline(ambient)))
    if isinstance(enum, dict):                    # UNIT-1 record-elided enumeration
        return [intro, "<p class='ql-note'>%s elements; the full ordered enumeration "
                "is in the machine record.</p>" % _esc(str(enum.get("length", "?")))]
    enum = list(enum or [])
    if not enum:
        return [intro, "<p class='ql-note'>the space is zero-dimensional.</p>"]
    items = "".join(
        "<li>%s</li>" % _esc(str(lbl).replace("->", "→").replace("(x)", "⊗"))
        for lbl in enum[:_REPS_ENUM_DISPLAY])
    out = [intro, "<ol class='ql-enum'>%s</ol>" % items]
    if len(enum) > _REPS_ENUM_DISPLAY:
        out.append("<p class='ql-note'>… %d more (the full ordered enumeration is in "
                   "the machine record).</p>" % (len(enum) - _REPS_ENUM_DISPLAY))
    return out


def _classes_html(classes, letter, n, chain_kind):
    """The explicit basis classes of one degree, each as its labeled term-sum AND its
    coordinate vector over the enumeration (UNIT 1's two coherent views). The
    (co)homology CLASS list is UNCAPPED (Marco 2026-08-02: "no limit for the bases of
    (co)homology") -- only the chain-space enumeration above it is capped."""
    if not classes:
        return ["<p class='ql-note'>no classes (the space is zero in this degree).</p>"]
    lis = []
    for i, cl in enumerate(classes, start=1):
        name = "%s^{%d}_{%d}" % (letter, n, i)
        term = _term_sum_text(cl.get("terms") or [], cl.get("kind") or chain_kind)
        lis.append("<li>%s = %s</li>" % (_math_inline(name), _esc(term)))
    return ["<p>Basis classes, each written over the ordered basis above:</p>",
            "<ul class='ql-classes'>%s</ul>" % "".join(lis)]


def _reps_differential_html(diff, dsym, n, letter, cyc, hh, has_classes):
    """The degree's annihilating differential as an indexed grid (or a stated note when
    UNIT 1 elided the body), plus the one-line verification sentence -- the self-cert
    the reader can carry out from the shipped vector + this matrix."""
    if diff is None:
        return []
    if hh == "HH^":
        symbol, arrow = ("%s^{%d}" % (dsym, n),
                         r"%s^{%d} : C^{%d} \to C^{%d}" % (dsym, n, n, n + 1))
    else:
        symbol, arrow = ("%s_{%d}" % (dsym, n),
                         r"%s_{%d} : C_{%d} \to C_{%d}" % (dsym, n, n, max(n - 1, 0)))
    out = [_math(arrow)]
    if diff.get("elided"):
        r, c = (diff.get("shape") or [0, 0])[:2]
        out.append("<p class='ql-note'>%s×%s matrix (body in the machine record; "
                   "rebuild: %s).</p>" % (_esc(str(r)), _esc(str(c)),
                                          _esc(str(diff.get("note", "")))))
    else:
        out.append(matrix_grid(diff.get("rows") or [], label=symbol))
    if not has_classes:
        return out
    if hh == "HH_" and n == 0:
        sentence = ("every 0-chain is a %s (%s vanishes), so each %s is a %s"
                    % (cyc, symbol, "z^{%d}_i" % n, cyc))
    else:
        sentence = ("each %s^{%d}_i is a %s: applying %s to its coordinate vector "
                    "gives 0" % (letter, n, cyc, symbol))
    out.append("<p class='ql-note'>Verification: %s.</p>" % _esc(sentence))
    return out


def _bar_term_basis_chunks(m, labels, n, collapsed_dim, kind):
    """The ordered basis of a bar HH (co)chain term at degree ``n``, reconstructed with
    UNIT 1's enumeration builders (Plan 35 UNIT 2, deliverable 2 -- HH worked-steps
    resolution narration). Rendering-only: the labels are NOT recomputed ad hoc -- they
    come from ``basis_reps.bar_chain_elements`` (the same builder the reps capture uses),
    and are shown ONLY when the reconstructed enumeration length equals the recorded
    term dimension, so the order provably matches the differential the step shows.

    Returns ``[]`` when the term is NOT the reconstructable bar (co)chain basis (e.g. a
    Chouhy-Solotar term, which carries ``corners`` instead and is named as a bimodule
    direct sum), so the section never mislabels."""
    from quiverlab.hochschild import basis_reps as BR
    expected = m * (m - 1) ** n
    if collapsed_dim != expected:               # not the bar (co)chain basis -- omit
        return []
    if expected <= _REPS_ENUM_DISPLAY:
        enum = BR.enumeration_labels(BR.bar_chain_elements(m, n, labels), kind)
    else:                                        # display-elided: point at the record
        enum = {"length": expected}
    long_name = "cochain" if kind == "cochain" else "chain"
    hh = "HH^" if kind == "cochain" else "HH_"
    return _enumeration_html(enum, long_name, hh, n)


def _cs_term_basis_chunks(cs_res, cs_labels, n, collapsed_dim, side_key, kind):
    """The ordered basis of a Chouhy-Solotar HH (co)chain term at degree ``n``, from the
    rebuilt CS resolution's ``res._basis(n, side)`` labels (Plan 35 UNIT 2 review fix --
    the CS resolution worked steps used to show the raw differential grid with no way to
    read its columns). Reuses UNIT 1's ``basis_reps.cs_elements`` (the SAME builder the CS
    product capture uses -- not ad hoc). The CS resolution is deterministic/canonical
    (Plan 17), so the rebuilt basis order matches the recorded run; shown ONLY when the
    reconstructed length equals the recorded term dimension, so the columns/rows the shown
    differential uses are the ones enumerated. Returns ``[]`` on any rebuild surprise."""
    from quiverlab.hochschild import basis_reps as BR
    try:
        elems = BR.cs_elements(cs_res, n, side_key, cs_labels)
    except Exception:
        return []
    if len(elems) != collapsed_dim:              # rebuilt basis != recorded term -- omit
        return []
    enum = BR.enumeration_labels(elems, kind)
    long_name = "cochain" if kind == "cochain" else "chain"
    hh = "HH^" if kind == "cochain" else "HH_"
    return _enumeration_html(enum, long_name, hh, n)


def product_degree_sections(basis_classes, chain_basis, differentials, anchor_prefix,
                            show_differential=True):
    """Per-degree explicit-representatives sub-sections for a product/Connes block
    (Plan 35 UNIT 2). ``basis_classes`` / ``chain_basis`` / ``differentials`` are the
    UNIT-1 ``{side: {str(degree): ...}}`` payloads; ``anchor_prefix`` namespaces the
    ``<prefix>-hh-<side>-deg-<n>`` anchors the product tables reference.

    ``show_differential`` (Marco 2026-07-31): the plain HH degree sections keep the
    annihilating differential; the PRODUCT sections drop it (they state only the
    products in terms of the basis classes -- the differentials already live in the HH
    degree sections). A degree whose (co)homology is ZERO collapses to a one-line
    ``HH^n = 0`` statement, keeping its anchor so the vanishing is still findable.

    Returns ``[]`` when the block carries no explicit-reps fields (a legacy/old-cache
    block) -- the caller then falls back to the naming-only legend (tolerance)."""
    if not basis_classes:
        return []
    out = []
    for side in ("coh", "hom"):
        by_deg = basis_classes.get(side)
        if not by_deg:
            continue
        long_name, hh, dsym, chain_kind, letter, cyc = _REPS_SIDE[side]
        for dkey in sorted(by_deg, key=lambda s: int(s)):
            n = int(dkey)
            anchor = "%s-hh-%s-deg-%d" % (anchor_prefix, side, n)
            classes = by_deg[dkey] or []
            out.append("<h4 id='%s'>Hochschild %s in degree %d</h4>"
                       % (anchor, long_name, n))
            if not classes:                     # zero space: one line, keep the anchor
                out.append(_math(r"%s{%d} = 0" % (hh, n)))
                continue
            enum = (chain_basis or {}).get(side, {}).get(dkey)
            out.extend(_enumeration_html(enum, long_name, hh, n))
            out.extend(_classes_html(classes, letter, n, chain_kind))
            if show_differential:
                diff = (differentials or {}).get(side, {}).get(dkey)
                out.extend(_reps_differential_html(diff, dsym, n, letter, cyc, hh,
                                                   bool(classes)))
    return out


def product_flat_classes_html(basis_classes):
    """ONE compact flat list of ALL (co)homology basis classes across degrees (Marco
    2026-08-02): cohomology classes ``α^n_i`` first, then -- for the cap -- the homology
    classes ``z^n_i``, ordered degree-major (the degree is the class' superscript). Each
    class is written as its representative term-sum over the chain basis ALREADY
    enumerated in the HH (co)homology sections above -- the enumeration, differential and
    per-degree headings are NOT repeated here (that is Marco's whole point: for the
    products just remind the (co)homology basis, all at once, then show the table). The
    flat (co)homology class list is UNCAPPED (Marco 2026-08-02: "no limit for the bases
    of (co)homology").

    Returns ``[]`` when the block carries no ``basis_classes`` (a legacy/old-cache block)
    -- the caller then shows the tables only (tolerance). Shared by the products
    worked-steps chapter (``render_html._products_html``) and the Computed-results product
    block (``results_html._product_tables_html``) -- one implementation, no drift."""
    if not basis_classes:
        return []
    items = []
    for side in ("coh", "hom"):
        by_deg = basis_classes.get(side)
        if not by_deg:
            continue
        chain_kind, letter = _REPS_SIDE[side][3], _REPS_SIDE[side][4]
        for dkey in sorted(by_deg, key=lambda s: int(s)):
            n = int(dkey)
            for i, cl in enumerate(by_deg[dkey] or [], start=1):
                name = "%s^{%d}_{%d}" % (letter, n, i)
                term = _term_sum_text(cl.get("terms") or [], cl.get("kind") or chain_kind)
                items.append("<li>%s = %s</li>" % (_math_inline(name), _esc(term)))
    if not items:
        return []
    return ["<p>The Hochschild (co)homology basis classes, over the chain bases "
            "enumerated in the sections above:</p>",
            "<ul class='ql-classes'>%s</ul>" % "".join(items)]


# --------------------------------------------------------------------------- #
# Plan 35 wave 3d: the PLAIN hh_cohomology / hh_homology blocks. The block carries a
# SINGLE-side `{str(degree): ...}` reps payload (like Ext / Tor / HC); the per-degree
# sections reuse `product_degree_sections` (the SAME HH cochain/chain renderer) by
# wrapping the single side, and the element-wise dictionary read-offs (central elements
# / derivations / deformation cochain / commutator residues) are read straight off the
# captured term-sums via `quiverlab.trace.interpretations`. Element-wise ONLY where reps
# are present; the framing sentence (results_html._dictionary_framing_html) covers the
# rest. Data is the capture layer's (hochschild.hh_reps); nothing is recomputed here.
# --------------------------------------------------------------------------- #
def hh_element_interpretation(kind, basis_classes, inner_dims):
    """The element-wise dictionary read-offs of the plain HH block, per degree: HH^0's
    central elements, HH^1's derivations ``D(arrow) = value`` (+ the inner-derivation
    subspace dimension ``rank δ^0``, shipped as ``inner_dims["1"]``), HH^2's deformation
    2-cocycle, HH_0's commutator residues. ``[]`` when the block carries no reps."""
    from quiverlab.trace.interpretations import element_heading, element_readoff
    if not basis_classes:
        return []
    letter = r"\alpha" if kind == "hh_cohomology" else "z"
    out = []
    for dkey in sorted(basis_classes, key=lambda s: int(s)):
        n = int(dkey)
        heading = element_heading(kind, n)
        classes = basis_classes.get(dkey) or []
        if heading is None or not classes:
            continue
        lis = []
        for i, cl in enumerate(classes, start=1):
            name = "%s^{%d}_{%d}" % (letter, n, i)
            lines = element_readoff(kind, n, cl.get("terms") or []) or []
            body = "; ".join(str(s).replace("->", "→") for s in lines) or "0"
            lis.append("<li>%s: %s</li>" % (_math_inline(name), _esc(body)))
        out.append("<p><i>%s</i></p>" % _esc(heading))
        out.append("<ul class='ql-interp'>%s</ul>" % "".join(lis))
        if kind == "hh_cohomology" and n == 1 and inner_dims is not None:
            out.append("<p class='ql-note'>inner derivations (the coboundaries "
                       "a ↦ ax − xa): dimension %s = rank δ⁰.</p>"
                       % _esc(str(inner_dims.get("1", "?"))))
    return out


def hh_reps_sections(kind, basis_classes, chain_basis, differentials, anchor_prefix):
    """Per-degree explicit-representatives sub-sections for a plain HH block, reusing the
    products renderer `product_degree_sections` by wrapping the single side (coh for
    hh_cohomology, hom for hh_homology). ``[]`` when the block carries no reps
    (legacy/old-cache tolerance)."""
    if not basis_classes:
        return []
    side = "coh" if kind == "hh_cohomology" else "hom"
    return product_degree_sections({side: basis_classes},
                                   {side: chain_basis} if chain_basis else None,
                                   {side: differentials} if differentials else None,
                                   anchor_prefix)


# --------------------------------------------------------------------------- #
# Plan 35 wave 3a: the per-degree EXPLICIT-REPRESENTATIVES layout for module Ext /
# Tor, the sibling of `product_degree_sections`. Ext / Tor blocks carry a SINGLE-side
# `{str(degree): ...}` payload (Ext is cohomological, Tor homological), so this reader
# is kind-scoped and never confused with the products `{side: {degree}}` shape or the
# module-resolution LIST-shaped `differentials`. Shared by the Computed-results block
# (results_html) and the module worked-steps chapter (an ExtReps event). Data is the
# capture layer's (modules.complex_reps); nothing is recomputed here.
# --------------------------------------------------------------------------- #

# per kind: (ambient long-name template, class-symbol letter, (co)cycle word,
# differential TeX symbol template, differential arrow template, verification tail).
_MODULE_REPS = {
    "ext": ("\\mathrm{Hom}_A(P_{%d}, N)", r"\alpha", "cocycle",
            r"\delta^{%d}", r"\delta^{%d} : \mathrm{Hom}(P_{%d},N) \to "
            r"\mathrm{Hom}(P_{%d},N)",
            "applying %s to its coordinate vector gives 0"),
    "tor": ("P_{%d} \\otimes_A N", "z", "cycle",
            r"d_{%d}", r"d_{%d} : P_{%d}\otimes_A N \to P_{%d}\otimes_A N",
            "applying %s to its coordinate vector gives 0"),
}


def _module_term_sum_text(terms, kind):
    """The class's labelled term-sum: Ext ``[g -> v]`` (the explicit hom sending the
    P_n generator ``g`` to the N-basis vector ``v``), Tor ``g (x) v`` (a tensor). A
    non-unit coefficient is shown. Owner of the module element label (single source)."""
    from quiverlab.modules.complex_reps import element_label
    pieces = []
    for coeff, gen, val in terms:
        neg, mag = _coeff_split(coeff)
        lab = element_label(gen, val, kind).replace("->", "→").replace("(x)", "⊗")
        pieces.append((neg, lab if mag == "1" else "%s %s" % (mag, lab)))
    return _signed_join(pieces)


def _module_classes_html(classes, letter, n, kind):
    """The explicit basis classes of one degree, each as its labelled term-sum (a
    self-contained representative -- e.g. Ext ``[g -> v]`` / Tor ``g (x) v``). The
    ordered Hom/tensor basis it indexes into is NOT re-listed here (Marco 2026-08-02:
    the complex-space enumerations live only in the JSON record); the class list keeps
    the ordinary :data:`DISPLAY_CAP` listing cap with a JSON pointer."""
    if not classes:
        return ["<p class='ql-note'>no classes (the group is zero in this degree).</p>"]
    lis = []
    for i, cl in enumerate(classes[:DISPLAY_CAP], start=1):
        name = "%s^{%d}_{%d}" % (letter, n, i)
        term = _module_term_sum_text(cl.get("terms") or [], kind)
        lis.append("<li>%s = %s</li>" % (_math_inline(name), _esc(term)))
    out = ["<p>Basis classes, each as its labelled representative:</p>",
           "<ul class='ql-classes'>%s</ul>" % "".join(lis)]
    if len(classes) > DISPLAY_CAP:
        out.append("<p class='ql-note'>… and %d more classes (see the JSON record).</p>"
                   % (len(classes) - DISPLAY_CAP))
    return out


def _module_reps_differential_html(diff, kind, n, letter, cyc):
    """The degree's annihilating differential as an indexed grid (or a stated note when
    the capture layer elided the body), plus the one-line verification sentence."""
    if diff is None:
        return []
    _amb, _l, _cyc, dsym_t, arrow_t, verify_t = _MODULE_REPS[kind]
    if kind == "ext":
        symbol, arrow = dsym_t % n, arrow_t % (n, n, n + 1)
    else:
        symbol, arrow = dsym_t % n, arrow_t % (n, n, max(n - 1, 0))
    out = [_math(arrow)]
    if diff.get("elided"):
        r, c = (diff.get("shape") or [0, 0])[:2]
        out.append("<p class='ql-note'>%s×%s matrix (body in the machine record; "
                   "rebuild: %s).</p>" % (_esc(str(r)), _esc(str(c)),
                                          _esc(str(diff.get("note", "")))))
    elif diff.get("shape") and diff["shape"][0] == 0:
        # a zero-row map (Tor d_0): every 0-chain is a cycle; state the note verbatim.
        note = diff.get("note")
        if note:
            out.append("<p class='ql-note'>%s.</p>" % _esc(str(note)))
        return out
    else:
        out.append(matrix_grid(diff.get("rows") or [], label=symbol))
    out.append("<p class='ql-note'>Verification: each %s^{%d}_i is a %s: %s.</p>"
               % (letter, n, cyc, verify_t % symbol))
    return out


def module_reps_sections(basis_classes, chain_basis, differentials, kind, anchor_prefix):
    """Per-degree explicit-representatives sub-sections for a module Ext / Tor block
    (Plan 35 wave 3a). ``basis_classes`` / ``chain_basis`` / ``differentials`` are the
    single-side ``{str(degree): ...}`` payloads; ``kind`` is ``"ext"`` / ``"tor"``;
    ``anchor_prefix`` namespaces the ``<prefix>-<kind>-deg-<n>`` anchors.

    Returns ``[]`` when the block carries no explicit-reps fields (a legacy/old-cache
    block) -- the caller then falls back to the dims table only (tolerance)."""
    if not basis_classes:
        return []
    from quiverlab.trace.interpretations import module_reps_label_note
    _ambient_t, letter, cyc, _dsym, _arrow, _v = _MODULE_REPS[kind]
    long_name = "Ext^{%d}" if kind == "ext" else "Tor_{%d}"
    # Marco 2026-08-02: the ordered Hom/tensor basis is NOT enumerated per degree here --
    # one pointer states it lives in the JSON; the class list + differentials stay.
    out = ["<p class='ql-note'>%s</p>" % _esc(module_reps_label_note(kind)),
           "<p class='ql-note'>The ordered basis of each Hom/tensor space (into which "
           "the coordinate vectors index) is recorded in the JSON.</p>"]
    for dkey in sorted(basis_classes, key=lambda s: int(s)):
        n = int(dkey)
        anchor = "%s-%s-deg-%d" % (anchor_prefix, kind, n)
        classes = basis_classes.get(dkey) or []
        out.append("<h4 id='%s'>%s in degree %d</h4>"
                   % (anchor, _math_inline(long_name % n), n))
        if not classes:                         # zero group: one line, keep the anchor
            out.append(_math(long_name % n + " = 0"))
            continue
        diff = (differentials or {}).get(dkey)
        out.extend(_module_classes_html(classes, letter, n, kind))
        out.extend(_module_reps_differential_html(diff, kind, n, letter, cyc))
    return out


# --------------------------------------------------------------------------- #
# Plan 35 wave 3c: the Yoneda exact-sequence INTERPRETATION of Ext classes. Each class
# of Ext^n(M, N) is shown as the constructed exact sequence 0 -> N -> Q -> ... -> M -> 0
# with the middle module Q as a full representation and its exactness self-certified.
# Data is the capture layer's (modules.complex_reps `interpretation`); nothing recomputes.
# --------------------------------------------------------------------------- #
def _yoneda_sequence_math(mods):
    """The exact sequence as one typeset line ``0 \\to N \\to E \\to M \\to 0`` using each
    module's display label. The labels are already display-ready (single symbols like
    ``N`` / ``E`` / ``M`` or TeX terms like ``P_{1} \\oplus P_{2}``)."""
    inner = " \\to ".join(m.get("label", "?") for m in mods)
    return "0 \\to %s \\to 0" % inner


def _dv_inline(dv):
    """A dimension vector dict {v: n} as ``(n_1, n_2, ...)`` in sorted vertex order."""
    if not dv:
        return "()"
    return "(" + ", ".join(str(dv[k]) for k in sorted(dv)) + ")"


def _yoneda_facts_sentence(facts):
    parts = []
    for f in facts or []:
        if f.get("fact") == "injective":
            parts.append("injective at %s (rank %s = dim %s)"
                         % (f["node"], f["rank"], f["dim"]))
        elif f.get("fact") == "surjective":
            parts.append("surjective at %s (rank %s = dim %s)"
                         % (f["node"], f["rank"], f["dim"]))
        elif f.get("fact") == "im=ker":
            parts.append("image = kernel at %s (%s + %s = %s)"
                         % (f["node"], f["rank_in"], f["rank_out"], f["dim"]))
    return "; ".join(parts)


def _yoneda_middle_html(mid):
    """The extension / middle module: named when it is a standard indecomposable, else
    its full per-arrow action (the no-code representation)."""
    label = mid.get("label", "Q")
    std = mid.get("standard")
    dv = _dv_inline(mid.get("dimvec") or {})
    if std is not None:
        sym = {"simple": "S", "projective": "P", "injective": "I"}.get(std["kind"], "?")
        return ["<p>%s is the standard indecomposable %s (dimension vector %s).</p>"
                % (_math_inline(label), _math_inline("%s_{%s}" % (sym, std["vertex"])),
                   _esc(dv))]
    out = ["<p>%s, the extension module, dimension vector %s -- its action:</p>"
           % (_math_inline(label), _esc(dv))]
    maps = (mid.get("module") or {}).get("maps") or {}
    out.extend(_arrow_maps_html(label, maps))
    if (mid.get("module") or {}).get("display_only"):
        out.append("<p class='ql-note'>display only — entries lie outside the "
                   "integer/fraction grammar (e.g. GF(p^n) elements).</p>")
    return out


def _yoneda_one_class_html(seq):
    name = seq.get("class_name", "?")
    if not seq.get("certified", False):
        return ["<h5>%s</h5>" % _math_inline(name),
                "<p class='ql-note'>this class's exact sequence could not be certified "
                "(%s); it is omitted rather than shown wrongly.</p>"
                % _esc(str(seq.get("error", "unknown")))]
    mods = seq.get("modules") or []
    out = ["<h5>%s</h5>" % _math_inline(name),
           _math(_yoneda_sequence_math(mods)),
           "<p class='ql-note'>dimension vectors: %s.</p>"
           % _esc(", ".join("%s %s" % (m.get("label", "?"), _dv_inline(m.get("dimvec") or {}))
                            for m in mods))]
    mid = next((m for m in mods if m.get("role") == "middle"), None)
    if mid is not None:
        out.extend(_yoneda_middle_html(mid))
    # connecting maps
    for mp in seq.get("maps") or []:
        arrow = "%s \\to %s" % (mp.get("from", "?"), mp.get("to", "?"))
        out.append(_math(arrow))
        if mp.get("elided"):
            r, c = (mp.get("shape") or [0, 0])[:2]
            out.append("<p class='ql-note'>%s×%s matrix (body in the machine record).</p>"
                       % (_esc(str(r)), _esc(str(c))))
        else:
            out.append(matrix_grid(mp.get("rows") or []))
    facts = _yoneda_facts_sentence(seq.get("facts"))
    if facts:
        out.append("<p class='ql-note'>Exactness verified: %s.</p>" % _esc(facts))
    return out


def ext_interpretation_sections(interpretation, anchor_prefix):
    """The Yoneda exact-sequence interpretation of every Ext class, per degree, with the
    dictionary framing sentence and the constructed + certified sequences.

    Returns ``[]`` for a legacy/old-cache block with no ``interpretation`` field
    (tolerance)."""
    from quiverlab.trace.interpretations import ext_degree
    if not interpretation:
        return []
    sequences = interpretation.get("sequences") or {}
    out = []
    for dkey in sorted(sequences, key=lambda s: int(s)):
        n = int(dkey)
        seqs = sequences.get(dkey) or []
        if not seqs:
            continue
        anchor = "%s-ext-yoneda-deg-%d" % (anchor_prefix, n)
        out.append("<h4 id='%s'>Interpretation: %s as exact sequences</h4>"
                   % (anchor, _math_inline("\\mathrm{Ext}^{%d}(M, N)" % n)))
        out.append("<p><i>%s</i></p>" % _esc(ext_degree(n)))
        for seq in seqs:
            out.extend(_yoneda_one_class_html(seq))
    return out


# --------------------------------------------------------------------------- #
# Plan 35 wave 3b: the per-degree EXPLICIT-REPRESENTATIVES layout for cyclic homology
# HC, the sibling of `module_reps_sections`. The HC block carries a SINGLE-side
# `{str(degree): ...}` payload (HC is homological) PLUS a `column_structure`. Marco
# 2026-08-02: the Tot column enumerations (the total-complex decomposition heading and
# the ordered Tot_n basis) are NOT re-listed per degree -- they live in the JSON; the
# section keeps the class list + total differential. Data is the capture layer's
# (hochschild.cyclic_reps); nothing is recomputed here.
# --------------------------------------------------------------------------- #
def _cyclic_term_sum_from_vector(vector, enum):
    """The HC class's term-sum built from its sparse vector + the column-annotated
    enumeration labels ('col C_d: v (x) w...'), so each term names its column. Returns
    None when the enumeration is display-elided (coordinate vector shown alone)."""
    if not isinstance(enum, list) or not enum:
        return None
    pieces = []
    for idx, coeff in vector:
        neg, mag = _coeff_split(coeff)
        lab = str(enum[idx]).replace("->", "→").replace("(x)", "⊗")
        pieces.append((neg, lab if mag == "1" else "%s %s" % (mag, lab)))
    return _signed_join(pieces)


def _cyclic_classes_html(classes, n, enum):
    if not classes:
        return ["<p class='ql-note'>no classes (HC_%d is zero).</p>" % n]
    lis = []
    for i, cl in enumerate(classes[:DISPLAY_CAP], start=1):
        name = "z^{%d}_{%d}" % (n, i)
        term = _cyclic_term_sum_from_vector(cl.get("vector") or [], enum)
        rhs = _esc(term) if term is not None else \
            "(class recorded in the JSON; the ordered basis is too large to display)"
        lis.append("<li>%s = %s</li>" % (_math_inline(name), rhs))
    out = ["<p>Basis classes, each as its labelled representative:</p>",
           "<ul class='ql-classes'>%s</ul>" % "".join(lis)]
    if len(classes) > DISPLAY_CAP:
        out.append("<p class='ql-note'>… and %d more classes (see the JSON record).</p>"
                   % (len(classes) - DISPLAY_CAP))
    return out


def _cyclic_reps_differential_html(diff, n, has_classes):
    if diff is None:
        return []
    arrow = r"D_{%d} : \mathrm{Tot}_{%d} \to \mathrm{Tot}_{%d}" % (n, n, max(n - 1, 0))
    out = [_math(arrow)]
    if diff.get("elided"):
        r, c = (diff.get("shape") or [0, 0])[:2]
        out.append("<p class='ql-note'>%s×%s matrix (body in the machine record; "
                   "rebuild: %s).</p>" % (_esc(str(r)), _esc(str(c)),
                                          _esc(str(diff.get("note", "")))))
    elif diff.get("shape") and diff["shape"][0] == 0:
        note = diff.get("note")                   # D_0: every 0-chain is a cycle
        if note:
            out.append("<p class='ql-note'>%s.</p>" % _esc(str(note)))
        return out
    else:
        out.append(matrix_grid(diff.get("rows") or [], label="D_{%d}" % n))
    if not has_classes:
        return out
    out.append("<p class='ql-note'>Verification: each z^{%d}_i is a cycle of the total "
               "complex: applying D_{%d} = b + B to its coordinate vector gives 0.</p>"
               % (n, n))
    return out


def cyclic_degree_sections(basis_classes, chain_basis, differentials, column_structure,
                           anchor_prefix):
    """Per-degree explicit-representatives sub-sections for a cyclic-homology block
    (Plan 35 wave 3b). Single-side ``{str(degree): ...}`` payloads (HC is homological)
    plus ``column_structure``. Marco 2026-08-02: the Tot column enumerations (the
    ``Tot_n = C_n (+) C_{n-2} (+) ...`` decomposition heading and the ordered Tot_n basis)
    are NOT re-listed -- one pointer states they live in the JSON; each degree keeps its
    class list -> total differential ``D_n = b + B`` + verification, under a stable anchor
    ``<prefix>-hc-deg-<n>``. Returns ``[]`` on a legacy/old-cache block (tolerance)."""
    if not basis_classes:
        return []
    # One pointer for the whole section (the Tot column structure + ordered basis are in
    # the JSON; the coordinate vectors index into it).
    out = ["<p class='ql-note'>The total complex Tot_n = C_n ⊕ C_{n-2} ⊕ … and the "
           "ordered basis of each Tot_n (into which the coordinate vectors index) are "
           "recorded in the JSON.</p>"]
    for dkey in sorted(basis_classes, key=lambda s: int(s)):
        n = int(dkey)
        anchor = "%s-hc-deg-%d" % (anchor_prefix, n)
        classes = basis_classes.get(dkey) or []
        enum = (chain_basis or {}).get(dkey)
        diff = (differentials or {}).get(dkey)
        out.append("<h4 id='%s'>%s in degree %d</h4>"
                   % (anchor, _math_inline("HC_{%d}" % n), n))
        if not classes:                         # zero space: one line, keep the anchor
            out.append(_math("HC_{%d} = 0" % n))
            continue
        out.extend(_cyclic_classes_html(classes, n, enum))
        out.extend(_cyclic_reps_differential_html(diff, n, bool(classes)))
    return out


_PRODUCT_TITLE = {
    "cup": "The cup product on Hochschild cohomology",
    "cap": "The cap product",
    "bracket": "The Gerstenhaber bracket",
    "connes_b": "The Connes differential B",
}
_PRODUCT_DEF = {
    "cup": (r"(f \cup g)(a_1 \otimes \cdots \otimes a_{p+q}) = "
            r"f(a_1 \otimes \cdots \otimes a_p) \cdot "
            r"g(a_{p+1} \otimes \cdots \otimes a_{p+q})"),
    "cap": (r"(f \cap z)(a_1 \otimes \cdots \otimes a_{n-p}) = "
            r"a_0\, f(a_1 \otimes \cdots \otimes a_p) \otimes "
            r"a_{p+1} \otimes \cdots \otimes a_n"),
    "bracket": (r"[f, g] = f \circ g - (-1)^{(p-1)(q-1)}\, g \circ f, \quad "
                r"HH^{p} \otimes HH^{q} \to HH^{p+q-1}"),
    "connes_b": (r"B : HH_n \to HH_{n+1}, \quad "
                 r"B([a_1 \otimes \cdots \otimes a_n]) = "
                 r"\sum_i (-1)^{ni}\, [1 \otimes a_i \otimes \cdots \otimes a_{i-1}]"),
}


def _products_html(events):
    """The HH-product chapter body: the prose definition (the chapter's StepNote), the
    typeset definitional formula for the kind, the per-degree EXPLICIT-REPRESENTATIVES
    sub-sections (ordered basis -> classes -> differential + verification, Plan 35
    UNIT 2), then one block per bidegree -- the structure-constant equation lines
    (cup/cap/bracket) or the induced Connes differential grid (connes_b), each headed
    by its typeset map label and referencing the degree classes above."""
    steps = [e for e in events if isinstance(e, ProductStep)]
    if not steps:
        return []
    kind = steps[0].kind
    from quiverlab.trace.products import balanced_rep_note
    out = []
    for e in events:                                  # the definitional StepNote(s)
        if isinstance(e, StepNote):
            out.append("<p><b>%s</b>%s</p>" % (
                _esc(e.text),
                "<br><i>%s</i>" % _esc(e.detail) if e.detail else ""))
    formula = _PRODUCT_DEF.get(kind)
    if formula:
        out.append(_math(formula))
    # The balanced-representative legend, once per family, when the constants are GF(p)
    # residues (Marco 2026-08-01). The tables show c-p for c > p/2; the JSON keeps raw.
    prime = next((s.prime for s in steps
                  if s.kind in ("cup", "cap", "bracket")
                  and getattr(s, "prime", None) is not None), None)
    if prime is not None:
        out.append("<p class='ql-note'>%s</p>" % _esc(balanced_rep_note(prime)))
    # The (co)homology basis classes as ONE flat list of ALL classes across degrees
    # (Marco 2026-08-02), then the tables/matrices right away -- no per-degree
    # sub-sections and no chain enumerations (those live ONLY in the HH cohomology/
    # homology sections). connes_b uses the SAME flat homology (z^n_j) class list, then
    # its induced B matrices per degree with rank lines. Driven by the ProductBasis event
    # the chapter emits; absent on a legacy object -> the section falls back to tables
    # only.
    pb = next((e for e in events if isinstance(e, ProductBasis)), None)
    if pb is not None:
        out.extend(product_flat_classes_html(pb.basis_classes))
    if kind == "connes_b":
        out.append("<p class='ql-note'>Each cycle class is written over the chain basis "
                   "enumerated in the Hochschild homology sections / JSON; the induced "
                   "matrices below act on these classes (rows index HH_{n+1}, columns "
                   "HH_n).</p>")
        for s in steps:
            if s.heading:
                out.append('<p class="ql-mlabel">%s</p>' % _math_inline(s.heading))
            if s.matrix is not None:
                out.append(matrix_grid(s.matrix))
            if s.note:
                out.append("<p class='ql-note'>%s</p>" % _esc(s.note))
        return out
    # cup / cap / bracket: ONE big degree-major Cayley table for the whole family
    # (Marco 2026-08-01 addendum). An all-vanishing family collapses to one line
    # (handled inside family_cayley_html so both report surfaces share the rule).
    table_steps = [s for s in steps if s.kind in ("cup", "cap", "bracket")]
    tables_data = [{"degrees": tuple(s.degrees), "out_degree": s.out_degree,
                    "dims": list(s.dims), "constants": s.constants}
                   for s in table_steps if s.constants is not None]
    if not tables_data:                               # legacy events: equation lines
        for s in table_steps:
            if s.heading:
                out.append('<p class="ql-mlabel">%s</p>' % _math_inline(s.heading))
            for line in (s.lines or ()):
                out.append(_math(line))
            if s.note:
                out.append("<p class='ql-note'>%s</p>" % _esc(s.note))
        return out
    prime = next((s.prime for s in table_steps
                  if getattr(s, "prime", None) is not None), None)
    out.extend(family_cayley_html(kind, tables_data, prime))
    return out


_FAMILY_HEADING = {
    "cup": r"HH^{*} \cup HH^{*} \to HH^{*}",
    "cap": r"HH^{*} \cap HH_{*} \to HH_{*}",
    "bracket": r"[HH^{*}, HH^{*}] \to HH^{*}",
}
_FAMILY_AXIS_NOTE = (
    "One row per left class and one column per right class, ordered degree-major "
    "(the degree is the class' superscript); a heavier rule marks each degree "
    "boundary.")


def family_cayley_html(kind, tables_data, prime):
    """Render a whole cup/cap/bracket family (given its per-bidegree table dicts) as ONE
    big degree-major Cayley table -- the family heading, the axis/degree-boundary note,
    the whole-region structural caption, the em-dash beyond-window legend, then the
    grid. Marco 2026-08-02: the table is UNCAPPED (product tables can be big) -- there is
    no per-bidegree fallback; the one combined table always renders. Shared by the
    worked-steps chapter and the Computed-results section (no drift). ``tables_data``
    empty -> ``[]``; an all-vanishing family keeps the one-line statement, never a grid
    of 0s (the beyond-window/zero semantics are unchanged)."""
    from quiverlab.trace.products import (
        beyond_window_note, combined_cayley, equation_lines)
    tables_data = [t for t in tables_data if t.get("constants") is not None]
    if not tables_data:
        return []
    # An all-vanishing family keeps the one-line statement (Marco), never a grid of 0s.
    if all(not equation_lines(kind, list(t["degrees"]), t["out_degree"],
                              list(t["dims"]), t["constants"]) for t in tables_data):
        name = {"cup": "cup products", "cap": "cap products",
                "bracket": "Gerstenhaber brackets"}.get(kind, "products")
        return ["<p class='ql-note'>All %s in the served bidegrees vanish.</p>" % name]
    out = ['<p class="ql-mlabel">%s</p>' % _math_inline(_FAMILY_HEADING[kind]),
           "<p class='ql-note'>%s</p>" % _esc(_FAMILY_AXIS_NOTE)]
    combined = combined_cayley(kind, tables_data, prime)
    if combined["note"]:
        out.append("<p class='ql-note'>%s</p>" % _esc(combined["note"]))
    if combined["has_beyond"]:
        out.append("<p class='ql-note'>%s</p>" % _esc(beyond_window_note()))
    out.append(cayley_grid_html(combined))
    return out


# Integer 24-point unit circle, scaled by 100 (index k = k * 15 degrees).
# Display-only layout data: the no-floats gate holds -- every coordinate below
# and every derived position is an int.
_CIRCLE24 = [(100, 0), (97, 26), (87, 50), (71, 71), (50, 87), (26, 97),
             (0, 100), (-26, 97), (-50, 87), (-71, 71), (-87, 50), (-97, 26),
             (-100, 0), (-97, -26), (-87, -50), (-71, -71), (-50, -87),
             (-26, -97), (0, -100), (26, -97), (50, -87), (71, -71),
             (87, -50), (97, -26)]


def _quiver_svg(algebra):
    """Inline SVG sketch of the presentation quiver (vertices on a circle,
    labeled arrows, loops), so the reader sees WHAT the example is before any
    mathematics. Integer arithmetic only. Returns '' without a presentation."""
    import math as _math_mod

    q = getattr(algebra, "quiver", None)
    if q is None or not getattr(q, "vertices", None):
        return ""
    verts = list(q.vertices)
    arrows = dict(getattr(q, "arrows", {}) or {})
    n = len(verts)
    cx, cy, R = 320, 165, 120
    pos = {}
    for i, v in enumerate(verts):
        ux, uy = _CIRCLE24[(i * 24 // max(n, 1)) % 24]
        pos[v] = (cx + ux * R // 100, cy - uy * R // 100)
    out = ["<svg viewBox='0 0 640 330' role='img' "
           "style='max-width:640px;width:100%;display:block;border:1px solid #ddd;"
           "border-radius:6px;background:#fafafa'>",
           "<defs><marker id='qlarr' viewBox='0 0 10 10' refX='9' refY='5' "
           "markerWidth='7' markerHeight='7' orient='auto-start-reverse'>"
           "<path d='M 0 0 L 10 5 L 0 10 z' fill='#444'/></marker></defs>"]
    pair_seen = {}
    for name in sorted(arrows):
        s, t = arrows[name]
        if s == t:                                   # loop: a small arc above
            x, y = pos[s]
            out.append("<path d='M %d %d C %d %d, %d %d, %d %d' fill='none' "
                       "stroke='#444' stroke-width='2' marker-end='url(#qlarr)'/>"
                       % (x - 10, y - 14, x - 26, y - 52, x + 26, y - 52,
                          x + 10, y - 14))
            out.append("<text x='%d' y='%d' font-size='15' font-style='italic' "
                       "text-anchor='middle' fill='#1c1c1c'>%s</text>"
                       % (x, y - 48, _esc(name)))
            continue
        (x1, y1), (x2, y2) = pos[s], pos[t]
        dx, dy = x2 - x1, y2 - y1
        L = max(_math_mod.isqrt(dx * dx + dy * dy), 1)
        # parallel arrows between the same unordered pair: fan them out with
        # alternating perpendicular offsets (0, -14, +14, -28, ...).
        pair = (min(s, t, key=str), max(s, t, key=str))
        k = pair_seen.get(pair, 0)
        pair_seen[pair] = k + 1
        shift = (k + 1) // 2 * 14 * (1 if k % 2 else -1) if k else 0
        px, py = -dy * shift // L, dx * shift // L        # perpendicular offset
        sx, sy = x1 + dx * 18 // L + px, y1 + dy * 18 // L + py
        ex, ey = x2 - dx * 18 // L + px, y2 - dy * 18 // L + py
        mx, my = (sx + ex) // 2 + px, (sy + ey) // 2 + py
        out.append("<path d='M %d %d Q %d %d %d %d' fill='none' stroke='#444' "
                   "stroke-width='2' marker-end='url(#qlarr)'/>"
                   % (sx, sy, mx, my, ex, ey))
        lx = (sx + ex) // 2 + px - (dy * 12 // L)
        ly = (sy + ey) // 2 + py + (dx * 12 // L) - 4
        out.append("<text x='%d' y='%d' font-size='15' font-style='italic' "
                   "text-anchor='middle' fill='#1c1c1c'>%s</text>"
                   % (lx, ly, _esc(name)))
    for v, (x, y) in pos.items():
        out.append("<circle cx='%d' cy='%d' r='15' fill='#fff' stroke='#3f51b5' "
                   "stroke-width='2'/>" % (x, y))
        out.append("<text x='%d' y='%d' font-size='14' text-anchor='middle' "
                   "fill='#1c1c1c'>%s</text>" % (x, y + 5, _esc(str(v))))
    out.append("</svg>")
    return "\n".join(out)


def _example_section(algebra):
    """'The example' header block: the quiver drawing + the defining relations +
    the algebra's one-line summary -- the reader knows INSTANTLY what is being
    computed (Marco 2026-07-28)."""
    if algebra is None:
        return []
    chunks = []
    svg = _quiver_svg(algebra)
    if svg:
        chunks.append(svg)
    rels = [str(r) for r in (getattr(algebra, "relations", None) or [])]
    if rels:
        chunks.append("<p><b>Relations:</b> <code>%s</code></p>"
                      % _esc(", ".join(rels)))
    else:
        arrows = getattr(getattr(algebra, "quiver", None), "arrows", None)
        if arrows is not None:
            chunks.append("<p><b>Relations:</b> none (hereditary).</p>")
    label = repr(algebra).splitlines()[0]
    chunks.append("<p>%s</p>" % _esc(label))
    return chunks


_JSON_NOTE = (
    "<p>This page pairs with a machine-readable record, <code>trace.json</code> "
    "(the “Report data (JSON)” download). Its structure:</p>"
    "<ul>"
    "<li><code>quiverlab_trace_schema</code> — integer version of this format;</li>"
    "<li><code>title</code> — the computation this record traces;</li>"
    "<li><code>references</code> — the same bibliography as the References "
    "section, as <code>[key, formatted]</code> pairs;</li>"
    "<li><code>events</code> — the complete ordered stream of worked steps. "
    "Each event is an object with a <code>type</code> field (e.g. "
    "<code>dispatch</code>, <code>resolution_term</code>, <code>rank_step</code>, "
    "<code>module_term</code>, <code>module_differential</code>, "
    "<code>ext_degree</code>, <code>step_note</code>, <code>result_dims</code>) "
    "plus that step's exact data — every matrix entry is exact "
    "(integers / rationals as strings), never floating point.</li>"
    "</ul>"
    "<p>The computation's summary results (one block per requested invariant, "
    "with dimensions, matrices, and citation keys) live in the separate "
    "<code>result.json</code> deliverable. HTML and JSON are the only report "
    "formats — this page is the complete human-readable record.</p>")


def _json_guide_html(json_guide):
    """The 'Reading the JSON record' appendix (Marco 2026-07-31): a table mapping each
    object THIS computation produced to its concrete path in ``result.json`` and a
    note. ``json_guide`` is a list of ``{object, path, note}`` dicts built per
    computation by the shared serializer; ``[]``/absent -> no section (tolerated for a
    pre-guide cached result)."""
    rows = [g for g in (json_guide or []) if isinstance(g, dict)]
    if not rows:
        return []
    out = ["<p><i>How to recover each object this computation produced from the "
           "machine-readable <code>result.json</code>. Every path is relative to the "
           "result object and uses only dot and <code>[\"key\"]</code> steps.</i></p>",
           "<table class='ql-guide'><tr><th>object</th><th>path</th>"
           "<th>note</th></tr>"]
    for g in rows:
        out.append("<tr><td>%s</td><td><code>%s</code></td><td>%s</td></tr>"
                   % (_esc(str(g.get("object", ""))), _esc(str(g.get("path", ""))),
                      _esc(str(g.get("note", "")))))
    out.append("</table>")
    return out


def _used_dispatches(events):
    """Collapse each CONSECUTIVE run of Dispatch events to its LAST entry --
    the resolution actually used. The engine may first record the bar attempt
    and then the reroute (auto -> Chouhy-Solotar); showing both read as a
    contradiction (Marco 2026-07-28: 'just print the one we use')."""
    used, run = [], None
    for e in events:
        if isinstance(e, Dispatch):
            run = e
        else:
            if run is not None:
                used.append(run)
                run = None
    if run is not None:
        used.append(run)
    return used


# --------------------------------------------------------------------------- #
# Fine-grained table of contents (Marco 2026-07-31): each top-level section lists,
# nested beneath it, the h3/h4 sub-headings it contains (a theory, a product, a
# degree). Built by scanning the already-rendered chunks for headings that carry a
# stable id -- render-side only, no data change.
# --------------------------------------------------------------------------- #
_SUBHEAD_RE = re.compile(r"<h([34]) id='([^']+)'>(.*?)</h\1>", re.DOTALL)
_MATH_BLOCK_RE = re.compile(r"<math\b[^>]*>.*?</math>", re.DOTALL)
_ANNOTATION_RE = re.compile(r"<annotation[^>]*>(.*?)</annotation>", re.DOTALL)
_TAG_RE = re.compile(r"<[^>]+>")


def _toc_label(heading_html):
    """A plain-text ToC label from a heading's inner HTML: a ``<math>`` run collapses
    to its x-tex annotation source, the remaining tags are stripped, entities decoded."""
    def _annot(m):
        a = _ANNOTATION_RE.search(m.group(0))
        return " %s " % a.group(1) if a else " "
    s = _MATH_BLOCK_RE.sub(_annot, heading_html)
    s = _TAG_RE.sub("", s)
    s = s.replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">")
    return " ".join(s.split())


def _sub_headings(section_html):
    """``[(level, id, label), ...]`` for every h3/h4 carrying an id, in document order."""
    return [(int(m.group(1)), m.group(2), _toc_label(m.group(3)))
            for m in _SUBHEAD_RE.finditer(section_html)]


def _nested_toc(subs):
    """A nested ``<ol>`` for a section: each h3 is a group, following h4s nest beneath
    the nearest preceding h3 (a stray h4 with no preceding h3 sits at the group level)."""
    out, i = ["<ol class='ql-toc-sub'>"], 0
    while i < len(subs):
        lvl, sid, label = subs[i]
        if lvl == 3:
            j, children = i + 1, []
            while j < len(subs) and subs[j][0] == 4:
                children.append(subs[j])
                j += 1
            out.append("<li><a href='#%s'>%s</a>" % (sid, _esc(label)))
            if children:
                out.append("<ol class='ql-toc-sub'>")
                out += ["<li><a href='#%s'>%s</a></li>" % (cid, _esc(clab))
                        for _l, cid, clab in children]
                out.append("</ol>")
            out.append("</li>")
            i = j
        else:
            out.append("<li><a href='#%s'>%s</a></li>" % (sid, _esc(label)))
            i += 1
    out.append("</ol>")
    return "".join(out)


def _build_toc(sections):
    body = ["<h2>Contents</h2><ol class='ql-toc'>"]
    for anchor, heading, chunks in sections:
        subs = _sub_headings("".join(chunks))
        body.append("<li><a href='#%s'>%s</a>" % (anchor, _esc(heading)))
        if subs:
            body.append(_nested_toc(subs))
        body.append("</li>")
    body.append("</ol>")
    return body


def render_html(events, title="", references=(), algebra=None, results=None,
                modules=(), json_guide=()):
    """The worked-steps report.

    ``results`` (optional) are the COMPUTED RESULT BLOCKS of the session -- every
    answer the page displayed, in either runner's shape. Passing them makes the
    saved report a complete record of what was computed, not only the worked steps
    of the one traced computation (Marco 2026-07-29). Omitted -> the page is
    byte-identical to before.

    ``modules`` (optional) is a sequence of ``(name, Module)`` the computation was
    about (``M``, and ``N`` when a second module was given): each is described in
    full -- Loewy series, top/socle, and every arrow's matrix. Descriptive only and
    individually guarded, so a module that cannot be described is skipped rather
    than sinking the report."""
    events = list(events)
    from quiverlab.errors import QuiverlabError
    from quiverlab.trace.events import ALL_EVENTS
    for e in events:
        if not isinstance(e, ALL_EVENTS):
            raise QuiverlabError(
                "render_html received a non-event object of type %r -- likely an "
                "unpacked (events, result) tuple from a trace_* helper; every "
                "element of the stream must be a trace event" % type(e).__name__)

    # ---- build the document as (anchor, heading, chunks) sections so a table
    # of contents can list, in order, everything the reader will scroll past.
    sections = []

    example = _example_section(algebra)
    if example:
        sections.append(("example", "The example", example))

    # The modules belong with the example: they are the PROBLEM statement, not an
    # answer. Then everything the session computed, then the worked steps that
    # justify it (Marco 2026-07-29).
    described = _describe_modules(modules)
    mod_sec = _modules_section_html(described)
    if mod_sec:
        sections.append(("modules", "The modules", mod_sec))

    from quiverlab.trace.results_html import results_section
    computed = results_section(results)
    if computed:
        sections.append(("computed", "Computed results", computed))

    used = _used_dispatches(events)
    if used:
        chunks = []
        for e in used:
            chunks.append("<p><b>Resolution used:</b> %s<br><i>%s</i><br>"
                          "defining relations: %d</p>"
                          % (_esc(e.route), _esc(e.reason), e.n_relations))
        sections.append(("resolution", "Resolution", chunks))

    objects, note = compute_algebra_objects(algebra)
    pi = _pi_section_html(objects, note)
    if pi:
        # _pi_section_html emits its own h2; keep its body, retitle via TOC.
        chunks = pi[1:]
        chunks.append("<p class='ql-cite'>Loewy series and the projectives/"
                      "injectives P<sub>v</sub>, I<sub>v</sub> follow [ASS2006] "
                      "(see References).</p>")
        sections.append(("projectives-injectives",
                         "The projectives and injectives of A", chunks))

    terms = {e.degree: e for e in events if isinstance(e, ResolutionTerm)}
    ranks = {e.degree: e for e in events if isinstance(e, RankStep)}
    if terms:
        chunks = ["<p><i>C<sub>n</sub> is the degree-n term of the resolution "
                  "named above; the matrices below are its differentials "
                  "(rows: target basis, columns: source basis).</i></p>"]
        echo = _MatrixEcho()
        # Plan 35 UNIT 2 (deliverable 2 + review fix): the ordered basis of each HH
        # (co)chain term, reconstructed with UNIT 1's builders. report_side (all RankSteps
        # of one HH run share it) picks cochain vs chain. A term with NO corners is a bar
        # (co)chain space (rebuilt from the unit-adapted algebra); a term WITH corners is a
        # Chouhy-Solotar term (rebuilt from the deterministic CS resolution, only when the
        # run actually used CS). Both are shown only when the reconstructed length equals
        # the recorded term dim, so the differential's columns/rows are the ones enumerated.
        report_side = next((r.side for r in ranks.values()), None)
        kind_side = (("cochain" if report_side == "cochain" else "chain")
                     if report_side else None)
        m_bar = labels_bar = cs_res = cs_labels = None
        if algebra is not None and kind_side is not None:
            try:
                from quiverlab.hochschild import basis_reps as BR
                _AU = algebra.unit_adapted()
                m_bar, labels_bar = _AU.dim, BR.labels_of(_AU)
            except Exception:
                m_bar = labels_bar = None
            # Build the CS resolution ONCE, only when a term carries corners (the run used
            # Chouhy-Solotar). The report re-derives it (the `_product_object` pattern); a
            # non-admissible / structure-constants algebra fails the build -> omit.
            if any(getattr(terms[k], "corners", None) is not None for k in terms):
                try:
                    from quiverlab.hochschild import basis_reps as BR
                    from quiverlab.resolutions_cs.build import reduction_system_of
                    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
                    cs_res = ChouhySolotarResolution(algebra, reduction_system_of(algebra),
                                                     max_degree=max(terms) + 2)
                    cs_labels = BR.labels_of(cs_res.ar.A)
                except Exception:
                    cs_res = cs_labels = None
        for n in sorted(terms):
            t = terms[n]
            chunks.append("<h3 id='resstep-deg-%d'>Degree %d</h3><p>Term with %d "
                          "generators (dim C<sub>%d</sub> = %d).</p>"
                          % (n, n, t.n_generators, n, t.collapsed_dim))
            summ = _term_summands_html(t, n)
            if summ:
                chunks.append(summ)
            corners = getattr(t, "corners", None)
            if kind_side and corners is None and m_bar is not None:
                chunks.extend(_bar_term_basis_chunks(m_bar, labels_bar, n,
                                                     t.collapsed_dim, kind_side))
            elif kind_side and corners is not None and cs_res is not None:
                side_key = "coh" if kind_side == "cochain" else "hom"
                chunks.extend(_cs_term_basis_chunks(cs_res, cs_labels, n,
                                                    t.collapsed_dim, side_key, kind_side))
            if n in ranks:
                rs = ranks[n]
                if rs.side == "cochain":
                    sym = "d^{%d}" % n
                    arrow = r"d^{%d} : C^{%d} \to C^{%d}" % (n, n, n + 1)
                else:
                    sym = "b_{%d}" % n
                    arrow = r"b_{%d} : C_{%d} \to C_{%d}" % (n, n, max(n - 1, 0))
                prior = echo.label_for(rs, sym)
                # The arrow declaration already names the map, so the rank clause
                # stays the bare "rank = k" both renderers share.
                chunks.append(_math(r"%s,\qquad \operatorname{rank} = %d"
                                    % (arrow, rs.rank)))
                if prior is not None:
                    chunks.append(_math(r"%s = %s" % (sym, prior)))
                    chunks.append("<p class='ql-note'>(the same matrix as above; "
                                  "not repeated)</p>")
                else:
                    chunks.append(_event_grid(rs, label=sym))
        sections.append(("resolution-steps", "Worked resolution steps", chunks))

    mod = _module_steps_html(events)
    if mod:
        sections.append(("module-steps", "Worked module steps", mod[1:]))

    # Plan 35: the HH-product chapter (cup / cap / bracket / connes_b). The stream
    # carries ProductStep events; the section names the product and shows the
    # definition + per-bidegree tables/matrices. The HH result (a ResultDims the
    # products builder injects) still renders below via the shared path.
    prod = _products_html(events)
    if prod:
        pkind = next(e.kind for e in events if isinstance(e, ProductStep))
        sections.append(("product-steps",
                         _PRODUCT_TITLE.get(pkind, "HH products"), prod))

    # The (co)homology Result: prefer the AUTHORITATIVE dims the engine returned (a
    # ResultDims event, injected by writer.py) so the line carries the engine's numbers
    # AND the correct HH^/HH_ variance; fall back to derive_dims when absent.
    # ...unless the Computed results section already carries that very table: it
    # would then be the SAME numbers printed twice, under a heading ("Result") that
    # says nothing, when everything on the page is a result (Marco 2026-07-29). The
    # section is named for what it holds -- Hochschild homology / cohomology.
    rd = next((e for e in events if isinstance(e, ResultDims)), None)
    kind = rd.kind if rd is not None else _dims_kind(events)
    dims = list(rd.dims) if rd is not None else derive_dims(events)
    if dims and not _results_carry_hh(results, kind):
        chunks = []
        if kind in ("HH^", "HH_"):              # typing statement at the top (Marco)
            chunks.append(hh_typing_html(kind, _route_of_dispatch(used)))
        chunks.append(_dims_table(_hh_row_label(kind), dims))
        if rd is not None and rd.note:
            chunks.append("<p><i>%s</i></p>" % _esc(rd.note))
        sections.append(("result", _hh_heading(kind), chunks))

    sections.append(("json-record", "The JSON record", [_JSON_NOTE]))

    # Addendum (Marco 2026-07-31): "Reading the JSON record" -- concrete path recipes
    # for every object THIS computation produced. Empty on a pre-guide cache (tolerated).
    guide = _json_guide_html(json_guide)
    if guide:
        sections.append(("json-guide", "Reading the JSON record", guide))

    if references:
        chunks = ["<ol>"]
        for key, entry in references:
            chunks.append("<li>[%s] %s</li>" % (_esc(key), _esc(entry)))
        chunks.append("</ol>")
        sections.append(("references", "References", chunks))

    # ---- assemble: header, contents, sections.
    body = ["<!doctype html><html><head><meta charset='utf-8'>",
            "<meta name='viewport' content='width=device-width, initial-scale=1'>",
            _STYLE,
            "<title>Worked steps: %s</title></head><body>" % _esc(title),
            "<h1>Worked steps: %s</h1>" % _esc(title),
            "<p class='ql-hint'>The complete worked-steps record, exactly as "
            "computed. Deliverables are this HTML page and the JSON records "
            "described at the end — mathematics is typeset as MathML.</p>"]
    if len(sections) > 1:
        body.extend(_build_toc(sections))
    for anchor, heading, chunks in sections:
        body.append("<h2 id='%s'>%s</h2>" % (anchor, _esc(heading)))
        body.extend(chunks)
    body.append("</body></html>")
    return "\n".join(body) + "\n"


def _esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


# --------------------------------------------------------------------------- #
# TeX-source math builders (the worked-steps math, embedded in the MathML x-tex
# ``<annotation>`` and typeset by the converter above). Formerly shared with the
# deleted LaTeX-document renderer; this module now owns them. Float-free: every
# number comes from event fields (ints/strings).
# --------------------------------------------------------------------------- #


def _pmatrix(rs):
    """The matrix as `pmatrix` TeX source, or a `\\text{...}` note when the recorder
    dropped the body (the 250k memory backstop). Works for any event carrying
    (elided, matrix, nrows, ncols, note): RankStep, ModuleDifferential, ExtDegree.

    A zero-DIMENSIONAL matrix (no rows or no columns -- e.g. every zero Ext/Tor
    differential) renders as the symbol ``0``: an empty
    ``\\begin{pmatrix}\\end{pmatrix}`` typesets as a stray ``()``."""
    if rs.elided or rs.matrix is None:
        return r"\text{%s}" % _tex_escape(rs.note)
    if rs.nrows == 0 or rs.ncols == 0:
        return "0"
    rows = r" \\ ".join(" & ".join(rs.matrix[i][j] for j in range(rs.ncols))
                        for i in range(rs.nrows))
    return r"\begin{pmatrix} %s \end{pmatrix}" % rows


def oplus_tex(summands, sym):
    """A direct-sum ``P_{1}^{2} \\oplus P_{3}`` from a vertex list with repetition;
    ``0`` for an empty term."""
    if not summands:
        return "0"
    groups = []
    for v in summands:
        for g in groups:
            if g[0] == v:
                g[1] += 1
                break
        else:
            groups.append([v, 1])
    return r" \oplus ".join(
        "%s_{%s}" % (sym, v) if c == 1 else "%s_{%s}^{%d}" % (sym, v, c)
        for v, c in groups)


def factor_stack_tex(dimvec):
    """One semisimple Loewy layer as a stacked ``\\begin{matrix} S_1 \\\\ S_2 ...``
    column of composition factors (with multiplicity)."""
    factors = []
    for v, m in dimvec.items():
        if m <= 0:
            continue
        factors.extend(["S_{%s}" % v] * m)
    if not factors:
        return "0"
    return r"\begin{matrix} %s \end{matrix}" % r" \\ ".join(factors)


def _dimvec_tex(dv):
    return "(" + ",\\ ".join("%s{:}%s" % (k, val) for k, val in dv.items()) + ")"


def _tex_escape(s):
    # The worked-steps math source is embedded verbatim in the MathML x-tex
    # annotation (and shown as source in the plain report), so every special char
    # must be escaped. Order matters: '\' is escaped FIRST to a placeholder (not
    # directly to \textbackslash{}, whose own braces would then be re-escaped by the
    # {}/} steps); the brace steps run BEFORE '^'/'~' (whose replacements introduce
    # their own {}); the placeholder is expanded LAST so nothing double-escapes.
    sentinel = "\x00"
    s = s.replace("\\", sentinel)
    for a, b in (("{", r"\{"), ("}", r"\}"), ("&", r"\&"), ("_", r"\_"),
                 ("#", r"\#"), ("%", r"\%"), ("$", r"\$"),
                 ("^", r"\textasciicircum{}"), ("~", r"\textasciitilde{}"),
                 ("<", r"\textless{}"), (">", r"\textgreater{}"),
                 ("|", r"\textbar{}")):
        s = s.replace(a, b)
    return s.replace(sentinel, r"\textbackslash{}")
