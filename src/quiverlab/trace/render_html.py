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
    ExtDegree, StepNote, ResultDims,
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
    ".ql-eq{margin:.7em 0;overflow-x:auto}"
    ".ql-hint{background:#eef3ff;border:1px solid #9db8e8;border-radius:6px;"
    "padding:.6em .85em;margin:0 0 1.4em;font-family:sans-serif;font-size:.92em;"
    "color:#123}"
    "ol,ul{padding-left:1.4em}"
    "@media print{.ql-hint{display:none}"
    "body{max-width:none;margin:0;font-size:11pt}"
    "h2,.ql-eq,math{break-inside:avoid;page-break-inside:avoid}"
    "h1,h2{break-after:avoid}}"
    "@page{margin:2cm}"
    "</style>")

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


def _math(expr):
    """A display equation: typeset presentation MathML with the LaTeX source kept
    verbatim in an x-tex ``<annotation>`` (self-contained, no JavaScript). On any
    conversion surprise, fall back to ``<mtext>`` of the escaped source so the
    document stays well-formed and the source is still readable."""
    try:
        pres = _tex_to_mathml_body(expr) or ("<mtext>%s</mtext>" % _esc(expr))
    except Exception:
        pres = "<mtext>%s</mtext>" % _esc(expr)
    return ('<div class="ql-eq"><math display="block"><semantics><mrow>%s</mrow>'
            '<annotation encoding="application/x-tex">%s</annotation>'
            '</semantics></math></div>' % (pres, _esc(expr)))


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
        _math(r"%s = %s" % (shown, _pmatrix(e))),
        _math(r"%s = %s,\qquad \dim = %s = %d"
              % (label, quotient, arith, e.result_dim)),
    ]


def _module_steps_html(events):
    mods = [e for e in events
            if isinstance(e, (ModuleTerm, ModuleDifferential, ExtDegree, StepNote))]
    if not mods:
        return []
    out = ["<h2>Worked module steps</h2>", "<p><i>%s</i></p>" % _esc(ELISION_PREAMBLE)]
    step_no = 0
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
            out.append(_math(r"%s : %s \to %s \qquad %s = %s"
                             % (e.symbol, dom, cod, e.symbol, _pmatrix(e))))
        elif isinstance(e, ExtDegree):
            out.extend(_ext_degree_html(e))
    # ALL Ext/Tor runs (not only the first): an Ext run then a Tor run both appear,
    # each with the correct subscript/superscript (Plan 34 MINOR).
    runs = ext_result_runs(events)
    if runs:
        out.append("<h2>Result</h2>")
        for op, dims in runs:
            sep = "_" if op == "Tor" else "^"
            out.append(_math(",\\quad ".join(
                r"\operatorname{%s}%s{%d} = %d" % (op, sep, i, d)
                for i, d in enumerate(dims))))
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


def render_html(events, title="", references=(), algebra=None):
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
        for n in sorted(terms):
            t = terms[n]
            chunks.append("<h3>Degree %d</h3><p>Term with %d generators "
                          "(dim C<sub>%d</sub> = %d).</p>"
                          % (n, t.n_generators, n, t.collapsed_dim))
            if n in ranks:
                rs = ranks[n]
                if rs.side == "cochain":
                    sym = "d^{%d}" % n
                    arrow = r"d^{%d} : C^{%d} \to C^{%d}" % (n, n, n + 1)
                else:
                    sym = "b_{%d}" % n
                    arrow = r"b_{%d} : C_{%d} \to C_{%d}" % (n, n, max(n - 1, 0))
                chunks.append(_math(r"%s,\qquad %s = %s \qquad "
                                    r"\operatorname{rank} = %d"
                                    % (arrow, sym, _pmatrix(rs), rs.rank)))
        sections.append(("resolution-steps", "Worked resolution steps", chunks))

    mod = _module_steps_html(events)
    if mod:
        sections.append(("module-steps", "Worked module steps", mod[1:]))

    # The (co)homology Result: prefer the AUTHORITATIVE dims the engine returned (a
    # ResultDims event, injected by writer.py) so the line carries the engine's numbers
    # AND the correct HH^/HH_ variance; fall back to derive_dims when absent.
    rd = next((e for e in events if isinstance(e, ResultDims)), None)
    if rd is not None:
        cells = ",\\quad ".join(r"%s{%d} = %d" % (rd.kind, i, d)
                                for i, d in enumerate(rd.dims))
        chunks = [_math(cells)]
        if rd.note:
            chunks.append("<p><i>%s</i></p>" % _esc(rd.note))
        sections.append(("result", "Result", chunks))
    else:
        dims = derive_dims(events)
        if dims:
            kind = _dims_kind(events)
            cells = ",\\quad ".join(r"%s{%d} = %d" % (kind, i, d)
                                    for i, d in enumerate(dims))
            sections.append(("result", "Result", [_math(cells)]))

    sections.append(("json-record", "The JSON record", [_JSON_NOTE]))

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
        body.append("<h2>Contents</h2><ol class='ql-toc'>")
        for anchor, heading, _ in sections:
            body.append("<li><a href='#%s'>%s</a></li>" % (anchor, _esc(heading)))
        body.append("</ol>")
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
