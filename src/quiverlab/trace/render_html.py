"""HTML worked-steps renderer -- NO JavaScript, NO external resources (Marco's
decision). LaTeX->PDF (pdflatex/tectonic) is the PRIMARY, typeset math output;
this HTML is the self-contained, offline fallback for when no LaTeX toolchain is on
PATH. Math is shown as readable *TeX source* (inside <pre><code>), NOT typeset --
there is no MathJax, no CDN <script>, and no external <link>, so the file renders
identically in any browser with the network off. Float-free: all numbers come from
event fields (ints/strings).

Shared helpers: `derive_dims` + `_dims_kind` are REUSED from Task 9's render_text,
and the matrix->pmatrix TeX helper `_pmatrix` is REUSED from render_latex -- so both
renderers emit identical `pmatrix` source and identical resulting dimensions from a
single definition, never duplicated here."""
from quiverlab.trace.events import (
    Dispatch, ResolutionTerm, RankStep, ModuleTerm, ModuleDifferential,
    ExtDegree, StepNote,
)
from quiverlab.trace.render_text import (
    derive_dims, _dims_kind, compute_algebra_objects, ext_result_dims,
    ELISION_PREAMBLE,
)
from quiverlab.trace.render_latex import (
    _pmatrix, oplus_tex, factor_stack_tex, _dimvec_tex,
)

# Inline-only styling (no external stylesheet); keeps the file fully self-contained.
_STYLE = ("<style>body{font-family:sans-serif}"
          "pre{background:#f4f4f4;padding:6px;overflow-x:auto}</style>")


def _math(expr):
    """Show the TeX SOURCE as text (no MathJax); escape HTML metachars only."""
    return "<pre><code>%s</code></pre>" % _esc(expr)


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


def _module_steps_html(events):
    mods = [e for e in events
            if isinstance(e, (ModuleTerm, ModuleDifferential, ExtDegree, StepNote))]
    if not mods:
        return []
    out = ["<h2>Worked module steps</h2>", "<p><i>%s</i></p>" % _esc(ELISION_PREAMBLE)]
    for e in mods:
        if isinstance(e, StepNote):
            out.append("<p>%s%s</p>" % (
                _esc(e.text),
                "<br><i>%s</i>" % _esc(e.detail) if e.detail else ""))
        elif isinstance(e, ModuleTerm):
            what = "Q_{%d}" % e.degree if e.sym == "P" else "E^{%d}" % e.degree
            out.append(_math(r"%s = %s \qquad \dim = %d"
                             % (what, oplus_tex(e.summands, e.sym), e.dim)))
        elif isinstance(e, ModuleDifferential):
            cod = "M" if e.cod_is_module else oplus_tex(e.cod_summands, e.sym)
            out.append(_math(r"%s : %s \to %s \qquad %s = %s"
                             % (e.symbol, oplus_tex(e.dom_summands, e.sym), cod,
                                e.symbol, _pmatrix(e))))
        elif isinstance(e, ExtDegree):
            out.append("<p>%s<sup>%d</sup>: dim Hom = %d, rank = %d (prev %d)</p>"
                       % (_esc(e.op), e.degree, e.space_dim, e.rank_here, e.rank_prev))
            out.append(_math(r"\delta^{%d} = %s \qquad %s^{%d} = %d"
                             % (e.degree, _pmatrix(e), e.op, e.degree, e.result_dim)))
    dims, op = ext_result_dims(events)
    if dims:
        out.append("<h2>Result</h2>" + _math(",\\quad ".join(
            r"%s^{%d} = %d" % (op, i, d) for i, d in enumerate(dims))))
    return out


def render_html(events, title="", references=(), algebra=None):
    events = list(events)
    body = ["<!doctype html><html><head><meta charset='utf-8'>", _STYLE,
            "<title>Worked steps: %s</title></head><body>" % _esc(title),
            "<h1>Worked steps: %s</h1>" % _esc(title),
            "<p><i>Math is shown as TeX source (no JavaScript); compile the PDF with "
            "pdflatex/tectonic for typeset output.</i></p>"]
    for e in events:
        if isinstance(e, Dispatch):
            body.append("<p><b>Chosen resolution:</b> %s<br><i>%s</i><br>"
                        "defining relations: %d</p>" % (_esc(e.route), _esc(e.reason), e.n_relations))
    objects, note = compute_algebra_objects(algebra)
    body.extend(_pi_section_html(objects, note))
    terms = {e.degree: e for e in events if isinstance(e, ResolutionTerm)}
    ranks = {e.degree: e for e in events if isinstance(e, RankStep)}
    for n in sorted(terms):
        t = terms[n]
        body.append("<h2>Degree %d</h2><p>Term with %d generators (dim C = %d).</p>"
                    % (n, t.n_generators, t.collapsed_dim))
        if n in ranks:
            rs = ranks[n]
            sym = "d^{%d}" % n if rs.side == "cochain" else "b_{%d}" % n
            body.append(_math(r"%s = %s \qquad \operatorname{rank} = %d"
                              % (sym, _pmatrix(rs), rs.rank)))
    body.extend(_module_steps_html(events))
    dims = derive_dims(events)
    if dims:
        kind = _dims_kind(events)
        cells = ",\\quad ".join(r"%s{%d} = %d" % (kind, i, d) for i, d in enumerate(dims))
        body.append("<h2>Result</h2>" + _math(cells))
    if references:
        body.append("<h2>References</h2><ol>")
        for key, entry in references:
            body.append("<li>[%s] %s</li>" % (_esc(key), _esc(entry)))
        body.append("</ol>")
    body.append("</body></html>")
    return "\n".join(body) + "\n"


def _esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
