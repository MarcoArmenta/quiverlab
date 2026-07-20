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
from quiverlab.trace.events import Dispatch, ResolutionTerm, RankStep
from quiverlab.trace.render_text import derive_dims, _dims_kind
from quiverlab.trace.render_latex import _pmatrix

# Inline-only styling (no external stylesheet); keeps the file fully self-contained.
_STYLE = ("<style>body{font-family:sans-serif}"
          "pre{background:#f4f4f4;padding:6px;overflow-x:auto}</style>")


def _math(expr):
    """Show the TeX SOURCE as text (no MathJax); escape HTML metachars only."""
    return "<pre><code>%s</code></pre>" % _esc(expr)


def render_html(events, title="", references=()):
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
