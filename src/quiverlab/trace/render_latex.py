"""LaTeX worked-steps renderer -> standalone article compiled to PDF (spec §3.8).

Matrices as pmatrix; resulting dims DERIVED from the recorded events (never echoed,
via render_text.derive_dims); References as thebibliography. Float-free: every
number comes from event fields (ints/strings).

Shared helpers: the resulting-dimension derivation (`derive_dims`) and the
cohomology/homology kind detection (`_dims_kind`) are REUSED from Task 9's
render_text -- they are not redefined here. The matrix->pmatrix TeX helper
(`_pmatrix`) is DEFINED here (LaTeX is the primary math output) and imported by
render_html so both renderers emit identical `pmatrix` source from one place."""
from quiverlab.trace.events import Dispatch, ResolutionTerm, RankStep
from quiverlab.trace.render_text import derive_dims, _dims_kind


def _pmatrix(rs):
    """The RankStep matrix as `pmatrix` TeX source (shared by render_latex and
    render_html), or a `\\text{...}` note when the body was elided."""
    if rs.elided or rs.matrix is None:
        return r"\text{%s}" % rs.note
    rows = r" \\ ".join(" & ".join(rs.matrix[i][j] for j in range(rs.ncols))
                        for i in range(rs.nrows))
    return r"\begin{pmatrix} %s \end{pmatrix}" % rows


def render_latex(events, title="", references=()):
    events = list(events)
    out = [r"\documentclass{article}", r"\usepackage{amsmath}",
           r"\begin{document}", r"\section*{Worked steps: %s}" % _tex_escape(title)]
    for e in events:
        if isinstance(e, Dispatch):
            out.append(r"\noindent\textbf{Chosen resolution:} %s\\" % _tex_escape(e.route))
            out.append(r"\textit{%s}\\" % _tex_escape(e.reason))
            out.append(r"defining relations: %d" % e.n_relations)
    terms = {e.degree: e for e in events if isinstance(e, ResolutionTerm)}
    ranks = {e.degree: e for e in events if isinstance(e, RankStep)}
    for n in sorted(terms):
        t = terms[n]
        out.append(r"\subsection*{Degree %d}" % n)
        out.append(r"Term with %d generators ($\dim C = %d$)." % (t.n_generators, t.collapsed_dim))
        if n in ranks:
            rs = ranks[n]
            sym = "d^{%d}" % n if rs.side == "cochain" else "b_{%d}" % n
            out.append(r"\[ %s = %s \qquad \operatorname{rank} = %d \]"
                       % (sym, _pmatrix(rs), rs.rank))
    dims = derive_dims(events)
    if dims:
        kind = _dims_kind(events)
        cells = r",\quad ".join(r"%s{%d} = %d" % (kind, i, d) for i, d in enumerate(dims))
        out.append(r"\subsection*{Result}")
        out.append(r"\[ %s \]" % cells)
    if references:
        out.append(r"\begin{thebibliography}{9}")
        for key, entry in references:
            # Registry/bibtex keys can carry '_' (e.g. "chouhy_solotar",
            # "Luebeck_ConwayPolynomials"), so sanitize to [A-Za-z0-9] before
            # interpolating into \bibitem{} -- these labels are never \cite'd, so
            # stripping is safe and keeps real PDF compilation robust.
            out.append(r"\bibitem{%s} %s" % (_bibkey(key), _tex_escape(entry)))
        out.append(r"\end{thebibliography}")
    out.append(r"\end{document}")
    return "\n".join(out) + "\n"


def _bibkey(key):
    """A \\bibitem label restricted to [A-Za-z0-9] (LaTeX-safe, never \\cite'd)."""
    return "".join(c for c in str(key) if c.isalnum() and c.isascii())


def _tex_escape(s):
    # LaTeX->PDF is the primary output, so every special char must be escaped.
    # Order matters: '\' is escaped FIRST to a placeholder (not directly to
    # \textbackslash{}, whose own braces would then be re-escaped by the {}/}
    # steps); the brace steps run BEFORE '^'/'~' (whose replacements introduce
    # their own {}); the placeholder is expanded LAST so nothing double-escapes.
    sentinel = "\x00"
    s = s.replace("\\", sentinel)
    for a, b in (("{", r"\{"), ("}", r"\}"), ("&", r"\&"), ("_", r"\_"),
                 ("#", r"\#"), ("%", r"\%"), ("$", r"\$"),
                 ("^", r"\textasciicircum{}"), ("~", r"\textasciitilde{}")):
        s = s.replace(a, b)
    return s.replace(sentinel, r"\textbackslash{}")
