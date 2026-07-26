"""LaTeX worked-steps renderer -> standalone article compiled to PDF (spec §3.8).

Matrices as pmatrix; resulting dims DERIVED from the recorded events (never echoed,
via render_text.derive_dims); References as thebibliography. Float-free: every
number comes from event fields (ints/strings).

Shared helpers: the resulting-dimension derivation (`derive_dims`) and the
cohomology/homology kind detection (`_dims_kind`) are REUSED from Task 9's
render_text -- they are not redefined here. The matrix->pmatrix TeX helper
(`_pmatrix`) is DEFINED here (LaTeX is the primary math output) and imported by
render_html so both renderers emit identical `pmatrix` source from one place."""
from quiverlab.trace.events import (
    Dispatch, ResolutionTerm, RankStep, ModuleTerm, ModuleDifferential,
    ExtDegree, StepNote,
)
from quiverlab.trace.render_text import (
    derive_dims, _dims_kind, compute_algebra_objects, ext_result_dims,
    ELISION_PREAMBLE,
)


def _pmatrix(rs):
    """The RankStep matrix as `pmatrix` TeX source (shared by render_latex and
    render_html), or a `\\text{...}` note when the body was elided. Works for any
    event carrying (elided, matrix, nrows, ncols, note): RankStep, ModuleDifferential,
    ExtDegree."""
    if rs.elided or rs.matrix is None:
        return r"\text{%s}" % _tex_escape(rs.note)
    rows = r" \\ ".join(" & ".join(rs.matrix[i][j] for j in range(rs.ncols))
                        for i in range(rs.nrows))
    return r"\begin{pmatrix} %s \end{pmatrix}" % rows


def oplus_tex(summands, sym):
    """A direct-sum ``P_{1}^{2} \\oplus P_{3}`` from a vertex list with repetition;
    ``0`` for an empty term. Shared by the LaTeX and (TeX-source) HTML renderers."""
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


def _pi_section_latex(objects, note):
    """The "projectives and injectives of A" subsection (Marco #4): each vertex's
    P_v and I_v with dimension vector and stacked Loewy layers (simples omitted)."""
    if not objects:
        if note:
            return [r"\subsection*{The projectives and injectives of $A$}",
                    r"\noindent (unavailable: %s)" % _tex_escape(note)]
        return []
    out = [r"\subsection*{The projectives and injectives of $A$}",
           r"\noindent Loewy layers stacked top to bottom; the simples $S_v$ are "
           r"omitted.\\"]
    for row in objects:
        v = _tex_escape(row["vertex"])
        for sym in ("P", "I"):
            d = row[sym]
            layers = r" \;\big|\; ".join(factor_stack_tex(L) for L in d["layers"]) \
                or "0"
            out.append(r"\[ %s_{%s} = %s \qquad \dim = %d,\ \underline{\dim} = %s \]"
                       % (sym, v, layers, d["dim"], _dimvec_tex(d["dimvec"])))
    return out


def _dimvec_tex(dv):
    return "(" + ",\\ ".join("%s{:}%s" % (k, val) for k, val in dv.items()) + ")"


def _module_steps_latex(events):
    """The module worked-steps (ModuleTerm / ModuleDifferential / ExtDegree /
    StepNote) as a LaTeX subsection, in emission order; empty when there are none."""
    mods = [e for e in events
            if isinstance(e, (ModuleTerm, ModuleDifferential, ExtDegree, StepNote))]
    if not mods:
        return []
    out = [r"\subsection*{Worked module steps}",
           r"\noindent %s\\" % _tex_escape(ELISION_PREAMBLE)]
    for e in mods:
        if isinstance(e, StepNote):
            out.append(r"\noindent %s\\" % _tex_escape(e.text))
            if e.detail:
                out.append(r"\textit{%s}\\" % _tex_escape(e.detail))
        elif isinstance(e, ModuleTerm):
            what = "Q_{%d}" % e.degree if e.sym == "P" else "E^{%d}" % e.degree
            out.append(r"\[ %s = %s \qquad \dim = %d \]"
                       % (what, oplus_tex(e.summands, e.sym), e.dim))
        elif isinstance(e, ModuleDifferential):
            cod = "M" if e.cod_is_module else oplus_tex(e.cod_summands, e.sym)
            out.append(r"\[ %s : %s \to %s \qquad %s = %s \]"
                       % (e.symbol, oplus_tex(e.dom_summands, e.sym), cod,
                          e.symbol, _pmatrix(e)))
        elif isinstance(e, ExtDegree):
            out.append(r"\noindent $%s^{%d}$: $\dim\operatorname{Hom} = %d$, "
                       r"$\operatorname{rank} = %d$ (prev %d).\\"
                       % (e.op, e.degree, e.space_dim, e.rank_here, e.rank_prev))
            out.append(r"\[ \delta^{%d} = %s \qquad %s^{%d} = %d \]"
                       % (e.degree, _pmatrix(e), e.op, e.degree, e.result_dim))
    dims, op = ext_result_dims(events)
    if dims:
        out.append(r"\subsection*{Result}")
        out.append(r"\[ %s \]" % r",\quad ".join(
            r"%s^{%d} = %d" % (op, i, d) for i, d in enumerate(dims)))
    return out


def render_latex(events, title="", references=(), algebra=None):
    events = list(events)
    out = [r"\documentclass{article}", r"\usepackage{amsmath}",
           r"\begin{document}", r"\section*{Worked steps: %s}" % _tex_escape(title)]
    for e in events:
        if isinstance(e, Dispatch):
            out.append(r"\noindent\textbf{Chosen resolution:} %s\\" % _tex_escape(e.route))
            out.append(r"\textit{%s}\\" % _tex_escape(e.reason))
            out.append(r"defining relations: %d" % e.n_relations)
    objects, note = compute_algebra_objects(algebra)
    out.extend(_pi_section_latex(objects, note))
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
    out.extend(_module_steps_latex(events))
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
