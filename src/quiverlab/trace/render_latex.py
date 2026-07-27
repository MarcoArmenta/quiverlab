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
)


MATRIX_COL_DEFAULT = 10          # amsmath's default MaxMatrixCols ceiling for pmatrix
# The PDF is a PAGE-BOUNDED homework document (Plan 34 artifact contract): the events
# carry every matrix in full (recorder), the HTML/JSON reports show them in full, and the
# PDF renderer -- HERE and in hpc/report -- is the only one that scales/elides wide
# matrices for the page. A matrix wider/taller than this in either dimension is elided to
# a stated shape (pointing at the complete HTML/JSON report); a matrix of 11..25 columns
# is scaled to the line width; up to 10 columns it is a plain pmatrix.
MATRIX_MAX_DIM = 25
# The one-sentence report policy, stated once in the module-steps preamble (Plan 34;
# replaces the old recorder-elision sentence, which no longer matches the contract).
_LATEX_MATRIX_POLICY = (
    "Matrices up to 25 columns are typeset in full (scaled to the page width when "
    "wider than the default), and a matrix past 25 rows or columns is elided here to "
    "its shape and shown in full in the accompanying HTML/JSON report; every step "
    "still appears.")
# A graphicx \resizebox macro that SHRINKS an over-wide matrix box to the line width but
# never ENLARGES a small one (the \ifdim guard) -- so a homework-scale matrix that already
# fits stays crisp, and an 11..25-column differential is scaled down onto the page instead
# of aborting the compile ("Extra alignment tab") or running off it. Defined once in the
# preamble of every worked-steps document (matrix_preamble_lines).
_QLMAT_DEF = (r"\newcommand{\qlmat}[1]{\resizebox{\ifdim\width>\linewidth"
              r"\linewidth\else\width\fi}{!}{$#1$}}")


def matrix_preamble_lines(max_cols):
    """The preamble additions that make a LaTeX document whose widest TYPESET matrix has
    ``max_cols`` columns compile AND stay on the page (spec §3.8, Plan 34):

      * ``\\usepackage{graphicx}`` + the ``\\qlmat`` shrink-to-width macro -- always,
        since even the projectives/injectives Loewy displays are wrapped in it;
      * ``\\setcounter{MaxMatrixCols}{max_cols}`` when ``max_cols`` exceeds amsmath's
        default of 10, or pdflatex/tectonic aborts with "Extra alignment tab has been
        changed to \\cr" the moment a pmatrix has an 11th column.

    ``max_cols`` is the widest matrix ACTUALLY typeset as a pmatrix (i.e. <= 25 columns;
    wider ones are elided to text, so they never need the raised ceiling). Shared by
    ``trace.render_latex`` and ``hpc.report`` so both LaTeX paths size matrices
    identically. Deterministic (constant strings + one int)."""
    lines = [r"\usepackage{graphicx}", _QLMAT_DEF]
    if max_cols > MATRIX_COL_DEFAULT:
        lines.append(r"\setcounter{MaxMatrixCols}{%d}" % max_cols)
    return lines


def _too_wide(ev):
    """A recorded (full) matrix too big to typeset on the page: more than
    ``MATRIX_MAX_DIM`` rows OR columns. It is STATED-elided in the PDF (shown in full in
    the HTML/JSON report), so it does not raise MaxMatrixCols and is not \\qlmat-scaled."""
    return (getattr(ev, "matrix", None) is not None and not getattr(ev, "elided", False)
            and ((getattr(ev, "ncols", 0) or 0) > MATRIX_MAX_DIM
                 or (getattr(ev, "nrows", 0) or 0) > MATRIX_MAX_DIM))


def _events_max_cols(events):
    """The widest matrix that will actually be TYPESET as a pmatrix (<= 25 columns) --
    what ``\\setcounter{MaxMatrixCols}`` must clear so the alignment never aborts.
    Record-elided (``\\text{...}``) and page-elided (>25) matrices contribute nothing."""
    best = 0
    for e in events:
        if (getattr(e, "matrix", None) is not None and not getattr(e, "elided", False)
                and not _too_wide(e)):
            best = max(best, getattr(e, "ncols", 0) or 0)
    return best


def _is_wide(ev):
    """A matrix scaled to the line width with ``\\qlmat``: present, not record-elided, not
    page-elided (>25), and 11..25 columns wide."""
    return (getattr(ev, "matrix", None) is not None and not getattr(ev, "elided", False)
            and not _too_wide(ev)
            and (getattr(ev, "ncols", 0) or 0) > MATRIX_COL_DEFAULT)


def _zero_map_note(ev):
    """A short parenthetical disambiguating a zero-DIMENSIONAL matrix (which _pmatrix
    prints as the symbol ``0``) as the zero map, not the scalar 0. Empty otherwise."""
    if (getattr(ev, "matrix", None) is not None and not getattr(ev, "elided", False)
            and (ev.nrows == 0 or ev.ncols == 0)):
        return r"\quad(\text{the zero map})"
    return ""


def _pmatrix(rs):
    """The matrix as `pmatrix` TeX source (shared by render_latex and render_html), or a
    `\\text{...}` note when the recorder dropped the body (the 250k memory backstop). Works
    for any event carrying (elided, matrix, nrows, ncols, note): RankStep,
    ModuleDifferential, ExtDegree.

    This renders the matrix IN FULL (the HTML report's complete view); the PDF's
    page-bounded scaling (\\qlmat) and >25-row/col elision are layered on top by the
    display helpers (``_matrix_eq_lines`` / ``_eq_lines`` / ``_page_elision_prose``),
    NOT here, so render_html keeps showing every matrix in full.

    A zero-DIMENSIONAL matrix (no rows or no columns -- e.g. every zero Ext/Tor
    differential) renders as the symbol ``0``: an empty ``\\begin{pmatrix}\\end{pmatrix}``
    typesets as a stray ``()`` in the PDF (Plan 34 BLOCKING-2)."""
    if rs.elided or rs.matrix is None:
        return r"\text{%s}" % _tex_escape(rs.note)
    if rs.nrows == 0 or rs.ncols == 0:
        return "0"
    rows = r" \\ ".join(" & ".join(rs.matrix[i][j] for j in range(rs.ncols))
                        for i in range(rs.nrows))
    return r"\begin{pmatrix} %s \end{pmatrix}" % rows


def _page_elision_prose(lhs, ev):
    """A page-elision as WRAPPING prose (not display math -- a long ``\\text{}`` in ``\\[
    \\]`` cannot line-break and overflows the page): ``$lhs$ is an r x c matrix over F,
    shown in full in the HTML/JSON report.`` The full matrix is in the complete report."""
    return (r"\noindent $%s$ is a $%d\times%d$ matrix over %s, shown in full in the "
            r"HTML/JSON report.\\" % (lhs, ev.nrows, ev.ncols, _tex_escape(ev.field)))


def _matrix_eq_lines(decl, lhs, ev):
    """Display line(s) for a differential: its map declaration then ``lhs = matrix``.
    A matrix past 25 rows/cols is page-elided to wrapping prose (pointing at the complete
    report); an 11..25-column one scales the whole equation to the line width on its OWN
    line (so the ``lhs =`` prefix scales with it, never overflowing); a narrow (<=10) or
    record-elided one keeps the compact single-line form."""
    if _too_wide(ev):
        return [r"\[ %s \]" % decl, _page_elision_prose(lhs, ev)]
    body = _pmatrix(ev)
    if _is_wide(ev):
        return [r"\[ %s \]" % decl, r"\[ \qlmat{%s = %s} \]" % (lhs, body)]
    return [r"\[ %s \qquad %s = %s%s \]" % (decl, lhs, body, _zero_map_note(ev))]


def _eq_lines(lhs, ev):
    """Display line(s) for a bare ``lhs = matrix`` (no separate declaration), wide-aware
    exactly like ``_matrix_eq_lines``."""
    if _too_wide(ev):
        return [_page_elision_prose(lhs, ev)]
    body = _pmatrix(ev)
    if _is_wide(ev):
        return [r"\[ \qlmat{%s = %s} \]" % (lhs, body)]
    return [r"\[ %s = %s%s \]" % (lhs, body, _zero_map_note(ev))]


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
            # A long uniserial has many Loewy layers side by side; \qlmat scales the
            # display to the line width so it never runs off the page (shrink-only, so
            # short P_v/I_v stay natural size).
            out.append(r"\[ \qlmat{%s_{%s} = %s \qquad \dim = %d,\ \underline{\dim} = %s} \]"
                       % (sym, v, layers, d["dim"], _dimvec_tex(d["dimvec"])))
    return out


def _dimvec_tex(dv):
    return "(" + ",\\ ".join("%s{:}%s" % (k, val) for k, val in dv.items()) + ")"


def _module_steps_latex(events):
    """The module worked-steps (ModuleTerm / ModuleDifferential / ExtDegree /
    StepNote) as a LaTeX subsection, in emission order; empty when there are none.

    Plan 34 (homework depth): a StepNote with ``heading=True`` opens a numbered worked
    step, rendered as a run-in ``\\paragraph{Step N. ...}`` heading whose narration
    (definition + justification) flows straight after it; plain StepNotes render as
    paragraphs of narration between the displayed matrices. Ext/Tor degrees spell the
    dimension count out in full (``dim = dim Hom - rank d - rank d = ... = result``)."""
    mods = [e for e in events
            if isinstance(e, (ModuleTerm, ModuleDifferential, ExtDegree, StepNote))]
    if not mods:
        return []
    out = [r"\subsection*{Worked module steps}", _tex_escape(_LATEX_MATRIX_POLICY), ""]
    step_no = 0
    for e in mods:
        if isinstance(e, StepNote):
            if e.heading:
                step_no += 1
                head = r"\paragraph{Step %d.~%s}" % (step_no, _tex_escape(e.text))
                out.append(head + (" " + _tex_escape(e.detail) if e.detail else ""))
            else:
                out.append(_tex_escape(e.text))
                if e.detail:
                    out.append(r"\emph{%s}" % _tex_escape(e.detail))
            out.append("")                       # paragraph break
        elif isinstance(e, ModuleTerm):
            what = "Q_{%d}" % e.degree if e.sym == "P" else "E^{%d}" % e.degree
            out.append(r"\[ %s = %s \qquad \dim = %d \]"
                       % (what, oplus_tex(e.summands, e.sym), e.dim))
        elif isinstance(e, ModuleDifferential):
            name = getattr(e, "mod_name", "M") or "M"
            dom_lbl = name if getattr(e, "dom_is_module", False) \
                else oplus_tex(e.dom_summands, e.sym)
            cod = name if e.cod_is_module else oplus_tex(e.cod_summands, e.sym)
            decl = r"%s : %s \to %s" % (e.symbol, dom_lbl, cod)
            out.extend(_matrix_eq_lines(decl, e.symbol, e))
        elif isinstance(e, ExtDegree):
            out.extend(_ext_degree_latex(e))
    from quiverlab.trace.render_text import ext_result_runs
    runs = ext_result_runs(events)
    if runs:
        out.append(r"\subsection*{Result}")
        for op, dims in runs:
            sep = "_" if op == "Tor" else "^"
            out.append(r"\[ %s \]" % r",\; ".join(
                r"\operatorname{%s}%s{%d} = %d" % (op, sep, i, d)
                for i, d in enumerate(dims)))
    return out


def _ext_degree_latex(e):
    """One Ext/Tor degree with the full rank arithmetic spelled out (homework depth):
    the collapsed Hom/tensor dimension, the two neighbouring differential ranks, the
    displayed differential matrix, and ``dim = space - rank - rank = ... = result``."""
    if e.op == "Tor":
        # homology: Tor_n = ker d_n / im d_{n+1}; the SHOWN matrix is d_{n+1} (the map
        # into degree n, rank_here); the outgoing d_n has rank_prev.
        shown, other = "d_{%d}" % (e.degree + 1), "d_{%d}" % e.degree
        space_word = r"\dim(Q_{%d}\otimes_A N)" % e.degree
        quotient = r"\ker d_{%d}/\operatorname{im}d_{%d}" % (e.degree, e.degree + 1)
        arith = "%d - %d - %d" % (e.space_dim, e.rank_prev, e.rank_here)
    else:  # Ext (cohomology): Ext^n = ker delta^n / im delta^{n-1}
        shown, other = r"\delta^{%d}" % e.degree, r"\delta^{%d}" % (e.degree - 1)
        space_word = r"\dim\operatorname{Hom}(Q_{%d},N)" % e.degree
        quotient = r"\ker\delta^{%d}/\operatorname{im}\delta^{%d}" % (e.degree, e.degree - 1)
        arith = "%d - %d - %d" % (e.space_dim, e.rank_here, e.rank_prev)
    label = (r"\operatorname{%s}_{%d}(M,N)" % (e.op, e.degree) if e.op == "Tor"
             else r"\operatorname{%s}^{%d}(M,N)" % (e.op, e.degree))
    lines = [
        r"\noindent $%s$: $%s = %d$, $\operatorname{rank}%s = %d$, "
        r"$\operatorname{rank}%s = %d$."
        % (label, space_word, e.space_dim, shown, e.rank_here, other, e.rank_prev)]
    lines += _eq_lines(shown, e)
    lines.append(r"\[ %s = %s,\qquad \dim = %s = %d \]"
                 % (label, quotient, arith, e.result_dim))
    lines.append("")
    return lines


def render_latex(events, title="", references=(), algebra=None):
    events = list(events)
    from quiverlab.errors import QuiverlabError
    from quiverlab.trace.events import ALL_EVENTS
    for e in events:
        if not isinstance(e, ALL_EVENTS):
            raise QuiverlabError(
                "render_latex received a non-event object of type %r -- likely an "
                "unpacked (events, result) tuple from a trace_* helper; every "
                "element of the stream must be a trace event" % type(e).__name__)
    out = [r"\documentclass{article}", r"\usepackage{amsmath}"]
    out.extend(matrix_preamble_lines(_events_max_cols(events)))
    out += [r"\begin{document}", r"\section*{Worked steps: %s}" % _tex_escape(title)]
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
            if _too_wide(rs):
                out.append(r"\[ \operatorname{rank}(%s) = %d \]" % (sym, rs.rank))
                out.append(_page_elision_prose(sym, rs))
            elif _is_wide(rs):
                out.append(r"\[ \qlmat{%s = %s} \]" % (sym, _pmatrix(rs)))
                out.append(r"\[ \operatorname{rank} = %d \]" % rs.rank)
            else:
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
                 ("^", r"\textasciicircum{}"), ("~", r"\textasciitilde{}"),
                 # '<', '>', '|' render as inverted punctuation / an em-dash in the
                 # default OT1 text font (so "A -> B" would print "A -¿ B"); map them to
                 # the base-LaTeX text commands. Added AFTER the brace steps so their own
                 # braces are not re-escaped (same discipline as ^ and ~).
                 ("<", r"\textless{}"), (">", r"\textgreater{}"),
                 ("|", r"\textbar{}")):
        s = s.replace(a, b)
    return s.replace(sentinel, r"\textbackslash{}")
