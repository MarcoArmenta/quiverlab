"""Plain-text worked-steps renderer (spec §3.8). The resulting dimensions are
DERIVED from the recorded ResolutionTerm + RankStep events (never echoed) -- that
is the binding discipline: the golden tests assert these derived numbers equal the
engine's own .dims, so a trace can never claim something the engine did not compute.

The kind-detection (`_dims_kind`) and matrix-format (`_matrix_block`) helpers are
written to be reusable by the LaTeX/HTML renderers (Tasks 10/11) -- they take a
RankStep / the event list and return plain data, with no text-format assumptions
baked in beyond the ASCII bracket layout pinned by the golden."""
from quiverlab.trace.events import (
    Dispatch, ResolutionTerm, RankStep, DifferentialEvent, LiftStep,
    AmbiguityEvent, ReductionStep, ModuleTerm, ModuleDifferential, ExtDegree,
    StepNote,
)
from quiverlab.trace.recorder import TERMS_ELISION, MATRIX_ELISION_CELLS

# The single sentence stating the matrix-elision threshold, emitted once per report
# (spec: "the threshold stated once in the preamble") whenever module matrices are
# present; per-matrix notes then stay terse.
ELISION_PREAMBLE = (
    "Matrices with more than %d entries are shown as [shape r x c, rank k -- "
    "elided]; every step still appears." % MATRIX_ELISION_CELLS)

# The direct-sum / composition-factor formatting shared by all three renderers.
# `oplus_text` groups a vertex list with repetition into `P_1^2 (+) P_3`; the
# LaTeX/HTML renderers reuse `oplus_tex` (defined in render_latex).


def oplus_text(summands, sym):
    """A direct-sum of indecomposables ``P_1^2 (+) P_3`` from a vertex list with
    repetition; ``0`` for an empty term. ASCII (plain-text renderer)."""
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
    return " (+) ".join("%s_%s" % (sym, v) if c == 1 else "%s_%s^%d" % (sym, v, c)
                        for v, c in groups)


def factor_stack_text(dimvec):
    """A semisimple composition-factor layer ``S_1 + S_2^2`` from a str-keyed
    multiplicity dict; ``0`` when empty."""
    parts = ["S_%s" % v if m == 1 else "S_%s^%d" % (v, m)
             for v, m in dimvec.items() if m > 0]
    return " + ".join(parts) if parts else "0"


def derive_dims(events):
    """Resulting (co)homology dims, DERIVED from the recorded events, side-aware:
      * cohomology (cochain side): HH^n = C_n - rank(d^n) - rank(d^{n-1});
      * homology   (chain side):   HH_n = C_n - rank(b_n) - rank(b_{n+1}).
    The engine emits RankStep(degree=n) <-> b_n for homology, so the neighbouring
    rank is rk[n+1], NOT rk[n-1]. Which pairing to use is read from the RankStep
    events' `.side` ("cochain" vs "chain"); `.get(...,0)` treats a missing degree
    as rank 0 (e.g. d^{-1} = 0, or b_0 = 0)."""
    cn = {e.degree: e.collapsed_dim for e in events if isinstance(e, ResolutionTerm)}
    rk = {e.degree: e.rank for e in events if isinstance(e, RankStep)}
    chain = any(getattr(e, "side", "") == "chain"
                for e in events if isinstance(e, RankStep))
    dims = []
    for n in sorted(cn):
        other = rk.get(n + 1, 0) if chain else rk.get(n - 1, 0)
        dims.append(cn[n] - rk.get(n, 0) - other)
    return dims


def _dims_kind(events):
    """"HH^" (cohomology) if any recorded differential is a cochain, else "HH_"
    (homology). Shared kind-detection for all three renderers."""
    return ("HH^" if any(getattr(e, "side", "") == "cochain"
                         for e in events if isinstance(e, RankStep)) else "HH_")


def _matrix_block(rs):
    """The ASCII lines for a RankStep's matrix (right-justified columns), or the
    one-line elision note when the body was dropped. Reused by the LaTeX/HTML
    renderers as the fallback plain rendering."""
    if rs.elided or rs.matrix is None:
        return ["  (%s)" % rs.note]
    widths = [max(len(rs.matrix[i][j]) for i in range(rs.nrows)) for j in range(rs.ncols)]
    out = []
    for i in range(rs.nrows):
        row = "  ".join(rs.matrix[i][j].rjust(widths[j]) for j in range(rs.ncols))
        out.append("    [ %s ]" % row)
    return out


def compute_algebra_objects(algebra):
    """``trace.modules.algebra_objects(algebra)`` guarded: a description of the
    projectives/injectives of A, or ``None`` (the report degrades to no P/I section)
    with a reason string. Descriptive-only, so a hiccup must not sink the bundle."""
    if algebra is None:
        return None, None
    try:
        from quiverlab.trace.modules import algebra_objects
        return algebra_objects(algebra), None
    except Exception as exc:                                   # pragma: no cover
        return None, "%s: %s" % (type(exc).__name__, exc)


def pi_section_lines(objects, note=None):
    """The "projectives and injectives of A" section as plain-text lines; empty when
    there is nothing to show."""
    if not objects:
        return ["The projectives and injectives of A",
                "  (unavailable: %s)" % note] if note else []
    lines = ["The projectives and injectives of A",
             "  (simples S_v omitted; Loewy layers listed top to bottom)"]
    for row in objects:
        v = row["vertex"]
        for sym in ("P", "I"):
            d = row[sym]
            stack = " / ".join(factor_stack_text(L) for L in d["layers"]) or "0"
            lines.append("  %s_%s: dim %d, dimvec %s"
                         % (sym, v, d["dim"], _fmt_dimvec(d["dimvec"])))
            lines.append("      Loewy layers: %s" % stack)
            lines.append("      top = %s, soc = %s"
                         % (factor_stack_text(d["top"]),
                            factor_stack_text(d["socle"])))
    lines.append("")
    return lines


def _fmt_dimvec(dv):
    return "(" + ", ".join("%s:%s" % (k, v) for k, v in dv.items()) + ")"


def ext_result_dims(events):
    """The Ext/Tor result dims, in degree order, from the ExtDegree events."""
    exts = sorted((e for e in events if isinstance(e, ExtDegree)), key=lambda e: e.degree)
    return [e.result_dim for e in exts], (exts[0].op if exts else None)


def _module_steps_lines(events):
    """The module worked-steps (ModuleTerm / ModuleDifferential / ExtDegree /
    StepNote) rendered in emission order; empty when there are none."""
    mods = [e for e in events
            if isinstance(e, (ModuleTerm, ModuleDifferential, ExtDegree, StepNote))]
    if not mods:
        return []
    lines = ["Worked module steps", "  " + ELISION_PREAMBLE, ""]
    for e in mods:
        if isinstance(e, StepNote):
            lines.append(e.text)
            if e.detail:
                lines.append("  " + e.detail)
        elif isinstance(e, ModuleTerm):
            what = "term Q_%d" % e.degree if e.sym == "P" else "term E^%d" % e.degree
            lines.append("%s = %s (dim %d, dimvec %s)"
                         % (what, oplus_text(e.summands, e.sym), e.dim,
                            _fmt_dimvec(e.dimvec or {})))
        elif isinstance(e, ModuleDifferential):
            cod = "M" if e.cod_is_module else oplus_text(e.cod_summands, e.sym)
            lines.append("  %s : %s -> %s, %d x %d over %s:"
                         % (e.symbol, oplus_text(e.dom_summands, e.sym), cod,
                            e.nrows, e.ncols, e.field))
            lines.extend(_matrix_block(e))
        elif isinstance(e, ExtDegree):
            lines.append("%s^%d: dim Hom = %d, rank = %d (prev %d) -> %s^%d = %d"
                         % (e.op, e.degree, e.space_dim, e.rank_here, e.rank_prev,
                            e.op, e.degree, e.result_dim))
            lines.append("  connecting map, %d x %d over %s:"
                         % (e.nrows, e.ncols, e.field))
            lines.extend(_matrix_block(e))
        lines.append("")
    dims, op = ext_result_dims(events)
    if dims:
        lines.append("Result: " + "   ".join("%s^%d = %d" % (op, i, d)
                                              for i, d in enumerate(dims)))
        lines.append("")
    return lines


def render_text(events, title="", references=(), algebra=None):
    events = list(events)
    lines = []
    if title:
        lines.append("Worked steps: " + title)
        lines.append("")
    for e in events:
        if isinstance(e, Dispatch):
            lines.append("Chosen resolution: " + e.route)
            lines.append("  reason: " + e.reason)
            lines.append("  defining relations: %d" % e.n_relations)
            lines.append("")
    objects, note = compute_algebra_objects(algebra)
    lines.extend(pi_section_lines(objects, note))
    terms = {e.degree: e for e in events if isinstance(e, ResolutionTerm)}
    ranks = {e.degree: e for e in events if isinstance(e, RankStep)}
    # Drive per-degree rendering over the union of recorded terms and ranks so a
    # RankStep is never dropped for lack of a co-recorded ResolutionTerm (e.g. an
    # isolated elided differential); the "term with ... generators" header is
    # only emitted when the ResolutionTerm was actually recorded.
    for n in sorted(set(terms) | set(ranks)):
        t = terms.get(n)
        if t is not None:
            lines.append("Degree %d: term with %d generators (dim C = %d)"
                         % (n, t.n_generators, t.collapsed_dim))
        if n in ranks:
            rs = ranks[n]
            side = "d^%d" % n if rs.side == "cochain" else "b_%d" % n
            lines.append("  differential %s (%s), %d x %d over %s:"
                         % (side, rs.side, rs.nrows, rs.ncols, rs.field))
            lines.extend(_matrix_block(rs))
            lines.append("  rank = %d" % rs.rank)
        lines.append("")
    # CS symbolic differentials / lifts / ambiguities (present in CS traces)
    for e in events:
        if isinstance(e, DifferentialEvent):
            lines.append("Symbolic differential (degree %d): %d term(s)%s"
                         % (e.degree, len(e.terms),
                            "" if len(e.terms) <= TERMS_ELISION
                            else " (%d shown)" % TERMS_ELISION))
        elif isinstance(e, LiftStep):
            lines.append("Lift/comparison step (degree %d): %s" % (e.degree, e.kind))
        elif isinstance(e, AmbiguityEvent):
            lines.append("Ambiguity chain (degree %d): %d word(s)"
                         % (e.degree, len(e.chain_words)))
    lines.extend(_module_steps_lines(events))
    dims = derive_dims(events)
    if dims:
        kind = _dims_kind(events)
        cells = "   ".join("%s%d = %d" % (kind, i, d) for i, d in enumerate(dims))
        lines.append("Result: " + cells)
        lines.append("")
    if references:
        lines.append("References:")
        for key, entry in references:
            lines.append("  [%s] %s" % (key, entry))
    return "\n".join(lines).rstrip("\n") + "\n"
