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
    AmbiguityEvent, ReductionStep,
)
from quiverlab.trace.recorder import TERMS_ELISION


def derive_dims(events):
    """HH^n / HH_n = collapsed_dim(n) - rank(n) - rank(n-1), from the events."""
    cn = {e.degree: e.collapsed_dim for e in events if isinstance(e, ResolutionTerm)}
    rk = {e.degree: e.rank for e in events if isinstance(e, RankStep)}
    dims = []
    for n in sorted(cn):
        dims.append(cn[n] - rk.get(n, 0) - rk.get(n - 1, 0))
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


def render_text(events, title="", references=()):
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
