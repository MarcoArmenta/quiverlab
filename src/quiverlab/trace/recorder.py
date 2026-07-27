"""Trace recorder: a bounded, list-compatible event buffer (spec §3.8 performance
guard). verbose must NOT blow up a long computation:

  * the buffer is capped at MAX_EVENTS (5000); beyond that, events are dropped and
    counted (one elision note), so memory is O(MAX_EVENTS) regardless of depth;
  * events carry the FULL differential matrix (Plan 34 artifact contract: the events
    are the complete record -- the HTML and trace.json reports show every matrix in
    full). Recording drops a matrix body only past MATRIX_ELISION_CELLS (250_000), a
    pure MEMORY backstop set orders of magnitude above the old 400-cell rule so a
    pathological deep run cannot exhaust memory; past it only (shape, rank) survive
    with an explicit elision note;
  * DifferentialEvent term lists over TERMS_ELISION (100) are the renderers'
    responsibility to truncate (Task 9/10) -- the recorder keeps the cap constant
    available to them.

Concretely: a top=40 monomial resolution records ~41 ResolutionTerm + ~41 RankStep
(each a full matrix under the 250k backstop) plus the capped construction
ReductionSteps -- well under MAX_EVENTS, bounded memory."""
from quiverlab.trace.events import RankStep, ModuleDifferential, ExtDegree, StepNote

MAX_EVENTS = 5000
# The event-stream note emitted when the buffer cap is hit (Plan 34 MINOR): a renderer
# only ever receives ``list(trace)``, so surfacing the elision purely in
# ``Trace.elision_notes`` would let the "every step appears" claim be silently violated.
# The recorder therefore also records ONE StepNote past the cap carrying this text; the
# renderers show it (they already render StepNotes) so the elision is never silent.
BUFFER_FULL_PREFIX = "Trace buffer full:"
# MEMORY backstop only (Plan 34): matrices up to this many cells are recorded IN FULL so
# the events are the complete record; past it, recording keeps just (shape, rank) so a
# pathological deep run stays bounded. The HTML and JSON reports show every recorded
# matrix in full.
MATRIX_ELISION_CELLS = 250_000
TERMS_ELISION = 100


class Trace:
    """A list-like event sink. `.append` is an alias for `.record` so this drops
    in wherever a plain list `trace` is expected (e.g. groebner.reduce_comb)."""

    def __init__(self, max_events=MAX_EVENTS):
        self.events = []
        self.max_events = max_events
        self.elided_events = 0
        self.elision_notes = []

    def record(self, event):
        if len(self.events) >= self.max_events:
            if self.elided_events == 0:
                self.elision_notes.append(
                    f"event buffer full at {self.max_events}; further steps elided")
                # Surface the elision IN the event stream (one reserved slot past the
                # cap) so a renderer -- which only receives list(self) -- shows it and
                # the "every step appears" claim is never silently violated.
                self.events.append(StepNote(
                    "%s recorded the first %d worked steps; later steps were elided "
                    "(the computed dimensions and Result are unaffected)."
                    % (BUFFER_FULL_PREFIX, self.max_events)))
            self.elided_events += 1
            return
        self.events.append(event)

    append = record  # list-compatible

    def __iter__(self):
        return iter(self.events)

    def __len__(self):
        return len(self.events)


def _cell_str(x, dom):
    """A single matrix entry rendered to a READABLE, exact string (Plan 34 MINOR).

    An int / fraction element keeps its canonical ``str`` (``"2"``, ``"-1/2"``) --
    byte-identical to the old ``str(x)`` for every GF(p) / CC-rational trace. An element
    whose ``str`` is NOT int/fraction-parseable (a GF(p^n) extension element, stored
    internally as a little-endian coefficient TUPLE like ``(0, 1)``) is rendered through
    the domain's own human ``to_str`` (``"x^1"``), so the trace never leaks the raw
    internal tuple the way plain ``str(x)`` did. Mirrors ``modules.qpa_module._json_entry``'s
    grammar, kept self-contained here to avoid a modules<-recorder import at load time."""
    s = str(x)
    try:
        int(s)
        return s
    except ValueError:
        pass
    try:
        from fractions import Fraction
        Fraction(s)
        return s
    except (ValueError, ZeroDivisionError):
        to_str = getattr(dom, "to_str", None)
        if callable(to_str):
            try:
                return to_str(x)
            except Exception:                      # pragma: no cover - defensive
                return s
        return s


def rankstep(degree, side, D, nrows, ncols, rank, dom):
    """Build a RankStep from a domain-element matrix D (list of lists, row-major
    `D[row][col]`). Param order (nrows, ncols) matches the RankStep field order and
    D's indexing. The full matrix is recorded; only past the MATRIX_ELISION_CELLS
    MEMORY backstop is the body dropped (shape + rank kept). Report-side readability
    scaling/elision of wide matrices is the LaTeX renderer's job, not the recorder's."""
    if nrows * ncols > MATRIX_ELISION_CELLS:
        return RankStep(
            degree=degree, side=side, nrows=nrows, ncols=ncols, rank=rank,
            field=dom.name, matrix=None, elided=True,
            note=(f"{nrows}x{ncols} matrix over {dom.name} elided "
                  f"(> {MATRIX_ELISION_CELLS} cells); rank {rank} recorded"))
    rendered = [[_cell_str(D[i][j], dom) for j in range(ncols)] for i in range(nrows)]
    return RankStep(
        degree=degree, side=side, nrows=nrows, ncols=ncols, rank=rank,
        field=dom.name, matrix=rendered, elided=False, note="")


def _render_or_elide(D, nrows, ncols, dom, rank=None):
    """(matrix_or_None, elided, note) for a domain-element matrix D, applying the same
    MATRIX_ELISION_CELLS MEMORY backstop as `rankstep`. Beyond the backstop the body is
    dropped and only the shape (+ rank, when supplied) is kept, in the spec's
    ``[shape r x c, rank k -- elided]`` form. (Wide-but-recordable matrices are kept in
    full here; the LaTeX renderer scales/elides them for the page.)"""
    if nrows * ncols > MATRIX_ELISION_CELLS:
        shape = "shape %dx%d" % (nrows, ncols)
        note = ("[%s, rank %d -- elided]" % (shape, rank) if rank is not None
                else "[%s -- elided]" % shape)
        return None, True, note
    rendered = [[_cell_str(D[i][j], dom) for j in range(ncols)] for i in range(nrows)]
    return rendered, False, ""


def module_differential(degree, kind, sym, symbol, dom_summands, cod_summands,
                        D, nrows, ncols, dom, cod_is_module=False, rank=None,
                        dom_is_module=False, mod_name="M"):
    """Build a ModuleDifferential from a domain-element matrix D (row-major
    `D[row][col]`), eliding the body past MATRIX_ELISION_CELLS. `dom_is_module` /
    `cod_is_module` render that endpoint as the traced module itself (self-maps),
    labeled by its actual name via ``mod_name``."""
    matrix, elided, note = _render_or_elide(D, nrows, ncols, dom, rank=rank)
    return ModuleDifferential(
        degree=degree, kind=kind, sym=sym, symbol=symbol,
        dom_summands=list(dom_summands), cod_summands=list(cod_summands),
        nrows=nrows, ncols=ncols, field=dom.name, cod_is_module=cod_is_module,
        dom_is_module=dom_is_module, mod_name=mod_name, matrix=matrix,
        elided=elided, note=note)


def ext_degree(degree, op, space_dim, rank_here, rank_prev, result_dim,
               D, nrows, ncols, dom):
    """Build an ExtDegree from the connecting-map matrix D, eliding the body past
    MATRIX_ELISION_CELLS (the rank is always retained in the note)."""
    matrix, elided, note = _render_or_elide(D, nrows, ncols, dom, rank=rank_here)
    return ExtDegree(
        degree=degree, op=op, space_dim=space_dim, rank_here=rank_here,
        rank_prev=rank_prev, result_dim=result_dim, nrows=nrows, ncols=ncols,
        field=dom.name, matrix=matrix, elided=elided, note=note)


def resolve_verbose(per_call, global_flag):
    """Per-call verbose (True/False) overrides the global flag; None defers to it."""
    return bool(global_flag) if per_call is None else bool(per_call)
