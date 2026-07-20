"""Trace recorder: a bounded, list-compatible event buffer with size-eliding rules
(spec §3.8 performance guard). verbose must NOT blow up a long computation:

  * the buffer is capped at MAX_EVENTS (5000); beyond that, events are dropped and
    counted (one elision note), so memory is O(MAX_EVENTS) regardless of depth;
  * a differential matrix larger than MATRIX_ELISION_CELLS (400 = 20x20) is not
    stored -- only its (shape, rank) survive, with an explicit elision note (so a
    deep resolution records at most one small ResolutionTerm + one RankStep per
    degree, i.e. O(top) small records, never a giant matrix dump);
  * DifferentialEvent term lists over TERMS_ELISION (100) are the renderers'
    responsibility to truncate (Task 9/10) -- the recorder keeps the cap constant
    available to them.

Concretely: a top=40 monomial resolution records ~41 ResolutionTerm + ~41 RankStep
(each RankStep either a <=400-cell matrix or a one-line elision note) plus the
capped construction ReductionSteps -- well under MAX_EVENTS, bounded memory."""
from quiverlab.trace.events import RankStep

MAX_EVENTS = 5000
MATRIX_ELISION_CELLS = 400
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
            self.elided_events += 1
            return
        self.events.append(event)

    append = record  # list-compatible

    def __iter__(self):
        return iter(self.events)

    def __len__(self):
        return len(self.events)


def rankstep(degree, side, D, nrows, ncols, rank, dom):
    """Build a RankStep from a domain-element matrix D (list of lists, row-major
    D[row][col]). Param order (nrows, ncols) matches the RankStep field order and
    D's indexing. Elide the matrix body (keep shape + rank) when it exceeds
    MATRIX_ELISION_CELLS."""
    cells = nrows * ncols
    if cells > MATRIX_ELISION_CELLS:
        return RankStep(
            degree=degree, side=side, nrows=nrows, ncols=ncols, rank=rank,
            field=dom.name, matrix=None, elided=True,
            note=(f"{nrows}x{ncols} matrix over {dom.name} elided "
                  f"(> {MATRIX_ELISION_CELLS} cells); rank {rank} recorded"))
    rendered = [[str(D[i][j]) for j in range(ncols)] for i in range(nrows)]
    return RankStep(
        degree=degree, side=side, nrows=nrows, ncols=ncols, rank=rank,
        field=dom.name, matrix=rendered, elided=False, note="")


def resolve_verbose(per_call, global_flag):
    """Per-call verbose (True/False) overrides the global flag; None defers to it."""
    return bool(global_flag) if per_call is None else bool(per_call)
