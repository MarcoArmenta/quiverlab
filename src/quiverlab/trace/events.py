"""The unified worked-steps step-event taxonomy (spec §3.8).

This module is the SINGLE import surface for trace events. It RE-EXPORTS the
event dataclasses shipped by Plan 03 (groebner) and Plan 04 (Chouhy-Solotar)
from their home modules -- it never redefines them -- and DEFINES the one new
event, RankStep, which carries a numeric differential matrix over a field plus
its rank (the bar/fast engines' notion of a differential; the CS
DifferentialEvent carries symbolic bimodule terms instead).

Reconciliation notes (see the plan's taxonomy table):
  * The requested `Differential` is Plan 04's DifferentialEvent.
  * AmbiguityEvent (Plan 04) is folded in.
  * Dispatch (Plan 03) is reused for BOTH the construction route (monomial vs
    groebner) AND the engine/resolution choice (bar/fast/bardzell/chouhy-solotar);
    same dataclass, distinguished by the value of `route`.
"""
from dataclasses import dataclass

from quiverlab.groebner.events import Dispatch, ReductionStep  # noqa: F401
from quiverlab.resolutions_cs.trace import (  # noqa: F401
    AmbiguityEvent, ResolutionTerm, DifferentialEvent, LiftStep,
)


@dataclass
class RankStep:
    """One rank computation of a differential matrix over a stated field.

    `matrix` is a list[list[str]] of domain-element string renderings (kept small),
    or None when `elided` is True (matrix larger than the elision threshold -- only
    shape + rank are retained). `side` is "cochain" (d^n : C^n -> C^{n+1}) or
    "chain" (b_n : C_n -> C_{n-1}). `nrows`/`ncols` are the matrix dimensions
    (= dim of the target/source cochain space)."""
    degree: int
    side: str
    nrows: int
    ncols: int
    rank: int
    field: str
    matrix: object = None
    elided: bool = False
    note: str = ""


__all__ = [
    "Dispatch", "ReductionStep", "AmbiguityEvent", "ResolutionTerm",
    "DifferentialEvent", "LiftStep", "RankStep",
]
