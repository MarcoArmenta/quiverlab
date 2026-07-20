"""quiverlab.trace: the worked-steps trace subsystem (spec §3.8, component 11).

Typed step events (events.py), a recording buffer with size-eliding rules
(recorder.py), engine->citation provenance (provenance.py), and three renderers
(render_text/render_latex/render_html) driven by writer.py. Default-on per D9
(quiverlab.verbose = True)."""
from quiverlab.trace.events import (  # noqa: F401
    Dispatch, ReductionStep, AmbiguityEvent, ResolutionTerm,
    DifferentialEvent, LiftStep, RankStep,
)

__all__ = [
    "Dispatch", "ReductionStep", "AmbiguityEvent", "ResolutionTerm",
    "DifferentialEvent", "LiftStep", "RankStep",
]
