"""quiverlab.groebner: noncommutative Groebner completion for path-algebra ideals
(spec §5, component 3). Emits the ReductionSystem consumed by Chouhy-Solotar
(Plan 04). Path composition is LEFT TO RIGHT throughout."""
from quiverlab.groebner.order import PathOrder, path_order  # noqa: F401
from quiverlab.groebner.events import Dispatch, ReductionStep  # noqa: F401
from quiverlab.groebner.reduction import ReductionRule  # noqa: F401
from quiverlab.groebner.overlap import Ambiguity  # noqa: F401
from quiverlab.groebner.system import ReductionSystem, build_reduction_system  # noqa: F401

__all__ = [
    "PathOrder", "path_order", "Dispatch", "ReductionStep",
    "ReductionRule", "Ambiguity", "ReductionSystem", "build_reduction_system",
]
