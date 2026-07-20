"""quiverlab.resolutions_cs: the Chouhy-Solotar general bimodule resolution
(arXiv:1406.2300 = J. Algebra 432 (2015)). Domain-generic; consumes the Plan-03
groebner.ReductionSystem and the Plan-01 core.Algebra. See the plan and
docs/internals/09-chouhy-solotar.md for the mathematics (S-sequence = Bardzell
associated paths of the tip monomial algebra; differentials = CS f_n leading map +
order-condition-pinned correction; two collapse maps for HH_• and HH^•).

Public surface (the frozen interface Plan 05 and later consume from this package;
see the Plan-04 T14 frozen-interface note). Re-exported here so the package root
exposes the documented API — the guards/errors themselves live canonically in
``quiverlab.errors``."""

from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
from quiverlab.resolutions_cs.engine_facade import CSResolution
from quiverlab.resolutions_cs.homology import (
    cs_cohomology_dims,
    cs_homology_dims,
    cs_hh_basis,
)
from quiverlab.resolutions_cs.build import reduction_system_of
from quiverlab.resolutions_cs.ambiguities import SSequence
from quiverlab.resolutions_cs.comparison import Comparison

__all__ = [
    "ChouhySolotarResolution",
    "CSResolution",
    "cs_cohomology_dims",
    "cs_homology_dims",
    "cs_hh_basis",
    "reduction_system_of",
    "SSequence",
    "Comparison",
]
