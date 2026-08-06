"""Spectral sequences over ``fields.linalg`` (Plan 42).

Bounded filtered complexes and double complexes, exact Weibel-section-5.4 pages
with byte-reproducible representatives, and a standing convergence certificate
(``E_inf`` totals == total-complex homology), plus the four presets: the Hochschild
``(b, B)`` bicomplex, the radical / associated-graded filtration, and the
Grothendieck / Cartan-Eilenberg change-of-rings sequence."""
from quiverlab.specseq.convergence import ConvergenceReport
from quiverlab.specseq.double import DoubleComplex
from quiverlab.specseq.filtered import FilteredComplex
from quiverlab.specseq.pages import Page, SpectralSequence, Subquotient
from quiverlab.specseq.presets import (
    cartan_eilenberg_ss,
    grothendieck_double_complex,
    hochschild_bB_ss,
    radical_filtration_ss,
)

__all__ = [
    "FilteredComplex", "DoubleComplex", "SpectralSequence", "Page",
    "Subquotient", "ConvergenceReport",
    "hochschild_bB_ss", "radical_filtration_ss", "cartan_eilenberg_ss",
    "grothendieck_double_complex",
]
