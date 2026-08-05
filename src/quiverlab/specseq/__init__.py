"""Spectral sequences over ``fields.linalg`` (Plan 42).

Bounded filtered complexes and double complexes, exact Weibel-section-5.4 pages
with byte-reproducible representatives, and a standing convergence certificate
(``E_inf`` totals == total-complex homology). The public surface + the four
presets are re-exported at the end of the plan (Task 10)."""
from quiverlab.specseq.filtered import FilteredComplex

__all__ = ["FilteredComplex"]
