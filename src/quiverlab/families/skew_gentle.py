"""Skew-gentle algebras as a first-class family (Plan 68 / R32).

A thin re-export of the :mod:`quiverlab.skewgentle` package so the split constructor is
discoverable through the family catalog (``families()``) and the top-level surface,
mirroring the ``BrauerGraphAlgebra`` non-scalar-constructor precedent.  The mathematics
lives in ``quiverlab.skewgentle`` (triple / split / modules / certificate / block)."""
from quiverlab.skewgentle import (SkewGentleAlgebra, SkewGentleTriple,  # noqa: F401
                                  is_skew_gentle_triple)

__all__ = ["SkewGentleAlgebra", "SkewGentleTriple", "is_skew_gentle_triple"]
