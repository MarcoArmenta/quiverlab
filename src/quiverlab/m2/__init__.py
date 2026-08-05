"""Macaulay2 oracle bridge (subprocess). Importing this package does NOT
probe the M2 binary -- everything is lazy, mirroring quiverlab.qpa.

``crosscheck`` is exposed as a SUBMODULE (not a package-level function): the
typed comparisons and the dispatcher live in ``quiverlab.m2.crosscheck`` and
are reached via ``from quiverlab.m2 import crosscheck`` -> ``crosscheck.crosscheck``
/ ``crosscheck.crosscheck_graded_dims`` / ``crosscheck.crosscheck_commutative_ext``.
The submodule is imported lazily on first ``from ... import crosscheck`` so a bare
``import quiverlab.m2`` never pulls the comparison layer (nor probes M2)."""
from quiverlab.m2.session import m2_available, require_m2, should_skip_m2

__all__ = ["m2_available", "require_m2", "should_skip_m2", "crosscheck"]
