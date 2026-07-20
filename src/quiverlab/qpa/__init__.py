"""quiverlab.qpa: optional GAP/QPA cross-check backend (spec §5 component 12).
Imported lazily by Algebra.crosscheck; importing this module does NOT import GAP
(that happens on first gap_available()/crosscheck call)."""
from quiverlab.qpa.session import gap_available, require_gap  # noqa: F401
from quiverlab.qpa.crosscheck import crosscheck, CrosscheckReport  # noqa: F401

__all__ = ["gap_available", "require_gap", "crosscheck", "CrosscheckReport"]
