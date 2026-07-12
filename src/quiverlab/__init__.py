"""quiverlab: quivers with relations and Hochschild theory, exactly."""

__version__ = "0.1.0.dev0"

from quiverlab.errors import (  # noqa: E402,F401
    QuiverlabError, ExactnessError, FieldError, RelationError,
    AdmissibilityError, NotFiniteDimensionalError, DepthLimitError,
)
from quiverlab.fields import GF  # noqa: E402,F401
from quiverlab.fields import CC, E  # noqa: E402,F401
from quiverlab.combinat import Quiver  # noqa: E402,F401
from quiverlab.core import Algebra  # noqa: E402,F401
from quiverlab.families import linear_path_algebra, truncated_polynomial  # noqa: E402,F401

__all__ = [
    "__version__",
    "QuiverlabError", "ExactnessError", "FieldError", "RelationError",
    "AdmissibilityError", "NotFiniteDimensionalError", "DepthLimitError",
    "GF", "CC", "E",
    "Quiver", "Algebra",
    "truncated_polynomial", "linear_path_algebra",
]
