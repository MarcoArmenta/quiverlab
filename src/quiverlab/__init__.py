"""quiverlab: quivers with relations and Hochschild theory, exactly."""

__version__ = "0.1.0.dev0"

from quiverlab.errors import (  # noqa: E402,F401
    QuiverlabError, ExactnessError, FieldError, RelationError,
    AdmissibilityError, NotFiniteDimensionalError, DepthLimitError,
)
from quiverlab.fields import GF  # noqa: E402,F401
from quiverlab.fields import CC, E  # noqa: E402,F401
