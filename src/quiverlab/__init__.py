"""quiverlab: quivers with relations and Hochschild theory, exactly."""

__version__ = "0.1.0"

# Worked-steps traces are ON by default (spec D9). Flip per-call via
# A.hochschild_cohomology(..., verbose=False) or globally via quiverlab.verbose.
verbose = True

from quiverlab.errors import (  # noqa: E402,F401
    QuiverlabError, ExactnessError, FieldError, RelationError,
    AdmissibilityError, NotFiniteDimensionalError, DepthLimitError,
    QpaUnavailableError,
)
from quiverlab.fields import GF  # noqa: E402,F401
from quiverlab.fields import CC, E  # noqa: E402,F401
from quiverlab.fields import QQi  # noqa: E402,F401
from quiverlab.combinat import Quiver  # noqa: E402,F401
from quiverlab.core import Algebra  # noqa: E402,F401
from quiverlab.families import linear_path_algebra, truncated_polynomial  # noqa: E402,F401
from quiverlab.families import (  # noqa: E402,F401
    NakayamaAlgebra, PathAlgebra, TruncatedPathAlgebra, RadicalSquareZero,
    IncidenceAlgebra, QuantumCI, ExteriorAlgebra, PreprojectiveAlgebra,
    TrivialExtension, TensorProduct, zoo, families,
)
from quiverlab.families import BrauerGraph, BrauerGraphAlgebra  # noqa: E402,F401
from quiverlab.citations import bibliography  # noqa: E402,F401
from quiverlab.invariants.sweep import sweep  # noqa: E402,F401
from quiverlab.modules.complexes import ChainComplex, ChainMap  # noqa: E402,F401

__all__ = [
    "__version__",
    "verbose",
    "QuiverlabError", "ExactnessError", "FieldError", "RelationError",
    "AdmissibilityError", "NotFiniteDimensionalError", "DepthLimitError",
    "QpaUnavailableError",
    "GF", "CC", "E", "QQi",
    "Quiver", "Algebra",
    "truncated_polynomial", "linear_path_algebra",
    "NakayamaAlgebra", "PathAlgebra", "TruncatedPathAlgebra", "RadicalSquareZero",
    "IncidenceAlgebra", "QuantumCI", "ExteriorAlgebra", "PreprojectiveAlgebra",
    "TrivialExtension", "TensorProduct", "zoo", "families",
    "BrauerGraph", "BrauerGraphAlgebra",
    "bibliography",
    "sweep",
    "ChainComplex", "ChainMap",
]
