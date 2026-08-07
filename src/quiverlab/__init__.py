"""quiverlab: quivers with relations and Hochschild theory, exactly."""

__version__ = "0.3.0"

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
    OnePointExtension, repetitive_slice, JacobianAlgebra, Potential, cyclic_derivative,
)
from quiverlab.families import BrauerGraph, BrauerGraphAlgebra  # noqa: E402,F401
from quiverlab.citations import bibliography  # noqa: E402,F401
from quiverlab.invariants.sweep import sweep  # noqa: E402,F401
from quiverlab.modules.complexes import ChainComplex, ChainMap  # noqa: E402,F401
from quiverlab.specseq import (  # noqa: E402,F401
    DoubleComplex, FilteredComplex, SpectralSequence)
from quiverlab.derived import (  # noqa: E402,F401
    hyper_hom_basis, tau_Db, tau_Db_minus, is_tilting_complex,
    end_algebra_of_complex, two_term_silting_from_presentation,
    derived_fingerprint, compare_fingerprints,
)
from quiverlab.surfaces import (  # noqa: E402,F401
    MarkedSurface, Triangulation, fan_triangulation, annulus_triangulation,
    once_punctured_torus, hexagon_with_internal_triangle,
    quiver_of, potential_of, jacobian_of, flip, certify_flip_mutation, surface_block,
)

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
    "OnePointExtension", "repetitive_slice", "JacobianAlgebra", "Potential",
    "cyclic_derivative",
    "bibliography",
    "sweep",
    "ChainComplex", "ChainMap",
    "FilteredComplex", "DoubleComplex", "SpectralSequence",
    # Plan 43 -- derived-category surface
    "hyper_hom_basis", "tau_Db", "tau_Db_minus", "is_tilting_complex",
    "end_algebra_of_complex", "two_term_silting_from_presentation",
    "derived_fingerprint", "compare_fingerprints",
    # Plan 48 -- marked surfaces -> triangulations -> gentle Jacobian algebras
    "MarkedSurface", "Triangulation", "fan_triangulation", "annulus_triangulation",
    "once_punctured_torus", "hexagon_with_internal_triangle",
    "quiver_of", "potential_of", "jacobian_of", "flip", "certify_flip_mutation",
    "surface_block",
]
