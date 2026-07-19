from quiverlab.families.basic import linear_path_algebra, truncated_polynomial  # noqa: F401
from quiverlab.families.dynkin import dynkin_quiver  # noqa: F401
from quiverlab.families.exterior import ExteriorAlgebra  # noqa: F401
from quiverlab.families.incidence import IncidenceAlgebra  # noqa: F401
from quiverlab.families.nakayama import NakayamaAlgebra  # noqa: F401
from quiverlab.families.path_algebra import PathAlgebra  # noqa: F401
from quiverlab.families.preprojective import PreprojectiveAlgebra  # noqa: F401
from quiverlab.families.quantum import QuantumCI  # noqa: F401
from quiverlab.families.radical_square_zero import RadicalSquareZero  # noqa: F401
from quiverlab.families.tensor import TensorProduct  # noqa: F401
from quiverlab.families.trivial_extension import TrivialExtension  # noqa: F401
from quiverlab.families.truncated import TruncatedPathAlgebra  # noqa: F401
from quiverlab.families.discover import (  # noqa: F401
    CATALOG, FamilyInfo, FamilyListing, families,
)


def zoo(dim_max=12):
    """Curated exact zoo of open (Han-conjecture) algebras (spec §3.4).

    PLACEHOLDER: the iterator itself is delivered by Plan 06 Task 11. Task 10
    wires only its catalog entry and public name, so ``families()`` and the
    top-level ``quiverlab`` exports are already complete; Task 11 replaces this
    stub with the real iterator. Loud until then.
    """
    raise NotImplementedError(
        "zoo() is delivered by Plan 06 Task 11; Task 10 wires only its catalog "
        "entry and public export."
    )
