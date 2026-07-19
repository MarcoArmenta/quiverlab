"""IncidenceAlgebra of a poset. Fixture N5 (diamond == commutative square)."""
import pytest

from quiverlab.errors import RelationError
from quiverlab.families import IncidenceAlgebra
from quiverlab.families.poset import Poset
from quiverlab.fields import CC, GF


DIAMOND = [("b", "x"), ("b", "y"), ("x", "t"), ("y", "t")]


def test_diamond_is_commutative_square_dim9():
    A = IncidenceAlgebra(DIAMOND, field=CC)
    assert A.dim == 9                                   # 4 trivial + 4 arrows + 1 [b,t]
    assert A.hochschild_cohomology(2).dims == [1, 0, 0]  # == Plan-03 commutative square


def test_diamond_char_independent():
    assert IncidenceAlgebra(DIAMOND, field=GF(32003)).dim == 9


def test_chain_poset_is_linear_path_algebra():
    A = IncidenceAlgebra([(1, 2), (2, 3)], field=CC)    # 1 < 2 < 3 chain
    assert A.dim == 6                                   # kA_3, dim n(n+1)/2


def test_directed_cycle_is_not_a_poset():
    with pytest.raises(RelationError):
        Poset([(1, 2), (2, 3), (3, 1)])                 # antisymmetry fails


def test_citations():
    assert "incidence" in IncidenceAlgebra(DIAMOND).citations()


def test_self_cover_rejected_at_construction():
    # A degenerate self-cover (x, x) would build a loop arrow and hang _all_paths;
    # the guard must fire loudly at Poset construction, never reaching that path.
    with pytest.raises(RelationError, match="1"):
        Poset([(1, 2), (1, 1)])
    with pytest.raises(RelationError, match="1"):
        IncidenceAlgebra([(1, 2), (1, 1)])
