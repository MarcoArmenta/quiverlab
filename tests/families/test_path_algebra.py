"""PathAlgebra (hereditary Dynkin/Euclidean). Fixture N2."""
import pytest

from quiverlab.errors import NotFiniteDimensionalError
from quiverlab.families import PathAlgebra
from quiverlab.fields import CC, GF


def test_D4_linear_orientation_dim9_hh_100():
    A = PathAlgebra("D4", field=CC)                 # arrows 1->2, 2->3, 2->4
    assert A.dim == 9                               # Fixture N2 (not 12)
    assert A.hochschild_cohomology(2).dims == [1, 0, 0]


def test_A5_is_triangular_number():
    assert PathAlgebra("A5", field=CC).dim == 15    # 5*6/2


def test_tree_hh_is_char_independent():
    assert PathAlgebra("D4", field=GF(7)).hochschild_cohomology(2).dims == [1, 0, 0]


def test_cyclic_euclidean_orientation_is_loud():
    with pytest.raises(NotFiniteDimensionalError):
        PathAlgebra("~A2", orientation={"e12": (1, 2), "e20": (2, 0), "e01": (0, 1)})


def test_citations():
    assert "path_algebra" in PathAlgebra("A3").citations()
