"""TruncatedPathAlgebra kQ/rad^r. Fixture N3."""
from quiverlab.combinat import Quiver
from quiverlab.families import TruncatedPathAlgebra
from quiverlab.fields import CC


def test_A3_rad2_dim5():
    assert TruncatedPathAlgebra("A3", 2, field=CC).dim == 5     # Fixture N3
    assert TruncatedPathAlgebra("A5", 2, field=CC).dim == 9     # spec example


def test_from_explicit_quiver():
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)})
    assert TruncatedPathAlgebra(Q, 2, field=CC).dim == 5
