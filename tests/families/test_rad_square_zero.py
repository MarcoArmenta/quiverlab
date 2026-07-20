"""RadicalSquareZero(Q) = kQ/rad^2. Fixture N4."""
from quiverlab.combinat import Quiver
from quiverlab.families import RadicalSquareZero
from quiverlab.fields import CC


def test_two_loop_is_finite_dim3():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})     # kQ infinite; rad^2 truncates
    A = RadicalSquareZero(Q, field=CC)
    assert A.dim == 3                               # e, x, y  (Fixture N4)
    assert A.hochschild_cohomology(0).dims == [3]   # commutative local: HH^0 = 3


def test_three_cycle_dim6_hh0_1():
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)})
    A = RadicalSquareZero(Q, field=CC)
    assert A.dim == 6                               # 3 vertices + 3 arrows
    assert A.hochschild_cohomology(0).dims == [1]
