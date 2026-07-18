import pytest
from quiverlab import CC, GF, Quiver
from quiverlab.errors import DepthLimitError


def _dual(field):
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    return Q.algebra(relations=["x^2"], field=field)


def test_dual_numbers_char0():
    # classical: HH^0 = 2 (commutative), HH^n = 1 for n >= 1 in char 0
    A = _dual(CC)
    t = A.hochschild_cohomology(4)
    assert t.dims == [2, 1, 1, 1, 1]
    assert t.kind == "HH^"


def test_dual_numbers_char2_pathology():
    # char 2: every differential vanishes -> HH^n = 2 for all n
    A = _dual(GF(2))
    assert A.hochschild_cohomology(4).dims == [2, 2, 2, 2, 2]


def test_dual_numbers_gf4_matches_gf2():
    A = _dual(GF(4))
    assert A.hochschild_cohomology(3).dims == [2, 2, 2, 2]


def test_semisimple_k_times_k():
    Q = Quiver(vertices=[1, 2], arrows={})
    A = Q.algebra()
    assert A.hochschild_cohomology(3).dims == [2, 0, 0, 0]


def test_hereditary_kA2():
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2)})
    A = Q.algebra()
    assert A.hochschild_cohomology(3).dims == [1, 0, 0, 0]


def test_guard_fails_loudly():
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x^4"], field=GF(3))
    with pytest.raises(DepthLimitError):
        A.hochschild_cohomology(30, max_cells=1000)


def test_repr_is_a_table():
    A = _dual(CC)
    s = repr(A.hochschild_cohomology(2))
    assert "HH^0" in s and "HH^2" in s
