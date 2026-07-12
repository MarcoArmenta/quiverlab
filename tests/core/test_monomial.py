import pytest
from quiverlab import CC, GF, Quiver
from quiverlab.errors import AdmissibilityError, NotFiniteDimensionalError


def test_dual_numbers_from_quiver():
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x^2"], field=GF(2))
    assert A.dim == 2                      # e_1, x
    assert "x" in (A.basis_labels or [])


def test_truncated_loop():
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x^4"])
    assert A.dim == 4                      # e, x, x^2, x^3


def test_a3_rad_square_zero():
    Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)})
    A = Q.algebra(relations=["a*b"])
    assert A.dim == 5                      # e1, e2, e3, a, b


def test_hereditary_no_relations():
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2)})
    A = Q.algebra()
    assert A.dim == 3                      # e1, e2, a


def test_loop_without_relations_is_infinite():
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    with pytest.raises(NotFiniteDimensionalError) as ei:
        Q.algebra()
    assert "x" in str(ei.value)            # the offending cycle is named


def test_two_loops_one_relation_still_infinite():
    # k<x, y>/(xy): the words y^j x^i are all irreducible -> infinite-dimensional,
    # and the automaton must name a cycle (a loop on x or on y).
    Q = Quiver(vertices=[1], arrows={"x": (1, 1), "y": (1, 1)})
    with pytest.raises(NotFiniteDimensionalError) as ei:
        Q.algebra(relations=["x*y"])
    msg = str(ei.value)
    assert ("x" in msg) or ("y" in msg)


def test_short_relation_not_admissible():
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2)})
    with pytest.raises(AdmissibilityError):
        Q.algebra(relations=["a"])


def test_nonmonomial_waits_for_plan03():
    Q = Quiver(vertices=[1, 2, 3],
               arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})
    with pytest.raises(NotImplementedError):
        Q.algebra(relations=["a*b - c"])


def test_multiplication_table_is_path_concatenation():
    Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)})
    A = Q.algebra(relations=[], field=CC)  # kA_3: e1,e2,e3,a,b,a*b -> dim 6
    assert A.dim == 6
    labels = A.basis_labels
    ia, ib, iab = labels.index("a"), labels.index("b"), labels.index("a*b")
    dom = A.domain
    va = [dom.one() if i == ia else dom.zero() for i in range(6)]
    vb = [dom.one() if i == ib else dom.zero() for i in range(6)]
    assert A.multiply(va, vb)[iab] == dom.one()
    assert all(dom.is_zero(c) for i, c in enumerate(A.multiply(vb, va)))
