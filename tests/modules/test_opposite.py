"""The opposite algebra A^op (Plan 23): reversed quiver, transposed structure
constants, reversed labels/relations, involutive. Monomial and non-monomial
(Groebner) presentations, single- and multi-vertex."""
import pytest

from quiverlab import Quiver, CC, GF, linear_path_algebra, truncated_polynomial
from quiverlab.errors import QuiverlabError


def _square(field=CC):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


@pytest.mark.parametrize("A", [linear_path_algebra(3), truncated_polynomial(3), _square(),
                               _square(GF(5))])
def test_opposite_transposes_structure_constants(A):
    Aop = A.opposite()
    m = A.dim
    dom = A.domain
    for i in range(m):
        for j in range(m):
            assert all(dom.is_zero(dom.sub(Aop.T[i][j][t], A.T[j][i][t])) for t in range(m))
    Aop._validate()                     # A^op is a genuine associative unital algebra


def test_opposite_reverses_arrows_and_labels():
    A = linear_path_algebra(3)          # 1 -a1-> 2 -a2-> 3
    Aop = A.opposite()
    assert Aop.quiver.arrows == {"a1": (2, 1), "a2": (3, 2)}
    # the length-2 path a1*a2 (1->3) reverses to a2*a1 (3->1) in A^op
    assert "a2*a1" in Aop.basis_labels
    assert "a1*a2" not in Aop.basis_labels


def test_opposite_is_cached_involution():
    A = _square()
    assert A.opposite().opposite() is A
    B = truncated_polynomial(4, field=GF(3))
    assert B.opposite().opposite() is B


@pytest.mark.parametrize("field", [CC, GF(7)])
def test_opposite_modules_are_genuine(field):
    A = _square(field=field)
    Aop = A.opposite()
    for v in Aop.quiver.vertices:
        ok, why = Aop.projective(v).check_module()
        assert ok, why
        ok, why = Aop.simple(v).check_module()
        assert ok, why


def test_opposite_needs_provenance():
    from quiverlab import Algebra
    T = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]
    A = Algebra.from_structure_constants(T, unit=[1, 0], field=CC)
    with pytest.raises(QuiverlabError):
        A.opposite()


def test_opposite_dim_and_cartan_transpose():
    # dim e_v A^op e_w = # A^op-paths v->w = # A-paths w->v = C_A[w][v]:  Cartan(A^op) = Cartan(A)^T
    A = _square()
    Aop = A.opposite()
    C = A.cartan_matrix()
    Cop = Aop.cartan_matrix()
    n = len(C)
    assert Cop == [[C[j][i] for j in range(n)] for i in range(n)]
