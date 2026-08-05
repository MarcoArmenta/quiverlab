"""One-point extension A[M] = [[k,M],[0,A]] (Plan 44 / C7). Certified per instance by
dim A[M] == 1 + dim M + dim A; the Cartan block form and the pd(S_w) = pd_A(M) + 1
identity are the literature oracles (ASS III / Happel)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.families.one_point import OnePointExtension

pytestmark = pytest.mark.oracle_literature


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def _new_vertex(B, A):
    """The extra source vertex of A[M] = the B-vertex not among A's vertices."""
    aset = set(A.quiver.vertices)
    return next(v for v in B.quiver.vertices if v not in aset)


def _strip_new_vertex(CB):
    """Strip the new source vertex's row/col (it is placed FIRST) and return the A-block."""
    return [row[1:] for row in CB[1:]]


def test_dim_certificate_and_quiver():
    # A[S1] for A = kA2: new source w, arrow w->1 (top S1 at vertex 1), and the relation
    # (w->1->2) = 0 because a acts as 0 on S1. Q' = w->1->2 (A3 shape) with one relation.
    A = _kA2()
    S1 = A.simple(1)
    B = OnePointExtension(A, S1)
    assert B.dim == 1 + S1.dim + A.dim == 5          # NOT kA3 (dim 6): the relation cuts 1
    assert len(list(B.quiver.vertices)) == 3
    assert len(B.relations) >= 1                      # the w->1->2 = 0 relation


def test_cartan_block_form():
    # C_{A[M]} = [[1, dimvec M], [0, C_A]] (new vertex FIRST; the Cartan equality arbitrates
    # the source-vertex placement + row/col convention).
    A = _kA2()
    S1 = A.simple(1)
    B = OnePointExtension(A, S1)
    CB, CA = B.cartan_matrix(), A.cartan_matrix()
    assert _strip_new_vertex(CB) == CA               # A-block matches C_A
    # the new source row is [1, dim-vector M] (w -> v paths = M_v), the column below it is 0.
    dv = S1.dimension_vector()
    assert CB[0] == [1] + [dv[v] for v in A.quiver.vertices]
    assert [CB[i][0] for i in range(1, len(CB))] == [0] * (len(CB) - 1)


def test_pd_of_new_simple_is_pd_M_plus_one():
    # pd_{A[M]}(S_w) = pd_A(M) + 1. For M = S1 over kA2, pd_A(S1) = 1, so pd = 2.
    A = _kA2()
    S1 = A.simple(1)
    B = OnePointExtension(A, S1)
    w = _new_vertex(B, A)                              # the source vertex of A[M]
    pd_new = B.simple(w).projective_resolution(8).pd()
    pd_M = S1.projective_resolution(8).pd()
    assert pd_new == pd_M + 1 == 2


@pytest.mark.oracle_selfcert
def test_projective_module_extension_has_no_relation():
    # M = S2 (= P2, projective, pd 0): A[S2] has quiver w->2, 1->2 and NO length-2
    # relation; dim 5, pd(S_w) = 1.
    A = _kA2()
    S2 = A.simple(2)
    B = OnePointExtension(A, S2)
    assert B.dim == 5
    w = _new_vertex(B, A)
    assert B.simple(w).projective_resolution(8).pd() == 1
