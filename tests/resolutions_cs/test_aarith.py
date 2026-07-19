import pytest
from quiverlab import Quiver, CC
from quiverlab.groebner import build_reduction_system
from quiverlab.resolutions_cs.aarith import AArith
from quiverlab.resolutions_cs.ambiguities import SSequence
pytest.importorskip("quiverlab.groebner")


def _square():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    A = Q.algebra(relations=["a*b - c*d"], field=CC)
    return A, build_reduction_system(Q, ["a*b - c*d"], CC)


def test_path_vec_reduces_cd_to_ab():
    A, rs = _square()
    ar = AArith(A, rs)
    assert ar.path_vec(("c", "d")) == ar.path_vec(("a", "b"))       # cd normal-forms to ab
    assert ar.path_vec(("a", "b")) != [ar.dom.zero()] * A.dim


def test_square_cohomology_corner_dims():
    A, rs = _square()
    ar, ss = AArith(A, rs), SSequence(rs, 3)
    assert [sum(len(ar.corner(c.o, c.t, "coh")) for c in ss.S(n)) for n in range(3)] == [4, 4, 1]


def test_square_homology_corner_dims():
    A, rs = _square()
    ar, ss = AArith(A, rs), SSequence(rs, 3)
    assert [sum(len(ar.corner(c.o, c.t, "hom")) for c in ss.S(n)) for n in range(3)] == [4, 0, 0]


def test_kx2_corner_is_full_algebra():
    Q = Quiver([1], {"x": (1, 1)})
    A = Q.algebra(relations=["x*x"], field=CC)
    ar = AArith(A, build_reduction_system(Q, ["x*x"], CC))
    assert len(ar.corner(1, 1, "hom")) == A.dim                      # e_v A e_v = A, dim 2
