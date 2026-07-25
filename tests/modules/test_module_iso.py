"""Exact module isomorphism test (Plan 23, modules/hom.py::is_isomorphic).

Positive answers carry an invertible-map certificate (never a false positive);
negatives are certain (distinct dim / dim vector, or exhaustive over small GF(p))."""
import pytest

from quiverlab import Quiver, CC, GF, linear_path_algebra


def _square(field=CC):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


@pytest.mark.parametrize("field", [CC, GF(2), GF(7)])
def test_reflexive_and_distinct(field):
    A = linear_path_algebra(3, field=field)
    for v in (1, 2, 3):
        assert A.simple(v).is_isomorphic(A.simple(v))
        assert A.projective(v).is_isomorphic(A.projective(v))
    assert not A.simple(1).is_isomorphic(A.simple(2))     # distinct dim vectors
    assert not A.simple(1).is_isomorphic(A.projective(1)) # distinct dims


def test_iso_same_dimvec_different_module():
    # In kA_2, S_1 (+) S_2 and P_1 share the dimension vector {1:1,2:1} but are NOT
    # isomorphic: P_1 is indecomposable (arrow acts nonzero), the sum is semisimple.
    A = linear_path_algebra(2)
    from quiverlab.modules.module import Module
    # S_1 (+) S_2: dimvec {1:1,2:1}, arrow a1 acts as ZERO (semisimple impostor).
    ss = Module.from_arrow_action(A, {1: 1, 2: 1}, {"a1": [[0, 0], [0, 0]]}, name="S1+S2")
    P1 = A.projective(1)
    assert ss.dimension_vector() == P1.dimension_vector()
    assert not ss.is_isomorphic(P1)                       # exact search refutes over CC


@pytest.mark.parametrize("field", [CC, GF(3)])
def test_iso_across_builders(field):
    # S_2 = P_2 and S_1 = I_1 in kA_2 (standard identifications).
    A = linear_path_algebra(2, field=field)
    assert A.simple(2).is_isomorphic(A.projective(2))
    assert A.simple(1).is_isomorphic(A.injective(1))


def test_iso_certificate_over_gfp_square():
    # I_v = D(A e_v): A e_v is the LEFT projective, D flips it to the right injective.
    # (Plan 24: D is side-aware; the pre-Plan-24 form A.opposite().projective(v).dualize()
    # now yields a left-tagged module, so the honest expression uses side="left".)
    A = _square(field=GF(2))
    for v in A.quiver.vertices:
        assert A.injective(v).is_isomorphic(A.projective(v, side="left").dualize())
