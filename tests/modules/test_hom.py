"""Hom/End of right A-modules over any Domain (spec §3.6)."""
from quiverlab import Quiver, CC, GF, linear_path_algebra


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"])


def test_hom_of_simples_is_kronecker_delta():
    A = linear_path_algebra(2)
    S1, S2 = A.simple(1), A.simple(2)
    assert A.hom(S1, S1) == 1
    assert A.hom(S2, S2) == 1
    assert A.hom(S1, S2) == 0
    assert A.hom(S2, S1) == 0


def test_hom_projective_to_module_reads_off_dimension_vector():
    # dim Hom(P_v, N) = dim (N e_v) = dimension_vector(N)[v]  (Yoneda for projectives)
    A = _square()
    N = A.projective(1)
    for v in (1, 2, 3, 4):
        assert A.hom(A.projective(v), N) == N.dimension_vector()[v]


def test_end_of_indecomposable_projective_is_local_dim_one_top():
    A = linear_path_algebra(2)
    from quiverlab.modules.hom import end_dim
    # End(P_1): P_1 indecomposable; endomorphisms are scalars + radical -> dim 1 here
    assert end_dim(A.projective(1)) == 1


def test_hom_over_gfp():
    A = linear_path_algebra(2, field=GF(7))
    assert A.hom(A.simple(1), A.simple(1)) == 1
