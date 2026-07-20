"""Hom/End of right A-modules over any Domain (spec §3.6)."""
from quiverlab import Quiver, CC, GF, linear_path_algebra
from quiverlab.modules.module import Module


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


def test_hom_projective_to_module_reads_off_NONUNIFORM_dimension_vector():
    # Strengthening of the test above: that one pins Hom(P_v, N) against a target N
    # (= P_1 over the square) whose dimvec is UNIFORM (all 1s), so it cannot tell the
    # Yoneda law "dim Hom(P_v, N) = dimvec(N)[v]" apart from a broken "always 1".
    # Here N has a genuinely non-uniform dimvec, with one coordinate = 2.
    #
    # Algebra: kA_3, the hereditary path algebra of 1 -a1-> 2 -a2-> 3 (no relations).
    # Target:  N = P_1 (+) S_2.
    #   * P_1 = e_1 A has k-basis {e_1, a1, a1*a2} sitting at vertices {1, 2, 3},
    #     so dimvec(P_1) = {1:1, 2:1, 3:1}.
    #   * S_2 is the simple at vertex 2, dimvec(S_2) = {1:0, 2:1, 3:0}.
    #   * hence dimvec(N) = {1:1, 2:2, 3:1}  -- two distinct values, the value 2 at v=2.
    #
    # We build N directly via from_arrow_action in the vertex-ordered basis
    #   b0 = e_1(P_1) [v1],  b1 = a1(P_1) [v2],  b2 = s2(S_2) [v2],  b3 = a1*a2(P_1) [v3].
    # Right action (m*b = action[b] @ m): a1 sends the P_1 generator b0 -> b1, and a2
    # sends b1 -> b3; both arrows kill the S_2 vector b2. from_arrow_action validates
    # that this is a genuine module (grading + multiplicativity) before returning.
    #
    # Yoneda: dim Hom(P_v, N) = dim(N e_v) = dimvec(N)[v], so the expected dims are
    #   Hom(P_1, N) = 1,  Hom(P_2, N) = 2,  Hom(P_3, N) = 1.
    # The v=2 assert (dim 2, not 1) is what a "always 1" bug would fail.
    A = linear_path_algebra(3)
    Z = [[0] * 4 for _ in range(4)]
    a1 = [row[:] for row in Z]; a1[1][0] = 1   # b0 -> b1
    a2 = [row[:] for row in Z]; a2[3][1] = 1   # b1 -> b3
    N = Module.from_arrow_action(A, {1: 1, 2: 2, 3: 1}, {"a1": a1, "a2": a2},
                                 name="P1(+)S2")
    dv = N.dimension_vector()
    assert dv == {1: 1, 2: 2, 3: 1}            # non-uniform, one coordinate >= 2
    expected = {1: 1, 2: 2, 3: 1}
    for v in (1, 2, 3):
        assert A.hom(A.projective(v), N) == expected[v] == dv[v]


def test_end_of_indecomposable_projective_is_local_dim_one_top():
    A = linear_path_algebra(2)
    from quiverlab.modules.hom import end_dim
    # End(P_1): P_1 indecomposable; endomorphisms are scalars + radical -> dim 1 here
    assert end_dim(A.projective(1)) == 1


def test_hom_over_gfp():
    A = linear_path_algebra(2, field=GF(7))
    assert A.hom(A.simple(1), A.simple(1)) == 1
