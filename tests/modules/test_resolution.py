"""Minimal projective resolutions of right modules, generalized from the bridge engine
(spec §5 component 7). Fixtures A & B; any vertex set, any Domain."""
from quiverlab import Quiver, CC, GF, linear_path_algebra


def _square(field=CC):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def test_resolution_of_s1_in_a2():
    # 0 -> P_2 -> P_1 -> S_1 -> 0
    A = linear_path_algebra(2)
    res = A.simple(1).projective_resolution(4)
    assert res.term(0) == [1]        # Q_0 = P_1
    assert res.term(1) == [2]        # Q_1 = P_2 (rad P_1 = S_2 = P_2)
    assert res.betti(0) == 1 and res.betti(1) == 1
    assert res.term(2) == []         # terminates: Omega_2 = 0
    assert res.pd() == 1


def test_projective_is_its_own_resolution():
    A = linear_path_algebra(2)
    res = A.projective(1).projective_resolution(3)
    assert res.term(0) == [1] and res.term(1) == []
    assert res.pd() == 0


def test_simple_s2_is_projective_in_a2():
    A = linear_path_algebra(2)
    assert A.simple(2).projective_resolution(3).pd() == 0    # S_2 = P_2


def test_square_s1_has_pd_two():
    A = _square()
    res = A.simple(1).projective_resolution(5)
    assert res.pd() == 2             # gl.dim(commutative square) = 2
    # d_n(Q_n) subset rad Q_{n-1}: minimality => betti numbers are the true Betti numbers
    assert res.betti(0) == 1
    assert res.term(3) == []


def test_resolution_over_gfp_matches_char0_shape():
    A = _square(field=GF(7))
    assert A.simple(1).projective_resolution(5).pd() == 2


def test_resolution_differentials_compose_to_zero():
    A = _square()
    res = A.simple(1).projective_resolution(4)
    from quiverlab.modules import linalg_mod as lm
    for n in range(1, res.length):
        Dn, Dn1 = res.differential(n), res.differential(n + 1)
        if not Dn or not Dn1 or len(Dn1[0]) == 0:
            continue
        prod = lm.matmul(Dn, Dn1, A.domain)
        assert all(A.domain.is_zero(x) for row in prod for x in row)
