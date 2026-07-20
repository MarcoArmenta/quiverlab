"""Minimal projective resolutions of right modules, generalized from the bridge engine
(spec §5 component 7). Fixtures A & B; any vertex set, any Domain."""
from quiverlab import Quiver, CC, GF, linear_path_algebra


def _square(field=CC):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def _kx3(field=CC):
    # k[x]/(x^3): one vertex, loop x, monomial relation x*x*x. Self-injective, so the
    # simple S has infinite projective dimension (periodic resolution).
    return Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x*x"], field=field)


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
    res = A.simple(1).projective_resolution(5)
    assert res.pd() == 2
    # d_n . d_{n+1} = 0 must hold over GF(7) too; the standalone compose-to-zero
    # test (test_resolution_differentials_compose_to_zero) only exercises CC.
    from quiverlab.modules import linalg_mod as lm
    for n in range(1, res.length):
        Dn, Dn1 = res.differential(n), res.differential(n + 1)
        if not Dn or not Dn1 or len(Dn1[0]) == 0:
            continue
        prod = lm.matmul(Dn, Dn1, A.domain)
        assert all(A.domain.is_zero(x) for row in prod for x in row)


def test_kx3_depth_limit_error_certifies_length():
    # k[x]/(x^3)'s simple S has syzygy Omega_1 = rad P_1 of dim 2, so any max_term_dim
    # below 2 must abort loudly with the certified length in the message (not silently).
    import pytest
    from quiverlab.errors import DepthLimitError
    from quiverlab.modules.resolution import minimal_resolution
    S = _kx3().simple(1)
    with pytest.raises(DepthLimitError) as exc:
        minimal_resolution(S, 6, max_term_dim=1)      # dim-2 syzygy overshoots the bound
    assert "certified through length 0" in str(exc.value)


def test_kx3_depth_limit_error_via_facade():
    # Fix 2: Module.projective_resolution now threads max_term_dim through to
    # minimal_resolution, so the DepthLimitError guard is reachable from the facade
    # (previously only from a direct minimal_resolution call). Default path unchanged.
    import pytest
    from quiverlab.errors import DepthLimitError
    S = _kx3().simple(1)
    with pytest.raises(DepthLimitError) as exc:
        S.projective_resolution(6, max_term_dim=1)     # dim-2 syzygy overshoots the bound
    assert "certified through length 0" in str(exc.value)
    # default path is unchanged (still resolves periodically to the requested length)
    res = S.projective_resolution(4)
    assert res.betti(0) == 1


def test_kx3_simple_has_infinite_pd():
    # k[x]/(x^3) is self-injective and non-semisimple: S resolves periodically forever.
    # pd() is None (never resolved), is_finite() is False, and every term is a single P_1.
    res = _kx3().simple(1).projective_resolution(8)
    assert res.pd() is None
    assert res.is_finite() is False
    assert all(res.betti(n) == 1 for n in range(res.length))


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
