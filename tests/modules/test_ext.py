"""Module Ext^n and global dimension (spec §3.5, §3.6). Literature oracles: ASS III /
Happel for hereditary A_n; Kunneth for the commutative square; bank diamond for Ext^2."""
import pytest
from quiverlab import Quiver, CC, GF, linear_path_algebra


def _square(field=CC):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def _mon_diamond():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b"])


def test_ext_of_simples_a3_hereditary():
    A = linear_path_algebra(3)          # 1 -> 2 -> 3
    S = {i: A.simple(i) for i in (1, 2, 3)}
    # Hom(S_a,S_b) = delta
    for a in (1, 2, 3):
        for b in (1, 2, 3):
            assert A.ext(S[a], S[b], 0) == (1 if a == b else 0)
    # Ext^1(S_a, S_{a+1}) = 1 ; all other Ext^1 = 0  (ASS III, Happel)
    assert A.ext(S[1], S[2], 1) == 1
    assert A.ext(S[2], S[3], 1) == 1
    assert A.ext(S[1], S[3], 1) == 0
    assert A.ext(S[1], S[1], 1) == 0
    # hereditary: Ext^{>=2} = 0
    assert A.ext(S[1], S[3], 2) == 0


def test_ext2_commutative_square_kunneth():
    A = _square()
    assert A.ext(A.simple(1), A.simple(4), 2) == 1      # kA_2 (x) kA_2, Kunneth
    assert A.ext(A.simple(1), A.simple(4), 1) == 0
    assert A.ext(A.simple(1), A.simple(4), 3) == 0


def test_ext2_monomial_diamond_bank_oracle():
    # bridge module_ext oracle: kD/(a*b) has Ext^2(S_1, S_4) = 1
    A = _mon_diamond()
    assert A.ext(A.simple(1), A.simple(4), 2) == 1


def test_projectives_have_no_higher_ext():
    A = _square()
    for v in (1, 2, 3, 4):
        assert A.ext(A.projective(v), A.simple(1), 1) == 0
        assert A.ext(A.projective(v), A.simple(1), 2) == 0


def test_global_dimension():
    assert int(linear_path_algebra(2).global_dimension()) == 1     # hereditary
    assert int(linear_path_algebra(3).global_dimension()) == 1
    gd = _square().global_dimension()
    assert int(gd) == 2 and gd.exact is True
    assert gd == 2


def test_global_dimension_selfinjective_is_infinite_bound():
    # k[x]/(x^2) is self-injective: pd(S) = infinity -> certified lower bound, not exact
    from quiverlab import truncated_polynomial
    gd = truncated_polynomial(2).global_dimension()
    assert gd.exact is False and gd.value >= 1


def test_is_selfinjective():
    from quiverlab import truncated_polynomial
    # k[x]/(x^n) is (Frobenius, hence) self-injective, over any field
    assert truncated_polynomial(2).is_selfinjective() is True
    assert truncated_polynomial(3).is_selfinjective() is True
    assert truncated_polynomial(3, field=GF(5)).is_selfinjective() is True
    # hereditary kA_2 is NOT self-injective (soc P_1 = S_2 = soc P_2, not a bijection)
    assert linear_path_algebra(2).is_selfinjective() is False
    # the commutative square is not self-injective either
    assert _square().is_selfinjective() is False


def test_selfinjective_agrees_with_frobenius_over_gfp():
    # over GF(p) the exact is_selfinjective must agree with the engine-backed is_frobenius
    from quiverlab import truncated_polynomial
    A = truncated_polynomial(2, field=GF(5))
    assert A.is_selfinjective() == A.is_frobenius() == True


def test_ext_over_gfp():
    A = _square(field=GF(7))
    assert A.ext(A.simple(1), A.simple(4), 2) == 1
