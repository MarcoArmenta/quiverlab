"""Minimal add(M)-approximations (Plan 44 / C7). Self-certifying: the approximation
property is checked by solve_columns surjectivity of Hom(M, M^k) -> Hom(M, C); minimality
by no-proper-summand-restriction. The projective-cover tie is the literature anchor."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.modules.approximations import (_left_approx_blocks, _prune,
                                              _right_approx_blocks,
                                              left_add_approximation,
                                              right_add_approximation)
from quiverlab.modules.morphism import direct_sum, hom_basis

pytestmark = pytest.mark.oracle_selfcert


def _kA2():
    # QQ so the decompose char scope (the minimal-approximation summand set) is rigorous.
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


def _prune_removes_nothing(blocks_fn, M, C):
    """Pruning-idempotence self-cert (Medium 1): prune the FULL approximation to its fixed
    point, then prune the PRUNED block set again -- it must remove NOTHING and stay an
    approximation. Replaces the old tautological 'k == k' minimality check, which compared
    the constructor's output to itself."""
    _types, blocks, is_approx = blocks_fn(M, C, 512)
    _prune(blocks, is_approx)
    n1 = len(blocks)
    _prune(blocks, is_approx)                              # prune the already-pruned set
    return len(blocks) == n1 and is_approx(blocks)


def _factors_every_hom(f, M, C):
    """Every h in Hom(M, C) equals f-after-(some M -> M^k): the approximation property.
    Reduce to a solve over the columns { (h_j composed into f's coordinates) }."""
    from quiverlab.modules import linalg_mod as lm
    dom = C.domain
    homs = hom_basis(M, C)
    # columns spanning im(Hom(M, M^k) -> Hom(M, C)) via f: each basis phi: M -> M^k
    # gives phi.then(f): M -> M^k -> C.  build the induced-map image.
    from quiverlab.modules.morphism import hom_basis as hb
    induced = [phi.then(f) for phi in hb(M, f.src)]        # M -> M^k -> C
    def vec(g):  # g: M -> C, flatten
        return [g.matrix[i][j] for j in range(M.dim) for i in range(C.dim)]
    B = lm.cols_to_matrix([vec(g) for g in induced]) if induced else []
    for h in homs:
        if not B or lm.solve_columns(B, lm.cols_to_matrix([vec(h)]), dom) is None:
            return False
    return True


def test_projective_cover_is_the_min_add_P_approximation():
    # add(P1) contains P(S1) = P1, so the min right add(P1)-approximation of S1 IS the
    # projective cover P1 ->> S1 (ASS I.5). Anchors approximation = cover.
    A = _kA2()
    P1, S1 = A.projective(1), A.simple(1)
    f = right_add_approximation(P1, S1)
    assert f.is_epi() and f.tgt is S1
    K, _ = f.kernel()
    assert K.dimension_vector() == P1.radical().dimension_vector()   # cover: ker in rad


def test_radical_inclusion_is_the_min_add_S2_approximation():
    # Hom(S2, P1) = k (S2 = rad P1 >-> P1); the min right add(S2)-approx of P1 is that
    # mono, NOT surjective (P1/S2 = S1 != 0). A genuine non-epi approximation.
    A = _kA2()
    S2, P1 = A.simple(2), A.projective(1)
    f = right_add_approximation(S2, P1)
    assert f.is_mono() and not f.is_epi() and f.src.dim == 1


@pytest.mark.parametrize("mv, cv", [(1, 1), (1, 2), (2, 1)])
def test_approximation_property_and_minimality_battery(mv, cv):
    A = _kA2()
    pool = {1: A.projective(1), 2: A.projective(2)}
    M, C = pool[mv], A.simple(cv)
    f = right_add_approximation(M, C)
    assert _factors_every_hom(f, M, C)                     # approximation property
    # minimality (Medium 1): a REAL assertion, not the old 'k == k' tautology --
    # pruning the pruned approximation removes nothing (irredundant = minimal).
    assert _prune_removes_nothing(_right_approx_blocks, M, C)


def test_decomposable_M_source_is_not_inflated():
    # THE HIGH bug: over a DECOMPOSABLE M the old code multiplied the WHOLE M by the
    # generator count, so it could never drop a superfluous summand of M. reg = P1 (+) P2:
    # the min right add(reg)-approx of S1 is the projective cover P1 (dv [1,1]), NOT reg
    # (dv [1,2]); of reg it is reg itself (dv [1,2]), NOT reg^2. QPA-arbitrated on the
    # decomposable-M cases in tests/qpa/test_tilting_qpa.py.
    A = _kA2()
    reg, _, _ = direct_sum(A.projective(1), A.projective(2))
    fS = right_add_approximation(reg, A.simple(1))
    assert fS.src.dimension_vector() == A.projective(1).dimension_vector()   # P1 cover
    assert _factors_every_hom(fS, reg, A.simple(1))
    fR = right_add_approximation(reg, reg)
    assert fR.src.dimension_vector() == reg.dimension_vector()               # reg, not reg^2
    assert _factors_every_hom(fR, reg, reg)
    # pruning-idempotence self-cert on the decomposable cases (Medium 1).
    assert _prune_removes_nothing(_right_approx_blocks, reg, A.simple(1))
    assert _prune_removes_nothing(_right_approx_blocks, reg, reg)
    assert _prune_removes_nothing(_left_approx_blocks, reg, A.simple(2))


def test_left_is_dual_of_right():
    A = _kA2()
    M, C = A.projective(1), A.simple(2)
    g = left_add_approximation(M, C)
    assert g.src is C                                      # C >-> M^k
    # every C -> M factors through g (dual approximation property)
    from quiverlab.modules.morphism import hom_basis as hb
    from quiverlab.modules import linalg_mod as lm
    dom = C.domain
    induced = [g.then(psi) for psi in hb(g.tgt, M)]        # C -> M^k -> M
    def vec(h):
        return [h.matrix[i][j] for j in range(C.dim) for i in range(M.dim)]
    B = lm.cols_to_matrix([vec(h) for h in induced]) if induced else []
    for h in hb(C, M):
        assert B and lm.solve_columns(B, lm.cols_to_matrix([vec(h)]), dom) is not None


def test_cross_algebra_refused():
    A, B = _kA2(), _kA2()
    from quiverlab.errors import QuiverlabError
    with pytest.raises(QuiverlabError):
        right_add_approximation(A.simple(1), B.simple(1))
