"""Recollement from an idempotent (Plan 47). The eAe-vs-subquiver trap is the pinned math;
the six functors' adjunction dim identities + the canonical exact sequences + j^*j_! ~= id
are the self-cert oracle (QPA white space). Char-clean (a GF(2) cell). Over kA3, the
commutative square, and a multi-vertex line algebra."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.modules.hom import hom_space, is_isomorphic
from quiverlab.modules.recollement import Recollement

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _a3(field=QQ):
    return linear_path_algebra(3, field=field)


def _square(field=QQ):
    # commutative square 1->2->4, 1->3->4 with ab = cd (one relation).
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=["a*b-c*d"], field=field)


def _line5(field=QQ):
    return linear_path_algebra(5, field=field)


def _a4(field=QQ):
    return linear_path_algebra(4, field=field)


@lit
def test_eae_is_not_the_subquiver_algebra():
    # S = {1,3}: the path ab:1->3 travels through 2 not in S but survives in eAe.
    A = _a3()
    R = Recollement(A, [1, 3])
    assert R.eAe.dim == 3                          # kA2, NOT the subquiver k x k (dim 2)
    assert len(list(R.eAe.quiver.vertices)) == 2


@lit
def test_worked_ka3_S2():
    A = _a3()
    R = Recollement(A, [2])
    assert R.eAe.dim == 1                          # e2 A e2 = k
    assert R.quotient.dim == 2                     # A/Ae2A = k x k on {1,3}
    P1 = A.projective(1)
    assert R.j_upper_star(P1).dim == 1             # P(1) e2 = vertex-2 part (1-dim)
    # i^*(P1) = P(1)/[2..3] = S_1 as an A/Ae2A-module = quotient.simple(1) (same category).
    assert is_isomorphic(R.i_upper_star(P1), R.quotient.simple(1))
    # and its inflation is the A-module S_1.
    assert is_isomorphic(R.i_star(R.i_upper_star(P1)), R.i_star(R.quotient.simple(1)))


@selfcert
def test_dim_certificates():
    A = _a3()
    for S in ([1], [2], [3], [1, 3], [1, 2]):
        R = Recollement(A, S)
        assert R.eAe.dim == R._corner_dim_sum()
        assert R.quotient.dim == A.dim - R._aeA_dim()


@selfcert
@pytest.mark.parametrize("build,S", [(_a3, [2]), (_square, [2, 3]), (_line5, [2, 4])])
def test_j_adjunctions_and_counit_iso(build, S):
    A = build()
    R = Recollement(A, S)
    verts = list(R.eAe.quiver.vertices) if R.eAe.quiver else []
    for xv in verts:
        X = R.eAe.simple(xv)
        for M in [A.projective(v) for v in A.quiver.vertices]:
            # (j_!, j^*):  dim Hom_A(j_! X, M) == dim Hom_eAe(X, j^* M)
            assert len(hom_space(R.j_shriek(X), M)) == len(hom_space(X, R.j_upper_star(M)))
            # (j^*, j_*):  dim Hom_A(M, j_* X) == dim Hom_eAe(j^* M, X)
            assert len(hom_space(M, R.j_star(X))) == len(hom_space(R.j_upper_star(M), X))
        # counit isos: j^* j_! X ~= X and j^* j_* X ~= X
        assert is_isomorphic(R.j_upper_star(R.j_shriek(X)), X)
        assert is_isomorphic(R.j_upper_star(R.j_star(X)), X)


@selfcert
def test_i_adjunctions():
    A = _a3()
    R = Recollement(A, [2])
    B = R.quotient
    for N in [B.simple(v) for v in B.quiver.vertices]:
        for M in [A.projective(v) for v in A.quiver.vertices]:
            # (i^*, i_*):  dim Hom_A(i^* -> ...): Hom_A(M, i_* N) == Hom_B(i^* M, N)
            assert len(hom_space(M, R.i_star(N))) == len(hom_space(R.i_upper_star(M), N))
            # (i_*, i^!):  dim Hom_A(i_* N, M) == dim Hom_B(N, i^! M)
            assert len(hom_space(R.i_star(N), M)) == len(hom_space(N, R.i_upper_shriek(M)))


@selfcert
@pytest.mark.parametrize("build,S", [(_a3, [2]), (_square, [2, 3]), (_line5, [2, 4]),
                                     (_a4, [2, 3])])
def test_canonical_sequences_exact(build, S):
    # BOTH BBD sequences, through the GENUINE natural maps (built from the actual functor
    # outputs, ModuleHom-certified as A-maps), across the algebra set + the new kA4 cell:
    #   counit:  j_! j^* M -> M -> i_* i^* M -> 0   (im(counit) = ker(unit_i), unit_i epi)
    #   unit:    0 -> i_* i^! M -> M -> j_* j^* M   (mono, im = ker(unit_j))
    A = build()
    R = Recollement(A, S)
    mods = [A.projective(v) for v in A.quiver.vertices]
    mods += [A.simple(v) for v in A.quiver.vertices]
    for M in mods:
        assert R.counit_sequence_exact(M), f"counit seq not exact for {M.name} in {build.__name__}"
        assert R.unit_sequence_exact(M), f"unit seq not exact for {M.name} in {build.__name__}"


@selfcert
def test_genuine_maps_land_where_hand_formulas_predict():
    # The genuine natural maps agree with the closed forms the OLD tautology used, but are
    # now built THROUGH the functors: im(counit_j) = M e A = M(AeA) = ker(unit_i).
    A = _square()
    R = Recollement(A, [2, 3])
    dom = A.domain
    from quiverlab.modules import linalg_mod as lm
    for M in [A.projective(v) for v in A.quiver.vertices]:
        eps, eta = R.counit_j(M), R.unit_i(M)
        Ie, _e, mono_e = eps.image()
        Ke, iota_k = eta.kernel()
        r_im = Ie.dim
        r_ker = Ke.dim
        # M(AeA) computed by hand for the arbiter (span, not identity of construction)
        hand = R._M_AeA_cols(M)
        assert r_im == r_ker == len(hand)
        assert eta.is_epi()                            # M ->> i_* i^* M


@selfcert
def test_degenerate_S_refused():
    # S = ALL vertices (e = 1) and S = EMPTY (e = 0) are degenerate recollements: refuse
    # loudly naming the idempotent, never a bare IndexError from a 0-vertex algebra build.
    A = _a3()
    allv = list(A.quiver.vertices)
    with pytest.raises(QuiverlabError, match="IDENTITY|full idempotent"):
        Recollement(A, allv)
    with pytest.raises(QuiverlabError, match="ZERO idempotent|empty"):
        Recollement(A, [])


@selfcert
def test_recollement_char_clean_gf2():
    A = _a3(GF(2))
    R = Recollement(A, [1, 3])
    assert R.eAe.dim == 3                           # no char caveat on the corner build
    # the j-side adjunction still holds over GF(2).
    X = R.eAe.simple(next(iter(R.eAe.quiver.vertices)))
    for M in [A.projective(v) for v in A.quiver.vertices]:
        assert len(hom_space(R.j_shriek(X), M)) == len(hom_space(X, R.j_upper_star(M)))
    assert is_isomorphic(R.j_upper_star(R.j_shriek(X)), X)
