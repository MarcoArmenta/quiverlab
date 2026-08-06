"""tau-rigidity + g-vectors (Plan 45 / C4). Self-cert: g^{P_v}=e_v, g additive,
tau-rigid one-call. Literature: over a hereditary algebra tau-rigid <=> rigid
(Ext^1(M,M)=0); the g-matrix of (A,0) is the identity (det 1)."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.rigid import g_matrix, g_vector, is_tau_rigid

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)


@selfcert
def test_projectives_are_tau_rigid_and_g_is_unit():
    A = _kA2()
    for v in (1, 2):
        Pv = A.projective(v)
        assert is_tau_rigid(Pv)                       # tau P = 0
        g = g_vector(Pv)
        assert g == {1: (1 if v == 1 else 0), 2: (1 if v == 2 else 0)}   # e_v


@selfcert
def test_g_vector_is_additive():
    A = _kA2()
    from quiverlab.modules.morphism import direct_sum
    D, _, _ = direct_sum(A.projective(1), A.simple(1))
    g1, gS1 = g_vector(A.projective(1)), g_vector(A.simple(1))
    assert g_vector(D) == {v: g1[v] + gS1[v] for v in (1, 2)}


@selfcert
def test_simple_at_source_g_vector():
    # S_1 over kA2: min presentation P_2 -> P_1 -> S_1 (P_1 = [1,2], rad = S_2 = P_2),
    # so 0 -> P_2 -> P_1 -> S_1 -> 0 ; g^{S_1} = e_1 - e_2 = {1: 1, 2: -1}.
    A = _kA2()
    assert g_vector(A.simple(1)) == {1: 1, 2: -1}


@lit
def test_hereditary_tau_rigid_iff_rigid():
    # over a hereditary algebra (kA3, linear) tau-rigid <=> Ext^1(M,M) = 0 (ASS/AIR).
    A = linear_path_algebra(3, field=QQ)
    for v in (1, 2, 3):
        M = A.simple(v)
        rigid = (A.ext(M, M, 1) == 0)
        assert is_tau_rigid(M) is rigid
    # and a genuine non-rigid check on a decomposable self-ext witness if one exists
    # (kA3 simples are bricks; the identity holds trivially -- keep the loop as the pin).


@lit
def test_initial_pair_g_matrix_is_identity():
    from quiverlab.tautilting.pairs import initial_pair
    A = linear_path_algebra(3, field=QQ)
    G = g_matrix(initial_pair(A))
    assert G == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]         # (A,0): columns e_1,e_2,e_3
