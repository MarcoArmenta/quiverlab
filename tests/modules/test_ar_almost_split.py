"""Almost-split sequences. Self-cert: exact + non-split + indecomposable ends.
Literature: the kA_n AR sequences have the mesh middle-term dimension vectors."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.modules.hom import is_isomorphic

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _a3():
    return linear_path_algebra(3, field=QQ)


@selfcert
def test_almost_split_is_exact_nonsplit_indecomposable_ends():
    A = _a3()
    M = A.simple(1)                                 # S1 = [1], indecomposable, non-projective
    assert M.is_indecomposable()
    ses = M.almost_split_sequence()                 # 0 -> tau M -> E -> M -> 0
    assert ses.N.is_isomorphic(M) or ses.N.dim == M.dim   # right end is M
    assert ses.M.dim == ses.L.dim + M.dim           # exact (P37 SES certifies)
    assert ses.is_split() is False                  # almost-split => non-split
    from quiverlab.modules.decompose import is_indecomposable
    assert is_indecomposable(ses.L)                 # tau M indecomposable
    assert is_indecomposable(M)


@selfcert
def test_projective_input_refused():
    A = _a3()
    with pytest.raises(QuiverlabError, match="projective"):
        A.projective(1).almost_split_sequence()     # no AR sequence ends at a projective


@selfcert
def test_decomposable_input_refused():
    A = _a3()
    from quiverlab.modules.morphism import direct_sum
    D, _, _ = direct_sum(A.simple(1), A.simple(2))
    with pytest.raises(QuiverlabError, match="indecomposable"):
        from quiverlab.modules.ar import almost_split_sequence
        almost_split_sequence(D)


@selfcert
def test_kx4_non_brick_middle_is_mesh_not_projective():
    # Devil's-advocate round (Finding B): the DISCRIMINATING non-brick case. Over k[x]/(x^4)
    # (GF(32003)), M = k[x]/(x^2) -- the length-2 module, x acting as the lower Jordan block
    # J_2 = [[0,0],[1,0]] -- is indecomposable but NOT a brick: End(M) = k[x]/(x^2) is local
    # of k-dimension 2, so dim Ext^1(M, tau M) = 2 and the three sanity guards (exact +
    # non-split + indecomposable ends) are necessary but NOT sufficient. A non-socle class of
    # the 2-dimensional Ext^1 has the PROJECTIVE middle k[x]/(x^4) = {4} (exact, non-split,
    # indecomposable ends) yet is NOT almost split. Only the socle-simplicity certificate
    # forces the true mesh middle {1, 3}.
    from quiverlab import truncated_polynomial
    from quiverlab.fields import GF
    from quiverlab.modules.ar import end_action_on_ext1
    A = truncated_polynomial(4, field=GF(32003))
    M = A.module({1: 2}, {"x": [[0, 0], [1, 0]]}, name="M")   # k[x]/(x^2)
    M.check_module()
    assert M.is_indecomposable()
    tauM = M.tau()
    _basis, action = end_action_on_ext1(M, tauM)
    assert len(action[0]) == 2                                 # dim Ext^1(M, tau M) = 2
    ses = M.almost_split_sequence()
    mid_dims = sorted(s.dim for s, mult in ses.M.decompose() for _ in range(mult))
    assert mid_dims == [1, 3]                                  # true AR mesh middle
    assert mid_dims != [4]                                     # NOT the projective middle


@lit
def test_ka3_middle_terms_match_the_mesh():
    # kA3 (1->2->3, linear) AR quiver is the ZA3 slice. The almost-split sequence ending
    # at the simple S2 = [2] is the mesh   0 -> [3] -> [2,3] -> [2] -> 0:
    #   tau(S2) = [3] = S3 = P3, middle E = [2,3] = P2 (indecomposable), dim 2.
    A = _a3()
    M = A.simple(2)                                 # interval [2]
    ses = M.almost_split_sequence()
    assert ses.L.dimension_vector() == {1: 0, 2: 0, 3: 1}          # tau S2 = [3]
    assert ses.M.dimension_vector() == {1: 0, 2: 1, 3: 1}          # E = [2,3]
    assert is_isomorphic(ses.M, A.projective(2))                   # E = P2, indecomposable
