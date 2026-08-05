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
