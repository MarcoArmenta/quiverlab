"""Classical hereditary exceptional sequences (Plan 65 / R28; Crawley-Boevey Ottawa 1992;
Ringel Contemp. Math. 171 (1994)). Recognizer (Hom/Ext orthogonality), braid mutation
sigma_i, Dynkin enumeration, closed-form counts. Hereditary + rep-finite scope; QQ (the
knit/decompose char caveat)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.families.basic import linear_path_algebra
from quiverlab.modules.exceptional import (c_matrix, is_exceptional_module,
                                           is_exceptional_sequence)

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


# --------------------------------------------------------------------------- #
# Task B1: recognizers + c-matrix
# --------------------------------------------------------------------------- #
@selfcert
def test_all_dynkin_indecs_are_exceptional():
    A = linear_path_algebra(3, field=QQ)
    for v in A.ar_quiver().vertices:
        assert is_exceptional_module(A, v["module"]) is True   # every A_n indec is a rigid brick


@lit
def test_a2_exceptional_pair_orthogonality():
    A = linear_path_algebra(2, field=QQ)                       # 1->2 : S1, S2, P1
    S1, S2 = A.simple(1), A.simple(2)
    P1 = A.projective(1)
    # CONVENTION: for i<j need Hom(Ej,Ei)=0 AND Ext^1(Ej,Ei)=0 (no BACKWARD maps).
    # Live-verified (QQ): Hom(S2,P1)=1 (socle inclusion S2=soc P1 -> P1), Hom(P1,S2)=0,
    # Ext^1(P1,S2)=0. So (S2, P1) IS exceptional (no backward Hom/Ext), while (P1, S2) is NOT.
    assert is_exceptional_sequence(A, [S2, P1]) is True
    assert is_exceptional_sequence(A, [P1, S2]) is False
    # a complete CES has length n=2; the three A2 CES are (S1,S2),(S2,P1),(P1,S1) -- Task B3.


@selfcert
def test_single_exceptional_module_is_a_sequence():
    A = linear_path_algebra(2, field=QQ)
    assert is_exceptional_sequence(A, [A.simple(1)]) is True
    assert is_exceptional_module(A, A.projective(1)) is True


@selfcert
def test_c_matrix_rows_are_dimvectors():
    A = linear_path_algebra(2, field=QQ)                       # vertex order [1, 2]
    S1, P1 = A.simple(1), A.projective(1)
    assert c_matrix(A, [P1, S1]) == [[1, 1], [1, 0]]


@selfcert
def test_non_hereditary_refused():
    from quiverlab.errors import QuiverlabError
    from quiverlab.families.radical_square_zero import RadicalSquareZero
    A = RadicalSquareZero(Quiver([1, 2, 3], {"a": (2, 1), "b": (3, 2)}), field=QQ)  # kA3/rad^2
    with pytest.raises(QuiverlabError):
        is_exceptional_sequence(A, [A.simple(1)])
    with pytest.raises(QuiverlabError):
        is_exceptional_module(A, A.simple(1))
    with pytest.raises(QuiverlabError):
        c_matrix(A, [A.simple(1)])
