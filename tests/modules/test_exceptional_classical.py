"""Classical hereditary exceptional sequences (Plan 65 / R28; Crawley-Boevey Ottawa 1992;
Ringel Contemp. Math. 171 (1994)). Recognizer (Hom/Ext orthogonality), braid mutation
sigma_i, Dynkin enumeration, closed-form counts. Hereditary + rep-finite scope; QQ (the
knit/decompose char caveat)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.families.basic import linear_path_algebra
from quiverlab.modules.exceptional import (braid_mutation, c_matrix,
                                           is_exceptional_module,
                                           is_exceptional_sequence)
from quiverlab.modules.hom import is_isomorphic

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


# --------------------------------------------------------------------------- #
# Task B2: braid mutation sigma_i
# --------------------------------------------------------------------------- #
def _seq_iso(A, s1, s2):
    return len(s1) == len(s2) and all(is_isomorphic(a, b) for a, b in zip(s1, s2))


def _local_ces(A):
    """A self-contained direct-orthogonality enumeration of the complete exceptional
    sequences of a rep-finite hereditary A (so the B2 tests do not depend on the B3
    library enumerator)."""
    from quiverlab.modules.ext import ext
    from quiverlab.modules.hom import end_dim, hom_dim
    n = len(list(A.quiver.vertices))
    mods = [v["module"] for v in A.ar_quiver().vertices
            if end_dim(v["module"]) == 1 and ext(A, v["module"], v["module"], 1) == 0]

    def ok(seq):
        for i in range(len(seq)):
            for j in range(i + 1, len(seq)):
                if hom_dim(seq[j], seq[i]) != 0 or ext(A, seq[j], seq[i], 1) != 0:
                    return False
        return True

    out = []

    def dfs(cur):
        if len(cur) == n:
            out.append(list(cur))
            return
        for m in mods:
            cur.append(m)
            if ok(cur):
                dfs(cur)
            cur.pop()

    dfs([])
    return out


@lit
def test_worked_braid_mutation_a2_kA2():
    """H4 hand-worked pin: sigma_1(P_1, S_1) = (S_2, P_1) via 0 -> S_2 -> P_1 -> S_1 -> 0
    (evaluation kernel, the clean surjective case). L_{P_1} S_1 = rad P_1 = S_2 = (0,1)."""
    A = linear_path_algebra(2, field=QQ)
    P1, S1, S2 = A.projective(1), A.simple(1), A.simple(2)
    out = braid_mutation(A, [P1, S1], 0, direction="left")
    assert _seq_iso(A, out, [S2, P1])
    assert is_exceptional_sequence(A, out) is True


@selfcert
def test_braid_mutation_is_an_involution_a3():
    """sigma_i . sigma_i^{-1} = id (termwise is_isomorphic) and each output is a complete CES.
    Runs over a spanning set of A_3 exceptional pairs embedded in complete sequences."""
    A = linear_path_algebra(3, field=QQ)
    seqs = _local_ces(A)
    for s in seqs[:12]:
        for i in range(2):
            left = braid_mutation(A, s, i, direction="left")
            assert is_exceptional_sequence(A, left) is True
            back = braid_mutation(A, left, i, direction="right")
            assert _seq_iso(A, back, s)                 # sigma_i^{-1} sigma_i = id
            right = braid_mutation(A, s, i, direction="right")
            assert is_exceptional_sequence(A, right) is True
            fwd = braid_mutation(A, right, i, direction="left")
            assert _seq_iso(A, fwd, s)                  # sigma_i sigma_i^{-1} = id


@selfcert
def test_braid_mutation_output_module_exceptional_a3():
    A = linear_path_algebra(3, field=QQ)
    s = _local_ces(A)[0]
    out = braid_mutation(A, s, 0, direction="left")
    for E in out:
        assert is_exceptional_module(A, E) is True


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
