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


# --------------------------------------------------------------------------- #
# Task B3: Dynkin enumeration + braid-orbit transitivity + closed-form count
# --------------------------------------------------------------------------- #
@lit
@pytest.mark.parametrize("n,expected", [(2, 3), (3, 16), (4, 125), (5, 1296)])  # (n+1)^(n-1)
def test_dynkin_An_ces_count(n, expected):
    # counts live-verified (QQ): A2/3/4/5 = 3/16/125/1296; A5 count leg ~ a few seconds.
    A = linear_path_algebra(n, field=QQ)
    rep = A.exceptional_sequences()
    assert rep.is_complete and rep.count == expected
    assert rep.closed_form_count == expected          # n! h^n / |W|
    assert rep.dynkin_type == f"A_{n}"
    if n <= 4:
        assert rep.transitive is True                 # braid-orbit BFS reached them all


@lit
@pytest.mark.slow   # opt-in: A5 braid-orbit transitivity certificate (heavier leg)
def test_dynkin_A5_transitive_optin():
    A = linear_path_algebra(5, field=QQ)
    rep = A.exceptional_sequences(budget=100_000, transitive=True)
    assert rep.count == 1296 and rep.closed_form_count == 1296
    assert rep.transitive is True                      # opt-in via transitive=True


@lit
def test_dynkin_D4_ces_count():
    A = Quiver([1, 2, 3, 4], {"a": (2, 1), "b": (3, 1), "c": (4, 1)}).algebra(field=QQ)  # subspace D4
    rep = A.exceptional_sequences()
    assert rep.is_complete and rep.count == 162 and rep.closed_form_count == 162
    assert rep.transitive is True
    assert rep.dynkin_type == "D_4"


@selfcert
def test_exceptional_count_equals_num_indecomposables_dynkin():
    """#exceptional = #indecomposable on Dynkin (every Dynkin indec is a rigid brick)."""
    for n in (2, 3, 4):
        A = linear_path_algebra(n, field=QQ)
        indecs = A.ar_quiver().vertices
        exc = [m for m in (v["module"] for v in indecs)
               if is_exceptional_module(A, m)]
        assert len(exc) == len(indecs)


@selfcert
def test_rep_infinite_hereditary_refused():
    Kron = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(field=QQ)   # 2-Kronecker (tame)
    rep = Kron.exceptional_sequences(budget=200)
    assert rep.is_complete is False and rep.status in ("budget", "unsupported")


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


# --------------------------------------------------------------------------- #
# Adversarial-review MINOR (a): c_matrix rows subset of the P45 wall normals.
# --------------------------------------------------------------------------- #
@xeng
@pytest.mark.parametrize("n", [2, 3])
def test_c_matrix_rows_subset_of_p45_wall_normals(n):
    """The c-matrix rows (the dim-vectors of a complete exceptional sequence's terms) are a
    SUBSET of the P45 tau-tilting WALL NORMALS -- the brick dim-vectors labelling the
    exchange-graph edges. This is the plan-promised crosscheck: the hereditary
    exceptional-sequence reading of the c-vectors overlaps the P45 wall-and-chamber surface
    (Plan 65 R28; the P45 exchange graph is on this branch, no P63 code needed)."""
    from quiverlab.tautilting.mutation import exchange_graph
    A = linear_path_algebra(n, field=QQ)
    verts = list(A.quiver.vertices)
    eg = exchange_graph(A, budget_pairs=2000)
    assert eg.status == "complete"
    walls = {tuple(int(lab["brick"].get(v, 0)) for v in verts)
             for lab in eg.arrows.values() if lab.get("brick") is not None}
    rep = A.exceptional_sequences()
    rows = {tuple(row) for seq in rep.sequences for row in c_matrix(A, seq)}
    assert rows                                          # non-empty (rep-finite hereditary)
    assert rows <= walls                                 # every c-matrix row is a wall normal
