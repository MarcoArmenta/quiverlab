"""Basic-ization + Gabriel-quiver recovery of a structure-constant algebra (Plan 44 / C7).
The Wedderburn crux: primitive idempotents via the trace-form radical + a semisimple
quotient + Newton lifting; the presentation via a length-lex kernel enumeration. Over QQ /
large-prime GF only (the trace-form radical needs char 0 or char > dim)."""
import pytest

from quiverlab import GF, Quiver
from quiverlab.core.algebra import Algebra
from quiverlab.core.basic import (basic_algebra, gabriel_quiver,
                                   idempotent_classes, presented_form,
                                   primitive_idempotents)
from quiverlab.fields import QQ


def _M2(field):
    """M_2(k) as structure constants in the matrix-unit basis (E11,E12,E21,E22):
    E_ij E_kl = delta_jk E_il ; unit = E11 + E22."""
    z = "0"
    def e(i):  # basis vector
        v = [z, z, z, z]; v[i] = "1"; return v
    # index: 0=E11,1=E12,2=E21,3=E22 ; pair (row,col): 0=(1,1),1=(1,2),2=(2,1),3=(2,2)
    rc = {0: (1, 1), 1: (1, 2), 2: (2, 1), 3: (2, 2)}
    idx = {v: k for k, v in rc.items()}
    T = [[["0"] * 4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        i, j = rc[a]
        for b in range(4):
            k, l = rc[b]
            if j == k:
                T[a][b][idx[(i, l)]] = "1"
    return Algebra.from_structure_constants(T, ["1", "0", "0", "1"], field=field)


@pytest.mark.oracle_literature
def test_matrix_algebra_basic_is_the_field():
    A = _M2(GF(32003))                       # char 32003 > dim 4
    idems = primitive_idempotents(A)
    assert len(idems) == 2                    # E11, E22 (both rank-1)
    classes = idempotent_classes(A)
    assert len(classes) == 1                  # conjugate: one iso class
    B = basic_algebra(A)
    assert B.dim == 1                          # e M_2 e = k
    Q = gabriel_quiver(A)
    assert len(list(Q.vertices)) == 1 and len(Q.arrows) == 0
    assert presented_form(A).dim == 1


@pytest.mark.oracle_selfcert
def test_primitive_idempotents_are_a_complete_orthogonal_set():
    A = _M2(QQ)
    idems = primitive_idempotents(A)
    dom = A.domain
    # orthogonal: e_i e_j = delta_ij e_i ; complete: sum e_i = unit.
    total = [dom.zero()] * A.dim
    for i, ei in enumerate(idems):
        assert A.multiply(ei, ei) == ei                       # idempotent
        for j, ej in enumerate(idems):
            prod = A.multiply(ei, ej)
            if i != j:
                assert all(dom.is_zero(x) for x in prod)      # orthogonal
        total = [dom.add(total[t], ei[t]) for t in range(A.dim)]
    assert total == list(A.unit)                              # complete


@pytest.mark.oracle_selfcert
def test_kA2_round_trip_identity():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)   # basic already
    B = presented_form(A)
    assert B.dim == A.dim == 3
    assert B.cartan_matrix() == A.cartan_matrix()             # recovered kA2 == kA2


@pytest.mark.oracle_crossengine
def test_end_of_regular_recovers_kA3():
    # End_A(A_A) ~ A (P37 endomorphism.py). Take A = kA3 (hereditary 1->2->3), build the
    # regular module's End as a presentation-LESS structure-constant algebra, and recover
    # kA3 from it -- tying to P37's regular_corner_dims (Cartan of End = Cartan of A).
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=[], field=QQ)
    from quiverlab.modules.endomorphism import end_algebra, regular_corner_dims
    from quiverlab.modules.morphism import direct_sum
    reg, _, _ = direct_sum(*[A.projective(v) for v in (1, 2, 3)])
    E = end_algebra(reg)                                       # quiver is None (presentation-less)
    assert E.quiver is None
    R = presented_form(E)
    assert R.dim == A.dim                                      # dim kA3 = 6
    assert R.cartan_matrix() == A.cartan_matrix()             # recovered kA3
    assert regular_corner_dims(A) == [[int(x) for x in row]   # P37 sided oracle tie
                                      for row in A.cartan_matrix()]


def _perm_similar(C, D):
    """True iff some permutation P has ``P^T C P == D`` (relabeling the vertices). Cartan
    matrices of the SAME algebra under different vertex orders are permutation-similar, not
    equal -- the cyclic case genuinely permutes (unlike the linear kA3 recovery, which is
    exactly equal)."""
    import itertools
    n = len(C)
    for perm in itertools.permutations(range(n)):
        # D[i][j] ?= C[perm[i]][perm[j]]  (i.e. P sends new i -> old perm[i])
        if all(D[i][j] == C[perm[i]][perm[j]] for i in range(n) for j in range(n)):
            return True
    return False


@pytest.mark.oracle_crossengine
def test_cyclic_presented_form_round_trip():
    # kZ_3 / J^2 (cyclic Nakayama, dim 6): rebuild it as a presentation-LESS structure-
    # constant algebra (End of its regular module, End(A_A) ~ A), then recover kQ/I by
    # presented_form. The cyclic symmetry permutes the vertex labels, so the recovered
    # Cartan is PERMUTATION-SIMILAR to the original (not equal, unlike the linear kA3 case).
    from quiverlab import NakayamaAlgebra
    from quiverlab.modules.endomorphism import end_algebra
    from quiverlab.modules.morphism import direct_sum
    N = NakayamaAlgebra(n=3, l=2, cyclic=True, field=QQ)      # kZ_3/J^2
    reg, _, _ = direct_sum(*[N.projective(v) for v in N.quiver.vertices])
    E = end_algebra(reg)                                      # presentation-less rebuild
    assert E.quiver is None and E.dim == N.dim
    R = presented_form(E)
    assert R.dim == N.dim                                     # dim recovered exactly (6)
    assert _perm_similar(N.cartan_matrix(), R.cartan_matrix())


@pytest.mark.oracle_selfcert
def test_char_caveat_refuses_loudly():
    from quiverlab.errors import QuiverlabError
    A = _M2(GF(2))                            # char 2 <= dim 4: trace-form radical unreliable
    with pytest.raises(QuiverlabError):
        primitive_idempotents(A)
