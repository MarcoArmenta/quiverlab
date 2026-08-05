"""AR-quiver knitting. Literature ground truth: kA_n has n(n+1)/2 indecomposables
(A2->3, A3->6, A4->10), D4 has 12; Nakayama is serial. Contract: a wild algebra
trips the budget LOUDLY (status 'budget', not a silent partial quiver)."""
import pytest

from quiverlab import (NakayamaAlgebra, Quiver, linear_path_algebra,
                       truncated_polynomial)
from quiverlab.fields import QQ
from quiverlab.families.dynkin import dynkin_quiver
from quiverlab.modules.ar import knit_ar_quiver

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@lit
@pytest.mark.parametrize("n, count", [(2, 3), (3, 6), (4, 10)])
def test_linear_a_n_indecomposable_count(n, count):
    A = linear_path_algebra(n, field=QQ)
    ar = knit_ar_quiver(A)
    assert ar.is_complete
    assert len(ar.vertices) == count


@lit
def test_d4_has_twelve_indecomposables():
    A = dynkin_quiver("D4").algebra(relations=[], field=QQ)
    ar = knit_ar_quiver(A)
    assert ar.is_complete and len(ar.vertices) == 12


@lit
def test_nakayama_serial_ar_quiver():
    # the linear Nakayama algebra kA3/rad^2 (n=3, l=2): its indecomposables are the
    # length-<=2 interval modules -- 3 simples + 2 length-2 = 5 (Kupisch series).
    A = NakayamaAlgebra(n=3, l=2, cyclic=False, field=QQ)
    ar = knit_ar_quiver(A)
    assert ar.is_complete
    assert len(ar.vertices) == 5


@selfcert
def test_mesh_relations_hold_on_a3():
    # every NON-projective vertex X (tau X != 0) is a mesh sink: at least one
    # irreducible map lands in it. (A projective source -- including the simple-
    # projective S_3 = P_3 -- may have no incoming arrows.)
    A = linear_path_algebra(3, field=QQ)
    ar = knit_ar_quiver(A)
    for j, vj in enumerate(ar.vertices):
        if vj["module"].tau().dim == 0:            # projective vertex: a source
            continue
        into = sum(m for (i, k), m in ar.arrows.items() if k == j)
        assert into >= 1                           # every non-projective is a mesh sink


@selfcert
@pytest.mark.parametrize("A", [
    truncated_polynomial(3, field=QQ),                 # k[x]/(x^3): self-injective, 3 indec
    NakayamaAlgebra(n=3, l=2, cyclic=True, field=QQ),  # kZ_3/rad^2: self-injective, 6 indec
])
def test_self_injective_knitting_refuses_loudly(A):
    # Devil's-advocate round (Finding A): the projective-seeded BFS DRAINS immediately on a
    # SELF-INJECTIVE algebra -- every indecomposable projective is also injective, so
    # tau^-(P) = 0 for each seed and the queue empties before discovering the rest of the
    # (periodic) component. It used to return status="complete" while grossly UNDERCOUNTING
    # (k[x]/(x^3): 1 vertex vs the true 3; cyclic kZ_3/rad^2: 3 vs 6). The projective-seeded
    # knitter must REFUSE LOUDLY, never emit a silently truncated "AR quiver".
    ar = knit_ar_quiver(A)
    assert ar.status == "unsupported"
    assert ar.is_complete is False
    assert ar.status != "complete"          # "complete" is UNREACHABLE for self-injective A
    assert "self-injective" in (ar.note or "")


@selfcert
def test_wild_algebra_trips_budget_loudly():
    # 3-Kronecker (wild): knitting cannot close -- must stop with status 'budget'.
    # A small budget_dim trips the loud cap fast (the preprojectives explode in dim).
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2), "c": (1, 2)}).algebra(
        relations=[], field=QQ)
    ar = knit_ar_quiver(A, budget_modules=40, budget_dim=20)
    assert ar.is_complete is False and ar.status == "budget"
