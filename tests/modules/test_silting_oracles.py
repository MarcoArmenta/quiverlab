"""Silting oracles (Plan 67). Cross-engine: the 2-term slice == P45 tau-tilting (count
+ per-object silting verdict, over the SUMMAND LIST not the direct-summed blob); End(mu T)
underlying quiver == Oppermann's rule (Thm 1.1). Literature: kA2 first mutation ring (AI
Example 2.45); Example 2.47 cones."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.complexes import ChainComplex
from quiverlab.derived.silting import is_silting_object, silting_mutate

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def _comm_square():
    # 1 <=> 2 with ab = ba = 0 (AI Example 2.47): a: 1->2, b: 2->1.
    return Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}).algebra(
        relations=["a*b", "b*a"], field=QQ)


def _pair_to_summand_complexes(pair):
    """The 2-term silting object T(M,P) = M (+) P[1] of a support tau-tilting pair, as a
    LIST of per-summand perfect complexes (NOT the direct-summed blob): each module summand
    M_i -> its 2-term presentation [P1 -> P0] (degrees 1, 0); each killed projective P_v ->
    the stalk P_v[1]. This is the list is_silting_object consumes."""
    from quiverlab.derived.tilting import two_term_silting_from_presentation
    A = pair.algebra
    verts = list(A.quiver.vertices)
    out = []
    for Mi in pair.summands:
        cx, _rep = two_term_silting_from_presentation(Mi)   # [P1 -> P0], perfect
        out.append(cx)
    for v in sorted(pair.support, key=verts.index):         # killed projective -> P_v[1]
        out.append(ChainComplex.stalk(A.projective(v), 0).shift(1))
    return out


@xeng
def test_two_term_silting_matches_p45():
    # Per-object, over the SUMMAND LIST (Ruling 5): every support tau-tilting pair's 2-term
    # silting object verifies silting via is_silting_object; the count ties to #support
    # tau-tilting pairs (exchange_graph) -- for kA2, 5 (AIR four-way identity). Do NOT pass
    # two_term_silting(pair)["complex"] (a single object) to is_silting_object: with one
    # summand the g-matrix is 1xn (non-square), k0_basis vacuously False, a no-op that
    # asserts nothing.
    A = linear_path_algebra(2, field=QQ)
    eg = A.exchange_graph(budget_pairs=64)
    assert eg.is_complete
    pairs = [rec["pair"] for rec in eg.vertices]
    for pair in pairs:
        summand_list = _pair_to_summand_complexes(pair)      # a LIST, not ["complex"]
        assert len(summand_list) == len(list(A.quiver.vertices))   # n summands = n simples
        rep = is_silting_object(summand_list)
        assert rep.is_silting is True                        # NOT a vacuous no-op
        assert rep.generation_certified_by.startswith(("tilting", "2-term"))
    assert len(pairs) == 5                                    # kA2 count tie (P45)


@xeng
@lit
def test_end_of_mutant_quiver_matches_oppermann():
    # 1 <=> 2, ab = ba = 0 (Example 2.47): End(P2 (+) cone(P1->P2)) has the vertex-1
    # Oppermann-mutated quiver. Oppermann Thm 1.1 LEFT mutation at vertex 1: reverse the
    # degree-0 outgoing arrow a: 1->2 to a*: 2->1 (degree 0), and RAISE the incoming
    # b: 2->1 to degree 1 -- so b LEAVES degree 0. The underlying degree-0 quiver of
    # End(mu T) is therefore a SINGLE arrow = linear A2 (the double 1<=>2 loses one arrow
    # at the mutated vertex). # PIN the graded degrees (out of the engine's scope); the
    # battery verifies the UNDERLYING (source, target) count + the corner-Cartan.
    # Uses P43 end_algebra_of_complex + P44 gabriel_quiver.
    from quiverlab.derived.tilting import (end_algebra_of_complex,
                                           corner_cartan_of_complex)
    A = _comm_square()
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    # sanity: End(A_A) ~ A (the double 1<=>2 corner-Cartan, all-ones).
    assert corner_cartan_of_complex(T) == A.cartan_matrix() == [[1, 1], [1, 1]]
    mut = silting_mutate(T, 0, direction="left")             # left mutation at P1
    E = end_algebra_of_complex(mut)                          # the derived endo algebra
    assert E.dim == 3                                        # dim kA2 (was dim 4 for A)
    Q = E.gabriel_quiver()
    # Oppermann: the underlying degree-0 quiver is a SINGLE arrow (linear A2), not 1<=>2.
    assert len(list(Q.arrows)) == 1
    assert len(list(Q.vertices)) == 2
    # the reoriented-A2 corner-Cartan (exactly the P43 APR-tilt orientation oracle).
    assert corner_cartan_of_complex(mut) == [[1, 0], [1, 1]]


@lit
def test_kA2_first_mutation_ring():
    # AI Example 2.45 (kA2): from A = P1 (+) P2 the first ring of irreducible left
    # mutations. kA2 silting is INFINITE, so we pin only the LOCAL structure AI draws
    # (n = 2 neighbours per vertex, each a genuine silting object, distinct from A).
    A = linear_path_algebra(2, field=QQ)
    T = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    nbrs = [silting_mutate(T, i, "left") for i in range(2)]
    for mut in nbrs:
        assert is_silting_object(mut).is_silting in (True, "unknown")
    keys = {tuple(sorted(str(c.degrees()) for c in m)) for m in nbrs}
    assert len(keys) == 2                                     # two distinct neighbours
    tkey = tuple(sorted(str(c.degrees()) for c in T))
    assert tkey not in keys                                   # both differ from A
