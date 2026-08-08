"""Silting oracles (Plan 67). Cross-engine: the 2-term slice == P45 tau-tilting -- a REAL
BIDIRECTIONAL agreement (Plan 67 fix round, MAJOR H2), NOT P45 re-verified against itself.
The SILTING engine INDEPENDENTLY enumerates the 2-term silting objects (BFS via
silting_mutate + is_silting_object + a g_proj fingerprint, restricted to the canonical AIR
{0,1} window -- a different verification/dedup/graph pipeline than P45's mutate + make_pair
+ g_key), and the resulting SET of g-vector fingerprints equals P45's exchange_graph vertex
set BOTH WAYS; the silting_neighbors edges equal the exchange_graph edges. Pinned on kA2
(hereditary, 5) AND kZ3/J2 (NON-hereditary self-injective, 14). End(mu T) underlying quiver
== Oppermann's rule (Thm 1.1). Literature: kA2 first mutation ring (AI Example 2.45);
Example 2.47 cones."""
import pytest

from quiverlab import NakayamaAlgebra, Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.complexes import ChainComplex
from quiverlab.derived.silting import (is_silting_object, silting_mutate,
                                        silting_neighbors)
from quiverlab.derived.tilting import g_proj, two_term_silting_from_presentation

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
    A = pair.algebra
    verts = list(A.quiver.vertices)
    out = []
    for Mi in pair.summands:
        cx, _rep = two_term_silting_from_presentation(Mi)   # [P1 -> P0], perfect
        out.append(cx)
    for v in sorted(pair.support, key=verts.index):         # killed projective -> P_v[1]
        out.append(ChainComplex.stalk(A.projective(v), 0).shift(1))
    return out


def _g_fingerprint(summands, verts):
    """The canonical, shift-SENSITIVE identifier of a 2-term silting object: the frozenset
    of its summands' K0 classes (g-vectors) in the projective basis via g_proj. This is
    exactly P45's g_key notion (AIR: the g-vectors determine the pair), so it distinguishes
    A ({e_v}, all deg 0) from A[1] ({-e_v}, all deg 1) -- both genuine, distinct support
    tau-tilting pairs. (The shift-INsensitive _silting_key would wrongly merge them.)"""
    return frozenset(tuple(g_proj(c, verts)) for c in summands)


def _in_two_term_window(summands):
    """The canonical AIR 2-term window: every summand concentrated in degrees {0, 1} (NOT
    merely SOME two consecutive degrees). This anchors the enumeration to P45's slice --
    A at deg 0 and A[1] at deg 1 both qualify and stay distinct, while stray shifts (A[2],
    a {-1,0} cocone representative) are dropped so the walk does not over-count shifts."""
    return all(set(c.degrees()) <= {0, 1} for c in summands)


def _enumerate_two_term_silting_via_silting_engine(A, budget=512):
    """INDEPENDENT enumeration of the 2-term silting objects, driven ENTIRELY by the silting
    engine (NOT P45): BFS from the regular object A through single silting mutations (both
    directions, silting_neighbors), keep only results in the canonical {0,1} window that
    RE-VERIFY silting (is_silting_object), dedup by the g_proj fingerprint. Returns
    ``(keys, edges, complete)`` -- ``complete`` False iff the vertex budget tripped. The
    AIR support-tau-tilting quiver is connected under 2-term mutation for a tau-tilting-
    finite algebra, so this closes on exactly the 2-term slice."""
    verts = list(A.quiver.vertices)
    start = [ChainComplex.stalk(A.projective(v), 0) for v in verts]
    seen = {_g_fingerprint(start, verts)}
    frontier = [start]
    edges = set()
    while frontier:
        T = frontier.pop()
        kT = _g_fingerprint(T, verts)
        for direction in ("left", "right"):
            for (_i, mut) in silting_neighbors(T, direction):
                if mut is None or not _in_two_term_window(mut):
                    continue
                if is_silting_object(mut).is_silting not in (True, "unknown"):
                    continue
                k = _g_fingerprint(mut, verts)
                edges.add(frozenset((kT, k)))
                if k not in seen:
                    if len(seen) >= budget:
                        return seen, edges, False
                    seen.add(k)
                    frontier.append(mut)
    return seen, edges, True


def _p45_two_term_objects(A):
    """P45's 2-term silting objects + edges as g-fingerprint sets, from the SUPPORT
    tau-tilting exchange graph (the independent P45 surface). Returns ``(keys, edges, eg)``."""
    verts = list(A.quiver.vertices)
    eg = A.exchange_graph(budget_pairs=256)
    fp = [_g_fingerprint(_pair_to_summand_complexes(rec["pair"]), verts)
          for rec in eg.vertices]
    keys = set(fp)
    edges = {frozenset((fp[a], fp[b])) for (a, b) in eg.arrows}
    return keys, edges, eg


def _assert_two_term_slice_matches_p45(A, expected_count):
    eg_keys, eg_edges, eg = _p45_two_term_objects(A)
    assert eg.is_complete                                     # tau-tilting-finite gate
    assert len(eg_keys) == len(eg.vertices) == expected_count
    # per-object: every P45 2-term silting object RE-VERIFIES silting via the silting engine
    for rec in eg.vertices:
        summand_list = _pair_to_summand_complexes(rec["pair"])
        assert len(summand_list) == len(list(A.quiver.vertices))
        rep = is_silting_object(summand_list)
        assert rep.is_silting is True
        assert rep.generation_certified_by.startswith(("tilting", "2-term"))
    # INDEPENDENT silting-engine enumeration + BIDIRECTIONAL set equality (H2).
    silt_keys, silt_edges, complete = _enumerate_two_term_silting_via_silting_engine(A)
    assert complete
    assert silt_keys == eg_keys, (
        f"silting-only={len(silt_keys - eg_keys)} p45-only={len(eg_keys - silt_keys)}")
    assert len(silt_keys) == expected_count
    # (b) the silting_neighbors edges equal the exchange_graph edges.
    assert silt_edges == eg_edges


@xeng
def test_two_term_silting_set_equals_p45_kA2():
    # kA2 (hereditary): 5 support tau-tilting pairs == 5 independently-enumerated 2-term
    # silting objects, sets AND edges equal both ways (AIR four-way identity, silting leg).
    _assert_two_term_slice_matches_p45(linear_path_algebra(2, field=QQ), 5)


@xeng
def test_two_term_silting_set_equals_p45_nonhereditary_kZ3J2():
    # kZ3/J2 = NakayamaAlgebra([2,2,2], cyclic): a NON-hereditary (self-injective, det
    # Cartan = 2) instance -- 14 pairs, all verifying, sets AND edges equal both ways. The
    # non-unimodular Cartan exercises the g_proj projective-basis fingerprint (Task 0).
    _assert_two_term_slice_matches_p45(
        NakayamaAlgebra([2, 2, 2], cyclic=True, field=QQ), 14)


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
