"""The tau-cluster morphism category W(A) (Plan 66 / R29, Buan-Marsh IMRN 2021 + Hanson-Igusa
Comm. Alg. 2021). Objects = tau-perpendicular wide subcategories (count = #wide, ties P64);
morphisms = support tau-rigid pairs of the source (out-degree = closed-star size = #faces(C_W));
the Hanson-Igusa CLASSIFYING-SPACE face vector f_0 = #wide, f_k = #(rank-k morphisms), chi =
sum (-1)^k f_k (NOT the g-fan SPHERE, reported separately as g_fan_face_vector with f_n = #sTt);
K(pi,1) theorem-anchored (Nakayama / hereditary Dynkin). QQ-scope (the P45 char caveat)."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.cluster_morphism import tau_cluster_category
from quiverlab.tautilting.congruence import wide_subcategories
from quiverlab.tautilting.torsion import bricks

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kA2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ)

def _nakayama_kA3_rad2():          # linear, non-hereditary, "beyond kA_n"
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=QQ)

def _nakayama_kZ3_rad2():          # cyclic, self-injective
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}
                  ).algebra(relations=["a*b", "b*c", "c*a"], field=QQ)


@xeng
@pytest.mark.parametrize("A_factory, nwide", [
    (lambda: linear_path_algebra(2, field=QQ), 5),
    (lambda: linear_path_algebra(3, field=QQ), 14),
    (_nakayama_kA3_rad2, 12),
    (_nakayama_kZ3_rad2, 14),
])
def test_object_count_equals_wide_count(A_factory, nwide):
    # OBJECTS = #wide via TWO INDEPENDENT ROUTES: P66's tau-perpendicular enumeration (each face
    # U -> J(U), deduped by iso-class simple-brick semibrick) vs P64's Enomoto core-label-order
    # size. NOT a tautology. Also NOT #semibricks/#torsion -- those coincide only rep-finite.
    A = A_factory()
    C = tau_cluster_category(A)
    assert C.is_complete and C.object_count == nwide
    assert C.object_count == wide_subcategories(A).size        # the cross-tie (ruling: .size)


@lit
@pytest.mark.parametrize("n, facevec, gfan", [
    (2, (5, 11, 5),          (1, 5, 5)),
    (3, (14, 49, 49, 14),    (1, 9, 21, 14)),
    (4, (42, 204, 326, 204, 42), (1, 14, 56, 84, 42))])
def test_face_vector_type_A(n, facevec, gfan):
    # face_vector = the Hanson-Igusa CLASSIFYING-SPACE cube complex (f_0 = #wide, f_k =
    # #(rank-k morphisms)), NOT the g-fan sphere. g_fan_face_vector = the SPHERE (f_n = #sTt).
    C = tau_cluster_category(linear_path_algebra(n, field=QQ))
    assert C.face_vector == facevec
    assert C.g_fan_face_vector == gfan
    assert C.face_vector[0] == C.object_count          # f_0 == #wide (discriminating self-cert)
    assert sum(C.face_vector) == C.morphism_count      # sum == total morphisms
    assert C.g_fan_face_vector[-1] == C.face_vector[-1]  # both = #sTt (top cells / top facets)


@selfcert
@pytest.mark.parametrize("n, chi", [(2, -1), (3, 0), (4, 2)])   # classifying-space chi (NOT (-1)^n)
def test_euler_characteristic(n, chi):
    C = tau_cluster_category(linear_path_algebra(n, field=QQ))
    assert C.euler_characteristic == chi
    assert C.euler_characteristic == sum((-1) ** k * f for k, f in enumerate(C.face_vector))
    # NB: chi is NOT (-1)^n -- that was the DISCARDED g-fan sphere's tautology (ruling H1).
    # For kA2 the classifying-space chi = -1 also equals the presentation Euler char 1 - 3 + 1.


@xeng
def test_morphisms_out_of_top_equals_g_fan_faces():
    # morphisms out of the terminal object mod A (U = empty) == #(all g-fan faces) ==
    # #(support tau-rigid pairs of A) == sum(g_fan_face_vector). Live: kA2->11, kA3->45.
    # NOTE (ruling H1): this is the SPHERE face count, NOT sum(face_vector) (= total morphisms).
    for A, top_out in [(_kA2(), 11), (linear_path_algebra(3, field=QQ), 45)]:
        C = tau_cluster_category(A)
        top = max(C.objects, key=lambda o: o[1])       # the rank-n object = mod A
        assert C.out_degree[top[0]] == top_out
        assert sum(C.g_fan_face_vector) == top_out      # closed star of the empty face = all faces
        assert sum(C.face_vector) != top_out            # the classifying space counts MORE (21/126)


@lit
@pytest.mark.parametrize("A_factory, total", [
    (_kA2, 21), (lambda: linear_path_algebra(3, field=QQ), 126)])
def test_total_morphism_count(A_factory, total):
    # total morphisms = sum over objects of the closed-star size (= #faces(C_W)). Hand-derived
    # kA2->21, kA3->126 (Verified pins). The per-object out-degree cross-checks against P65's
    # C_W face count at implementation (the Jasso link identity).
    C = tau_cluster_category(A_factory())
    assert C.morphism_count == total
    assert C.morphism_count == sum(C.out_degree.values())


@selfcert
def test_out_degree_is_closed_star_and_matches_reduction():
    # For every object W = J(U): out_degree[W] == |closedStar(U)| (star in A's g-fan) AND ==
    # #faces(C_W) computed from P65's tau_perpendicular_reduction(A, U).reduced exchange graph.
    # This is the Jasso link identity -- a genuine cross-engine self-cert of the morphism count.
    A = linear_path_algebra(3, field=QQ)
    C = tau_cluster_category(A)
    # spot-check the rank-2 objects: 4 must have #faces(C_W)=11 (C_W ~ kA2), 2 must have 9 (kxk)
    rank2 = [o for o in C.objects if o[1] == 2]
    assert len(rank2) == 6
    stars = sorted(C.out_degree[o[0]] for o in rank2)
    assert stars == [9, 9, 11, 11, 11, 11]              # 2 commutation-wides (9) + 4 kA2-wides (11)


@selfcert                                              # fix-round ruling 3: NOT xeng -- definitional in P65
def test_factorization_count_ties_P65():
    # #(complete signed exceptional sequences) = n! * #sTt.  P65 DEFINES signed_count :=
    # n! * |exchange_graph.vertices|, and P66's #sTt = face_vector[-1] = g_fan_face_vector[-1] =
    # |eg.vertices| off the SAME graph -- so this equality is DEFINITIONALLY true, an internal
    # consistency gate (oracle_selfcert), NOT a cross-engine oracle. (If P65 ships MATERIALISED
    # signed sequences, upgrade to len(materialised) == n!*#sTt for a genuine xeng check.)
    from math import factorial
    from quiverlab.tautilting.exceptional import tau_exceptional_sequences
    A = linear_path_algebra(3, field=QQ)
    C = tau_cluster_category(A)
    stt = C.face_vector[-1]                             # #sTt = top-rank facets (= g_fan top cell)
    assert factorial(3) * stt == tau_exceptional_sequences(A).signed_count


@lit
def test_kA3_vs_kZ3rad2_share_g_fan_but_split_differs():
    # The anti-tautology pin (see also test_picture_group_split): kA3 and kZ3/rad^2 agree on every
    # COARSE invariant AND on the g-fan SPHERE, but differ on the Ext-driven relation split -- and,
    # consequently, on the CLASSIFYING-SPACE face vector (ruling H1: the sphere sees nothing, the
    # classifying space counts morphisms).
    kA3 = tau_cluster_category(linear_path_algebra(3, field=QQ))
    kZ3 = tau_cluster_category(_nakayama_kZ3_rad2())
    assert kA3.object_count == kZ3.object_count == 14
    assert kA3.g_fan_face_vector == kZ3.g_fan_face_vector == (1, 9, 21, 14)   # identical spheres
    assert len(bricks(linear_path_algebra(3, field=QQ))) == len(bricks(_nakayama_kZ3_rad2())) == 6
    # the classifying spaces DIFFER (49 vs 48): the bonus discriminator from the relation split
    assert kA3.face_vector == (14, 49, 49, 14)
    assert kZ3.face_vector == (14, 48, 48, 14)
    assert kA3.face_vector != kZ3.face_vector


@lit
def test_tau_tilting_infinite_refuses_no_partial_category():
    # 2-Kronecker: tau-tilting-INFINITE (DIJ) -> W(A) is infinite -> NO category emitted.
    K = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ)
    C = tau_cluster_category(K, budget=40)
    assert C.is_complete is False and C.status in ("budget", "error")
    assert C.objects == () and C.morphism_count is None and C.face_vector is None and C.note
