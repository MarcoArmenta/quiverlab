"""The picture group presentation pi_1(|W(A)|) (Plan 66 / R29, Igusa-Todorov-Weyman 1609.02636;
Hanson-Igusa). Generators = bricks; relations = rank-2 wides (commutation for k x k, atom/
pentagon for connected). kA2 -> (3 gens, 1 atom); kA3 -> (6 gens, 4 atom + 2 comm). The
DISCRIMINATOR: kA3 and kZ3/rad^2 share every coarse count but split 4+2 vs 3+3. QQ-scope."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.cluster_morphism import picture_group
from quiverlab.tautilting.torsion import bricks

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _kZ3_rad2():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}
                  ).algebra(relations=["a*b", "b*c", "c*a"], field=QQ)

def _kA3_rad2():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=QQ)


@lit
def test_kA2_picture_group():
    G = picture_group(Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ))
    assert G.num_generators == 3 and G.num_relations == 1
    assert G.num_atom == 1 and G.num_commutation == 0          # single connected rank-2 wide
    assert G.abelianization_rank == 2                          # Z^{3-1}, 1 distinct ext brick


@lit
def test_kA3_picture_group():
    G = picture_group(linear_path_algebra(3, field=QQ))
    assert G.num_generators == 6 and G.num_relations == 6
    assert G.num_atom == 4 and G.num_commutation == 2          # 4 kA2-type + 2 kxk
    # H2 correction: NOT Z^{6-4}=Z^2. The 4 atoms share extension bricks -> only 3 DISTINCT
    # ext bricks ({S1,M23} and {M12,S3} both extend to the top brick M123), so rank = 3.
    assert G.abelianization_rank == 3                          # Z^{6 - rank(rel matrix)} = Z^{6-3}


@lit
def test_kA4_picture_group_abelianization():
    # H2 counterexample the draft LACKED: #bricks-#atom = 10-10 = 0 would give the ABSURD Z^0.
    # The 10 atoms have only 6 DISTINCT extension bricks, so rank(rel matrix) = 6 and G^ab = Z^4.
    G = picture_group(linear_path_algebra(4, field=QQ))
    assert G.num_generators == 10 and (G.num_atom, G.num_commutation) == (10, 10)
    assert G.abelianization_rank == 4                          # Z^{10-6}, NOT Z^{10-10}=Z^0


@xeng
@pytest.mark.parametrize("A_factory, ngen", [
    (lambda: linear_path_algebra(2, field=QQ), 3),
    (lambda: linear_path_algebra(3, field=QQ), 6),
    (lambda: linear_path_algebra(4, field=QQ), 10),
    (_kA3_rad2, 5), (_kZ3_rad2, 6)])
def test_generators_are_bricks(A_factory, ngen):
    A = A_factory()
    assert picture_group(A).num_generators == len(bricks(A)) == ngen


@lit
def test_kA3_vs_kZ3rad2_relation_split_discriminates():
    # THE anti-tautology pin: identical (#bricks, #wide, g_fan_face_vector) but DIFFERENT relation
    # split -- kA3 = 4 atom + 2 comm, kZ3/rad^2 = 3 atom + 3 comm. A construction that ignores
    # the Ext^1 direction/existence between the two simple bricks would give the same split (and
    # the same classifying-space face vector -- which actually differs 49 vs 48, see Task A).
    kA3 = picture_group(linear_path_algebra(3, field=QQ))
    kZ3 = picture_group(_kZ3_rad2())
    assert (kA3.num_generators, kA3.num_relations) == (kZ3.num_generators, kZ3.num_relations) == (6, 6)
    assert (kA3.num_atom, kA3.num_commutation) == (4, 2)
    assert (kZ3.num_atom, kZ3.num_commutation) == (3, 3)       # the discriminator


@xeng
def test_relation_count_equals_rank2_wide_count():
    # #relations == #(rank-2 wides) == #(size-2 all-semibricks). Cross-check via the category's
    # rank-2 objects (Task A) -- the two enumerations must agree.
    from quiverlab.tautilting.cluster_morphism import tau_cluster_category
    A = linear_path_algebra(3, field=QQ)
    G = picture_group(A); C = tau_cluster_category(A)
    assert G.num_relations == sum(1 for o in C.objects if o[1] == 2) == 6


@selfcert
def test_atom_relation_has_extension_brick():
    # Each atom relation names an EXTENSION brick (the third brick of the connected wide);
    # each commutation relation has ext_brick_id None. Structural, catches type mislabelling.
    G = picture_group(linear_path_algebra(3, field=QQ))
    atoms = [r for r in G.relations if r[0] == "atom"]
    comms = [r for r in G.relations if r[0] == "commutation"]
    assert len(atoms) == 4 and all(r[2] is not None for r in atoms)
    assert len(comms) == 2 and all(r[2] is None for r in comms)


@lit
@pytest.mark.parametrize("A_factory, kpi1", [
    (_kA3_rad2, True), (_kZ3_rad2, True),                        # Nakayama -> K(pi,1) (Hanson-Igusa)
    (lambda: linear_path_algebra(3, field=QQ), True)])           # hereditary Dynkin -> K(pi,1)
def test_kpi1_theorem_anchored(A_factory, kpi1):
    assert picture_group(A_factory()).is_kpi1 is kpi1


@selfcert
def test_kpi1_not_claimed_beyond_theorems():
    # A tau-tilting-finite, NON-Nakayama, NON-hereditary algebra: is_kpi1 must be None (honest
    # "not certified"), NEVER True -- P66 never asserts asphericity beyond the cited theorems.
    # (A commutative-square-with-relation kQ/I: tau-tilting-finite, not Nakayama, not hereditary.)
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (1, 3), "c": (2, 4), "d": (3, 4)})
    A = Q.algebra(relations=["a*c-b*d"], field=QQ)               # commutative square
    G = picture_group(A)
    if G.is_complete:
        assert G.is_kpi1 is None and "not certified" in G.kpi1_reason
