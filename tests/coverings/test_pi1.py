"""pi1(Q,I): presentation + abelianization by exact SNF. Oracles: the CORRECTED
commutative square (WITH relation => 1, WITHOUT => Z), Le Meur Example 1, trees,
single loop, multi-loop free ranks, the annulus. Refs: Le Meur math/0503302;
Martinez-Villa-de la Pena 1983; Assem-de la Pena 1996."""
import pytest
from quiverlab import GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.errors import QuiverlabError
from quiverlab.invariants.coverings import fundamental_group, intrinsic_fundamental_group

pytestmark = pytest.mark.oracle_literature


def _square(rel):
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    return Q.algebra(relations=rel, field=GF(7))


def test_commutative_square_with_relation_is_trivial():
    A = _square(["a*b - c*d"])                 # commutativity: identifies the parallel paths
    g = fundamental_group(A)
    assert g.free_rank == 0 and g.invariant_factors == ()   # pi1^ab = 0
    assert g.abelianization_repr() == "0"


def test_commutative_square_without_relation_is_Z():
    A = _square([])                            # hereditary square: a hole survives
    g = fundamental_group(A)
    assert g.free_rank == 1 and g.invariant_factors == ()   # pi1^ab = Z
    assert g.abelianization_repr() == "Z"


def test_le_meur_example_1_monomial_vs_binomial():
    # Le Meur Ex.1 concretely + ADMISSIBLY (all relations in rad^2). d*a (1->2->3) and
    # d*c*b (1->2->4->3) are parallel 1->3 paths. I = <d*a> monomial (n=1) => NO
    # identification => pi1^ab = Z; J = <d*a - d*c*b> binomial minimal => identifies
    # the parallel paths => pi1^ab = 0.
    Q = Quiver([1, 2, 3, 4], {"d": (1, 2), "a": (2, 3), "c": (2, 4), "b": (4, 3)})
    Amono = Q.algebra(relations=["d*a"], field=GF(7))
    Abino = Q.algebra(relations=["d*a - d*c*b"], field=GF(7))
    assert fundamental_group(Amono).free_rank == 1        # I = <da>     -> Z
    assert fundamental_group(Abino).free_rank == 0        # J = <da-dcb> -> 0


def test_two_independent_pairs_not_over_glued():
    # DISCRIMINATING ORACLE (adversarial H1). Four parallel routes 1->3 and a
    # NON-minimal generating set I = <u1-u2, (u1-u2)+(u3-u4)>. The exact
    # I/(rad.I + I.rad) linear algebra recovers the TWO minimal relations u1-u2 and
    # u3-u4 => glue u1~u2 and u3~u4 ONLY => pi1^ab = Z (one route-class comparison
    # survives). The OLD greedy/minimality-refusal design would REFUSE g2 (proper
    # sub-relation g1 in I) or naively glue all four routes => pi1^ab = 0 (WRONG).
    Q = Quiver([1, 3, "A", "B", "C", "D"],
               {"a1": (1, "A"), "a2": ("A", 3), "b1": (1, "B"), "b2": ("B", 3),
                "c1": (1, "C"), "c2": ("C", 3), "d1": (1, "D"), "d2": ("D", 3)})
    A = Q.algebra(relations=["a1*a2 - b1*b2", "a1*a2 - b1*b2 + c1*c2 - d1*d2"],
                  field=GF(7))
    g = fundamental_group(A)
    assert g.free_rank == 1 and g.invariant_factors == ()     # Z, NOT 0


def test_tree_is_trivial():
    A = linear_path_algebra(5, field=GF(7))    # kA5: underlying graph a tree
    g = fundamental_group(A)
    assert g.free_rank == 0 and g.invariant_factors == () and g.components == 1


def test_single_loop_monomial_is_Z():
    A = truncated_polynomial(4, field=GF(7))   # k[x]/(x^4): one loop, monomial relation
    g = fundamental_group(A)
    assert g.free_rank == 1 and g.invariant_factors == ()   # presentation pi1 = Z (NOT ZxC_p)


def test_multi_loop_free_rank():
    A = Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "y*y", "x*y", "y*x"], field=GF(5))   # all monomial
    g = fundamental_group(A)
    assert g.free_rank == 2 and g.invariant_factors == ()      # free on 2 loops


def test_disconnected_is_direct_sum():
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (3, 4)}).algebra(field=GF(5))
    g = fundamental_group(A)
    assert g.components == 2 and g.free_rank == 0               # two trees


def test_torsion_pi1ab_triangular():
    # END-TO-END TORSION pin (adjudicated Fix 2): a TRIANGULAR (acyclic) presentation
    # whose pi1^ab has torsion. Two parallel arrows 1->2 and two parallel 2->3, with
    # the two admissible binomial relations a*c - b*d and a*d - b*c. The minimal
    # relations glue ac~bd and ad~bc separately; over the two non-tree generators {b,d}
    # the integer SNF of [[-1,-1],[-1,1]] is diag(1,2) => pi1^ab = Z/2. This also
    # exercises the char-divisibility branch of hom_to_additive_dim end-to-end (was
    # only synthetic before): dim Hom(pi1, k+) = 1 over char 2, 0 over char 3.
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (1, 2), "c": (2, 3), "d": (2, 3)})
    A = Q.algebra(relations=["a*c - b*d", "a*d - b*c"], field=GF(7))
    g = fundamental_group(A)
    assert g.free_rank == 0 and g.invariant_factors == (2,)   # pi1^ab = Z/2
    assert g.abelianization_repr() == "Z/2"
    assert g.hom_to_additive_dim(2) == 1                       # char 2 | 2
    assert g.hom_to_additive_dim(3) == 0                       # char 3 does not divide 2
    assert A.quiver.is_acyclic()                               # triangular


def test_intrinsic_refused_loudly():
    A = truncated_polynomial(3, field=GF(7))
    with pytest.raises(QuiverlabError, match="intrinsic"):
        intrinsic_fundamental_group(A)


def test_minimal_relation_counts_self_cert():
    # ADDENDUM (P62 contract): the total per-(x,y) minimal-relation count equals
    # sum over blocks; hereditary => empty; the commutative square WITH relation has
    # exactly one minimal relation (1 -> 4).
    from quiverlab.invariants.coverings import minimal_relation_counts, _relation_data
    A = _square(["a*b - c*d"])
    counts = minimal_relation_counts(A)
    assert counts == {(1, 4): 1}
    assert sum(counts.values()) == sum(_relation_data(A).counts.values())
    assert minimal_relation_counts(_square([])) == {}          # hereditary => empty


def test_minimal_relation_counts_presentationless_refuses():
    from quiverlab import Algebra
    from quiverlab.invariants.coverings import minimal_relation_counts
    # k[x]/(x^2) as bare structure constants (no quiver presentation).
    T = [[["1", "0"], ["0", "1"]], [["0", "1"], ["0", "0"]]]
    B = Algebra.from_structure_constants(T, ["1", "0"], field=GF(7))
    assert B.quiver is None
    with pytest.raises(QuiverlabError, match="presentation"):
        minimal_relation_counts(B)
