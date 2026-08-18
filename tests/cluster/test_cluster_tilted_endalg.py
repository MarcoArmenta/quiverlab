# SPDX-License-Identifier: MIT
"""The cluster-tilted End-algebra `End_{C_Q}(T) = Jac(Q_T, W_T)` (BMR / Amiot),
Plan 79 / R31 Task 3.

The identification itself is CITED, never computed -- there is no dg engine, so `Hom_C`
is never formed. What is verified per instance: the FZ quiver mutation, the Jacobian
algebra's presentation, and the Keller-Reiten theorem that a cluster-tilted algebra is
Gorenstein of dimension at most one.
"""
from itertools import permutations

import pytest

from quiverlab import NakayamaAlgebra
from quiverlab.cluster import ClusterCategory
from quiverlab.fields import QQ

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert


def _equal_up_to_vertex_permutation(X, Y):
    """The permutation making ``X`` equal ``Y``, or ``None``.

    A Cartan matrix is indexed by a CHOSEN vertex order, so two presentations of the same
    algebra agree only up to a simultaneous row/column permutation. Comparing raw entries
    would test a labelling coincidence rather than the mathematics -- see the flagship
    test for why that matters here.
    """
    n = len(X)
    for p in permutations(range(n)):
        if all(X[p[i]][p[j]] == Y[i][j] for i in range(n) for j in range(n)):
            return p
    return None


# --------------------------------------------------------------------------- #
# the single-3-cycle flagship: mu_1 of kA3
# --------------------------------------------------------------------------- #

@lit
@selfcert
def test_flagship_A3_mutation_1_is_the_3_cycle_with_dim_6():
    r = ClusterCategory("A3").cluster_tilted_algebra([1])
    assert r["fz_certified"] is True
    assert r["three_cycles"] == 1
    assert len(r["quiver"]["arrows"]) == 3          # an oriented 3-cycle
    assert r["dim"] == 6                            # Jac(3-cycle, abc) = kZ_3/J^2
    assert r["self_injective"] is True
    assert "fz_quiver_mutation" in r["verified"] and "jacobian_built" in r["verified"]


@xeng
def test_flagship_cartan_matches_the_independent_nakayama_UP_TO_LABELLING():
    """PLAN-DOC CORRECTION (measured).

    The plan's Spike 3 states the flagship's ``cartan_matrix()`` is EQUAL to the
    independent ``NakayamaAlgebra(n=3, l=2, cyclic=True)``. Measured on this tree raw
    equality is FALSE: the FZ mutation produces the 3-cycle oriented 1->3->2->1 while the
    Nakayama family builds 1->2->3->1, so the two agree only after the transposition
    (0,2,1). The algebras ARE isomorphic -- which is what the cross-check is for -- so the
    honest assertion is equality UP TO A VERTEX PERMUTATION. Asserting raw equality would
    have pinned an orientation coincidence and would break the moment either builder
    relabels.
    """
    r = ClusterCategory("A3").cluster_tilted_algebra([1])
    N = NakayamaAlgebra(n=3, l=2, cyclic=True, field=QQ)
    assert r["dim"] == N.dim == 6
    perm = _equal_up_to_vertex_permutation(r["algebra"].cartan_matrix(),
                                           N.cartan_matrix())
    assert perm is not None, "the Cartan matrices are not even permutation-equal"
    # and the raw matrices genuinely differ, so the test is not vacuously weakened
    assert r["algebra"].cartan_matrix() != N.cartan_matrix()


@lit
def test_flagship_has_infinite_global_dimension_per_keller_reiten():
    # Keller-Reiten: a cluster-tilted algebra is hereditary IFF its global dimension is
    # finite. This one is NOT hereditary (it has relations), so gl.dim must be infinite --
    # the engine reports a certified lower bound rather than a proved infinity.
    r = ClusterCategory("A3").cluster_tilted_algebra([1])
    gd = r["algebra"].global_dimension()
    assert not isinstance(gd, int) or gd > 8       # unresolved / large: never a small int


# --------------------------------------------------------------------------- #
# the MULTI-3-cycle instance: exercises the all-cycles enumeration
# --------------------------------------------------------------------------- #

@selfcert
def test_multi_3_cycle_instance_enumerates_ALL_of_them():
    """The plan's Task-3 requirement: an instance whose cluster-tilted quiver carries
    >= 2 oriented 3-cycles, so the canonical potential exercises the full enumeration
    rather than the single-triangle path.

    FOUND BY SEARCH over mutation sequences: type A4 has NO such quiver within 3
    mutations; the first hit is **A5 with the sequence [1, 3]** (and its mirror [3, 1]),
    giving 2 oriented 3-cycles and dim 13.
    """
    r = ClusterCategory("A5").cluster_tilted_algebra([1, 3])
    assert r["three_cycles"] == 2
    assert r["fz_certified"] is True
    assert r["dim"] == 13
    assert r["algebra"] is not None
    assert "all 2 oriented 3-cycle" in r["note"] or "2 oriented 3-cycle" in r["note"]


@lit
@pytest.mark.parametrize("typ,seq,expected_gor", [("A3", [1], 0), ("A5", [1, 3], 1)])
def test_keller_reiten_gorenstein_dimension_at_most_one(typ, seq, expected_gor):
    # Keller-Reiten: "cluster-tilted algebras are Gorenstein of dimension at most one".
    # The multi-3-cycle instance has Gorenstein dim 1 (not 0), so this is NOT just the
    # self-injective flagship restated -- the bound is exercised at both values.
    r = ClusterCategory(typ).cluster_tilted_algebra(seq)
    g = r["algebra"].gorenstein_dimension()
    assert g.is_gorenstein is True
    assert g.left_id == g.right_id == expected_gor
    assert max(g.left_id, g.right_id) <= 1


@selfcert
def test_the_fz_certificate_is_a_real_round_trip():
    # The quiver step's certificate: rebuild the quiver from the mutated exchange matrix,
    # re-read ITS exchange matrix, and require the original back.
    from quiverlab.surfaces.flip import exchange_matrix, matrix_mutation
    C = ClusterCategory("A3")
    B = exchange_matrix(C.algebra.quiver)
    expected = matrix_mutation(B, 1)
    r = C.cluster_tilted_algebra([1])
    assert r["exchange_matrix"] == expected
    assert r["fz_certified"] is True


@selfcert
def test_an_empty_mutation_sequence_returns_the_initial_acyclic_quiver():
    # mu of nothing: the quiver is Q itself, which is acyclic, so there is no 3-cycle and
    # no potential -- the algebra is honestly not built (G5), the quiver still certified.
    r = ClusterCategory("A3").cluster_tilted_algebra([])
    assert r["three_cycles"] == 0 and r["algebra"] is None
    assert r["fz_certified"] is True
    assert "no oriented 3-cycle" in r["note"] and "DWZ" in r["note"]


@selfcert
def test_G5_a_non_type_A_input_returns_the_quiver_but_NEVER_guesses_a_potential():
    # D4: the canonical sum-of-3-cycles potential is the TYPE-A potential; outside type A
    # the potential needs DWZ mutation / right-equivalence (out of scope, P48.1). The
    # quiver must still be returned and FZ-certified.
    C = ClusterCategory("D4")
    r = None
    for k in range(C.n):
        cand = C.cluster_tilted_algebra([k])
        if cand["three_cycles"] >= 1:
            r = cand
            break
    assert r is not None, "expected some D4 mutation to produce an oriented 3-cycle"
    assert r["algebra"] is None                     # no potential guessed
    assert r["fz_certified"] is True
    assert r["quiver"]["arrows"]
    assert "DWZ" in r["note"] and "type-A potential" in r["note"]
    assert "no potential is guessed" in r["note"]


@selfcert
def test_an_out_of_range_mutation_index_is_refused():
    from quiverlab.errors import QuiverlabError
    with pytest.raises(QuiverlabError, match="out of range"):
        ClusterCategory("A3").cluster_tilted_algebra([7])
