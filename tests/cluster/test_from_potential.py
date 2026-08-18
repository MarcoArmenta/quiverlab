# SPDX-License-Identifier: MIT
"""`from_potential`: the Jacobi-finite (Q,W) ALGEBRA-level certificate, and the narrow
certified reduction to a hereditary model (Plan 79 / R31, Task 5).

Amiot: when (Q,W) is Jacobi-finite, C_{(Q,W)} carries a cluster-tilting object whose
endomorphism algebra IS Jac(Q,W) -- cited, not computed. Keller-Reiten: cluster-tilted
algebras are Gorenstein of dimension at most one. Both are certified here at the ALGEBRA
level; the CATEGORY invariants are served only through a certified hereditary model.
"""
import pytest

from quiverlab import Quiver
from quiverlab.cluster import ClusterCategory
from quiverlab.errors import NotFiniteDimensionalError, QuiverlabError
from quiverlab.families.jacobian import Potential
from quiverlab.fields import QQ

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _three_cycle():
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)})
    return Q, Potential(Q, [(1, ("a", "b", "c"))])


@lit
def test_the_flagship_3_cycle_is_kZ3_mod_J2():
    Q, W = _three_cycle()
    C = ClusterCategory.from_potential(Q, W, field=QQ)
    cert = C.jacobian_certificate()
    assert cert["dim"] == 6                        # Jac(3-cycle, abc) = kZ_3/J^2
    assert cert["self_injective"] is True
    assert cert["finite_dimensional"] is True


@lit
def test_keller_reiten_gorenstein_at_most_one():
    Q, W = _three_cycle()
    cert = ClusterCategory.from_potential(Q, W, field=QQ).jacobian_certificate()
    assert cert["gorenstein"] is True
    assert cert["gorenstein_at_most_one"] is True
    assert cert["gorenstein_dimension"] == 0       # self-injective => Gorenstein dim 0


@selfcert
def test_the_flagship_is_certified_reducible_and_agrees_with_the_A3_route():
    # The whole point of the narrow recognition: (3-cycle, abc) IS the cluster-tilted
    # algebra of A3, so its cluster category is C_{A3} and the combinatorics must MATCH
    # the algebra route exactly -- computed twice, by different entry points.
    Q, W = _three_cycle()
    C = ClusterCategory.from_potential(Q, W, field=QQ)
    assert C.jacobian_certificate()["reducible_to_hereditary"] is True
    assert C.jacobian_certificate()["model"] == "A3"
    direct = ClusterCategory("A3")
    assert C.num_indecomposables() == direct.num_indecomposables() == 9
    assert (C.num_cluster_tilting()["count"]
            == direct.num_cluster_tilting()["count"] == 14)


@selfcert
def test_G6_a_jacobi_infinite_potential_raises_the_shipped_refusal():
    # One vertex, two loops, W = 0: Jac = the free algebra on two generators, infinite.
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    W = Potential(Q, [])
    with pytest.raises(NotFiniteDimensionalError):
        ClusterCategory.from_potential(Q, W, field=QQ)


@selfcert
def test_G7_a_non_reducible_QW_refuses_CATEGORY_invariants_but_keeps_the_algebra_one():
    # A 4-cycle with its own cyclic potential: Jacobi-finite, but NOT certified reducible
    # to a hereditary model. The algebra-level Amiot/Keller-Reiten certificate still
    # stands; every category invariant refuses, naming D^b(Gamma) as what is missing.
    Q = Quiver([1, 2, 3, 4],
               {"a": (1, 2), "b": (2, 3), "c": (3, 4), "d": (4, 1)})
    W = Potential(Q, [(1, ("a", "b", "c", "d"))])
    C = ClusterCategory.from_potential(Q, W, field=QQ)
    cert = C.jacobian_certificate()
    assert cert["reducible_to_hereditary"] is False and cert["model"] is None
    assert cert["dim"] > 0 and cert["finite_dimensional"] is True
    for call in (C.num_indecomposables, C.num_cluster_tilting):
        with pytest.raises(QuiverlabError, match=r"D\^b\(Gamma\)"):
            call()


@selfcert
def test_the_certificate_says_the_amiot_identification_is_CITED():
    Q, W = _three_cycle()
    note = ClusterCategory.from_potential(Q, W, field=QQ).jacobian_certificate()["note"]
    assert "Amiot" in note and "Keller-Reiten" in note
    assert "CITED, not computed" in note


@selfcert
def test_jacobian_certificate_is_refused_on_an_acyclic_input():
    # It is only defined for a from_potential construction; an acyclic C_Q has no (Q,W).
    with pytest.raises(QuiverlabError, match="from_potential"):
        ClusterCategory("A3").jacobian_certificate()


@selfcert
def test_the_reduction_is_certified_not_guessed():
    # A 3-cycle whose potential is EMPTY has Jac = kQ/(0) -- infinite-dimensional, so it
    # never reaches the recognizer at all. The recognition is gated on the built algebra's
    # dimension (6), not merely on the quiver's shape.
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)})
    with pytest.raises(NotFiniteDimensionalError):
        ClusterCategory.from_potential(Q, Potential(Q, []), field=QQ)
