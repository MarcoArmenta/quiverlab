# SPDX-License-Identifier: MIT
"""The cluster category's finite fundamental domain (Plan 79 / R31, Task 1).

`ind(C_Q) = ind(mod kQ) |_| {P_v[1]}` (BMRRT), so `#indec(C_Q) = #ind(mod kQ) + n` --
the almost-positive roots of the underlying diagram. Deep bucket (tests/cluster is in
`_DEEP_DIRS`): the AR knit is seconds-scale on D4/D5.
"""
import pytest

from quiverlab import Quiver
from quiverlab.cluster import ClusterCategory
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@lit
@pytest.mark.parametrize("typ,expected", [
    ("A2", 5), ("A3", 9), ("A4", 14), ("A5", 20),
    ("D4", 16),
    # D5's #indec CLOSES (20 positive roots + 5) even though its cluster-tilting count
    # is a separate certificate -- the two must not be conflated.
    ("D5", 25),
])
def test_num_indecomposables_is_the_almost_positive_root_count(typ, expected):
    C = ClusterCategory(typ)
    assert C.num_indecomposables() == expected
    assert C.almost_positive_roots_count() == expected


@lit
@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_type_A_matches_the_closed_form_n_times_n_plus_3_over_2(n):
    # #almost-positive-roots of A_n = n(n+3)/2.
    assert ClusterCategory(f"A{n}").num_indecomposables() == n * (n + 3) // 2


@lit
@pytest.mark.parametrize("n", [4, 5])
def test_type_D_matches_the_closed_form_n_squared(n):
    # D_n has n(n-1) positive roots, plus n shifted projectives = n^2.
    assert ClusterCategory(f"D{n}").num_indecomposables() == n * n


@selfcert
@pytest.mark.parametrize("typ", ["A3", "D4"])
def test_the_domain_splits_into_modules_plus_exactly_n_shifts(typ):
    C = ClusterCategory(typ)
    ind = C.indecomposables()
    kinds = [k for k, _ in ind]
    assert len(ind) == C.num_indecomposables()
    assert kinds.count("shift") == C.n
    assert kinds.count("module") == C.num_indecomposables() - C.n
    # the shifted entries carry the PROJECTIVES themselves, one per vertex
    shifted = [M for k, M in ind if k == "shift"]
    assert len(shifted) == len(set(M.name for M in shifted)) == C.n


@selfcert
def test_the_identity_indec_equals_ind_mod_A_plus_n():
    from quiverlab.modules.ar import knit_ar_quiver
    for typ in ("A3", "A4", "D4"):
        C = ClusterCategory(typ)
        knit = knit_ar_quiver(C.algebra)
        assert knit.status == "complete"
        assert C.num_indecomposables() == len(knit.vertices) + C.n


@selfcert
def test_G2_a_quiver_with_an_oriented_cycle_is_refused():
    # kQ is infinite-dimensional and C_Q is undefined; the refusal must come from the
    # cluster layer's own acyclicity gate, named as such.
    Q = Quiver([1, 2], {"a": (1, 2), "b": (2, 1)})
    with pytest.raises(QuiverlabError, match="oriented cycle"):
        ClusterCategory(Q, field=QQ)


@selfcert
def test_G2_a_loop_counts_as_an_oriented_cycle():
    Q = Quiver([1], {"x": (1, 1)})
    with pytest.raises(QuiverlabError, match="oriented cycle"):
        ClusterCategory(Q, field=QQ)


@selfcert
def test_G1_a_non_hereditary_algebra_is_refused():
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=QQ)
    with pytest.raises(QuiverlabError, match="hereditary"):
        ClusterCategory(A)


@selfcert
def test_G3b_a_representation_infinite_quiver_refuses_as_INFINITE():
    # The affine Kronecker ~A_1 (two parallel arrows) is tame-infinite: by Gabriel the
    # non-Dynkin diagrams are exactly the representation-infinite ones, so the claim is
    # made off the TYPE CERTIFICATE -- never off a truncated enumeration.
    Q = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)})
    C = ClusterCategory(Q, field=QQ)
    assert C.is_representation_finite is False
    with pytest.raises(QuiverlabError, match="infinitely many indecomposables"):
        C.num_indecomposables()


@selfcert
def test_G3b_message_names_gabriel_and_the_diagram_not_a_budget():
    Q = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)})
    C = ClusterCategory(Q, field=QQ)
    try:
        C.num_indecomposables()
    except QuiverlabError as exc:
        msg = str(exc)
    assert "Gabriel" in msg and "~A1" in msg
    # an infiniteness claim must NEVER be dressed up as a budget stop
    assert "budget" not in msg.lower()


@selfcert
def test_a_dynkin_string_a_quiver_and_an_algebra_agree():
    from quiverlab.families.path_algebra import PathAlgebra
    by_string = ClusterCategory("A3")
    by_quiver = ClusterCategory(
        Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}), field=QQ)
    by_algebra = ClusterCategory(PathAlgebra("A3", field=QQ))
    counts = {c.num_indecomposables() for c in (by_string, by_quiver, by_algebra)}
    assert counts == {9}


@selfcert
def test_presentation_less_algebra_is_refused():
    from quiverlab import Algebra
    one = QQ.one()
    with pytest.raises(QuiverlabError, match="quiver presentation"):
        ClusterCategory(Algebra(QQ, [[[one]]], [one]))
