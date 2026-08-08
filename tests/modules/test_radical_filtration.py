"""The radical filtration of mod A (Plan 57 / R37). Self-cert: rad^{n+1} subset
rad^n, rad^N = 0 while rad^{N-1} != 0, composition closure. Literature: kA_n layer
tables + nilpotency index N(kA_n)=n; kA_3/J^2 index 3 (Thm 1.3). Cross-engine: the
mesh closed form == route (ii); Thm-1.3 index == route (ii). Honest scope:
self-injective refused (status 'unsupported'); rep-infinite -> window, no verdict."""
import pytest

from quiverlab import Quiver, linear_path_algebra, NakayamaAlgebra
from quiverlab.fields import QQ, GF
from quiverlab.modules.ar import knit_ar_quiver
from quiverlab.modules.hom import hom_dim
from quiverlab.modules.radical import (radical_filtration, _thm13_index,
                                       _mesh_layer_dim)

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _name_idx(rf):
    return {nm: i for i, nm in enumerate(rf.names)}


def _total(rf, n):
    r = len(rf.indecs)
    return sum(rf.layer_dim(i, j, n) for i in range(r) for j in range(r))


@lit
@pytest.mark.parametrize("n", [2, 3, 4])
def test_nilpotency_index_of_linear_a_n_is_n(n):
    rf = radical_filtration(linear_path_algebra(n, field=QQ))
    assert rf.is_complete and rf.status == "complete"
    assert rf.nilpotency_index == n            # Thm 1.5(a); hand-derived n=2,3
    assert rf.rad_infinity_zero is True         # Auslander witness = N


@lit
def test_ka3_layer_table_matches_the_mesh():
    A = linear_path_algebra(3, field=QQ)
    rf = radical_filtration(A)
    ni = _name_idx(rf)
    # ZA_3 mesh (route (i), hand-derived): total dim rad^1 = 9 (all radical maps),
    # rad^2 = 3 (the length-2 composites S_3->P_1, P_2->I_2, P_1->S_1), rad^3 = 0.
    assert _total(rf, 1) == 9
    assert _total(rf, 2) == 3
    assert _total(rf, 3) == 0
    # the DESCENDING pair sequences of the three depth-2 maps: [1, 1] (in rad^1 and rad^2)
    for src, tgt in (("P_1", "S_1"), ("S_3", "P_1"), ("P_2", "I_2")):
        assert rf.layer[(ni[src], ni[tgt])] == [1, 1], (src, tgt)
    # an irreducible arrow that is NOT a longer composite: depth 1 only
    assert rf.layer[(ni["P_1"], ni["I_2"])] == [1]
    # no map 2->1 direction: S_2 -> S_1 is zero at every layer
    assert rf.layer.get((ni["S_2"], ni["S_1"]), []) == []


@lit
def test_ka3_radsq_truncated_nakayama_index_is_3():
    # kA_3/J^2 (1->2->3, ab=0): NOT self-injective, 5 indecomposables, index 3.
    A = NakayamaAlgebra(n=3, l=2, cyclic=False, field=QQ)
    rf = radical_filtration(A)
    assert rf.is_complete and len(rf.indecs) == 5
    assert rf.nilpotency_index == 3            # Thm 1.3: max{2,3,2}


@selfcert
def test_layers_are_a_descending_filtration_and_terminate():
    rf = radical_filtration(linear_path_algebra(3, field=QQ))
    r = len(rf.indecs)
    N = rf.nilpotency_index
    for i in range(r):
        for j in range(r):
            seq = rf.layer[(i, j)]
            for n in range(1, len(seq)):
                assert seq[n] <= seq[n - 1]            # rad^{n+1} subset rad^n
    # rad^{N-1} != 0 somewhere, rad^N == 0 everywhere
    assert any(rf.layer_dim(i, j, N - 1) > 0 for i in range(r) for j in range(r))
    assert all(rf.layer_dim(i, j, N) == 0 for i in range(r) for j in range(r))


@selfcert
def test_composition_closure_rad_next_inside_rad():
    # rad^{n+1}(X,Y) subset rad^n(X,Y): every basis vector of rad^{n+1} lies in rad^n.
    rf = radical_filtration(linear_path_algebra(3, field=QQ))
    r, N = len(rf.indecs), rf.nilpotency_index
    for n in range(1, N - 1):
        for i in range(r):
            for j in range(r):
                for vec in rf._layer_basis(i, j, n + 1):
                    assert rf._in_layer(i, j, n, vec) is True


@selfcert    # NOT cross-engine: the Thm-1.3 walk and route (ii) both read the SAME
def test_thm13_index_equals_route_ii():
    # knit (one engine, two readings). This certifies the index computation is
    # internally consistent with the Chaio-Guazzelli formula on our own knit.
    for A in (linear_path_algebra(3, field=QQ),
              NakayamaAlgebra(n=3, l=2, cyclic=False, field=QQ)):
        rf = radical_filtration(A)
        assert _thm13_index(A) == rf.nilpotency_index


@xeng
def test_mesh_layer_dims_match_route_ii_on_ka3():
    # GENUINELY INDEPENDENT route (i): the ZA_3 mesh closed form (pure interval
    # combinatorics on dimension vectors -- no Hom, no matmul) reproduces route (ii)'s
    # dim rad^n(X,Y) on kA_3 (standard directed => the mesh count is exact).
    A = linear_path_algebra(3, field=QQ)
    rf, ar = radical_filtration(A), knit_ar_quiver(A)
    r, N = len(rf.indecs), rf.nilpotency_index
    for i in range(r):
        for j in range(r):
            for n in range(1, N + 1):
                assert _mesh_layer_dim(ar, i, j, n) == rf.layer_dim(i, j, n)


@xeng
def test_mesh_layer_dims_match_route_ii_on_ka4():
    A = linear_path_algebra(4, field=QQ)
    rf, ar = radical_filtration(A), knit_ar_quiver(A)
    r, N = len(rf.indecs), rf.nilpotency_index
    for i in range(r):
        # kA_4 has a non-interval-free knit? no: every kA_n indecomposable is an
        # interval, so the mesh form applies to every vertex.
        for j in range(r):
            for n in range(1, N + 1):
                assert _mesh_layer_dim(ar, i, j, n) == rf.layer_dim(i, j, n)


@selfcert
def test_self_injective_refused_as_unsupported():
    # cyclic Nakayama kZ_3/J^2 is self-injective -> knit unsupported -> no verdict.
    A = NakayamaAlgebra(n=3, l=2, cyclic=True, field=QQ)
    rf = radical_filtration(A)
    assert rf.is_complete is False and rf.status == "unsupported"
    assert rf.nilpotency_index is None and rf.rad_infinity_zero is None


@selfcert
def test_representation_infinite_is_window_only_no_verdict():
    # 3-Kronecker (wild): a small budget_dim trips the knit cap FAST (the
    # preprojectives explode in dim -- the P41 test_ar_knit budget pattern), giving
    # window layers with no index/rad^inf verdict.
    A = Quiver([1, 2], {"a": (1, 2), "b": (1, 2), "c": (1, 2)}).algebra(
        relations=[], field=QQ)
    rf = radical_filtration(A, budget_modules=40, budget_dim=20)
    assert rf.is_complete is False and rf.status in ("budget", "error")
    assert rf.nilpotency_index is None and rf.rad_infinity_zero is None
    assert "window" in rf.note.lower()          # loud: NOT rad(mod A)


@selfcert
def test_layer1_equals_hom_off_diagonal_field_parity():
    for field in (QQ, GF(32003)):
        rf = radical_filtration(linear_path_algebra(3, field=field))
        # rad(X_i, X_j) = Hom(X_i, X_j) for i != j (all maps radical)
        for i, Xi in enumerate(rf.indecs):
            for j, Xj in enumerate(rf.indecs):
                if i != j:
                    assert rf.layer_dim(i, j, 1) == hom_dim(Xi, Xj)


@selfcert
def test_generalized_standard_true_in_rep_finite_scope():
    rf = radical_filtration(linear_path_algebra(3, field=QQ))
    assert rf.is_generalized_standard() is True
    # pair_layers matches layer_dim through the module accessor
    X, Y = rf.indecs[0], rf.indecs[1]
    assert rf.pair_layers(X, Y) == rf.layer[(0, 1)]
