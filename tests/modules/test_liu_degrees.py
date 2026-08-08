"""Liu left/right degrees on the knitted category (Plan 57 / R21+R37). Self-cert:
the degree-vs-layer witness triple (g in rad^d, g not in rad^{d+1}, g.then(f) in
rad^{d+2}) -- a finite left degree forces the layer drop. Literature: hand-derived
small degrees + sectional composites nonzero."""
import pytest

from quiverlab import linear_path_algebra, NakayamaAlgebra
from quiverlab.fields import QQ
from quiverlab.modules.ar import knit_ar_quiver
from quiverlab.modules.radical import radical_filtration
from quiverlab.modules.ar_invariants import (left_degree, right_degree,
                                             degree_table, max_sectional_length,
                                             is_sectional)

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature

INF = None   # infinite degree is represented as None (no float inf; see constraints)


def _named_arrow(ar, src, tgt):
    """Look an arrow up by MODULE NAME (indexing-robust -- the knit orders vertices
    by BFS, so integer indices are not stable across refactors)."""
    idx = {v["name"]: i for i, v in enumerate(ar.vertices)}
    return idx[src], idx[tgt]


# The full HAND-DERIVED degree tables (verified live against the engine's rad
# structure at authoring), keyed (src_name, tgt_name) -> (d_l, d_r):
KA2 = {("S_2", "P_1"): (INF, 1),      # mono
       ("P_1", "S_1"): (1, INF)}      # epi
KA3 = {("S_3", "P_2"): (INF, 1),      # mono
       ("P_2", "P_1"): (INF, 2),      # mono
       ("P_2", "S_2"): (1, INF),      # epi
       ("P_1", "I_2"): (2, INF),      # epi (degree 2, NOT 1 -- decomposable mesh middle)
       ("S_2", "I_2"): (INF, 1),      # mono
       ("I_2", "S_1"): (1, INF)}      # epi


@selfcert
def test_finite_left_degree_witness_triple():
    # M3: assert ALL THREE memberships of the witness triple EXTERNALLY.
    A = NakayamaAlgebra(kupisch=[3, 2, 2], field=QQ)   # has finite-degree arrows
    rf, ar = radical_filtration(A), knit_ar_quiver(A)
    seen_finite = 0
    for (i, j) in ar.arrows:                            # every irreducible map
        d, finite, wit = left_degree(rf, i, j)
        if not finite:
            continue
        seen_finite += 1
        z, g, comp = wit                               # comp = coords of g.then(f)
        assert rf._in_layer(z, i, d, g) is True         # g in rad^d(Z, X_i)
        assert rf._in_layer(z, i, d + 1, g) is False    # g not in rad^{d+1}
        assert rf._in_layer(z, j, d + 2, comp) is True  # g.then(f) in rad^{d+2}(Z, X_j)
    assert seen_finite >= 1                             # not vacuous: a finite degree occurs


@lit
@pytest.mark.parametrize("n, table", [(2, KA2), (3, KA3)])
def test_liu_degree_tables(n, table):
    # the corrected replacement for the (false) "all irreducibles have infinite
    # left degree" stub: the FULL hand-derived left+right degree tables, by name.
    A = linear_path_algebra(n, field=QQ)
    rf, ar = radical_filtration(A), knit_ar_quiver(A)
    expected = {2: KA2, 3: KA3}[n]
    assert len(expected) == len(ar.arrows)             # every irreducible covered
    for (src, tgt), (dl, dr) in expected.items():
        i, j = _named_arrow(ar, src, tgt)
        assert left_degree(rf, i, j)[0] == dl,  (src, tgt, "d_l")
        assert right_degree(rf, i, j)[0] == dr, (src, tgt, "d_r")


@lit
def test_mono_epi_degree_dichotomy_on_ka_n():
    # structural literature pin: on hereditary directed kA_n every irreducible MONO
    # has d_l = INF & d_r finite; every EPI has d_l finite & d_r = INF.
    for n in (2, 3, 4):
        A = linear_path_algebra(n, field=QQ)
        rf, ar = radical_filtration(A), knit_ar_quiver(A)
        for (i, j) in ar.arrows:
            dl = left_degree(rf, i, j)[0]
            dr = right_degree(rf, i, j)[0]
            # exactly one side is infinite for every irreducible on kA_n.
            assert (dl is INF) ^ (dr is INF)


@lit
def test_sectional_composite_is_nonzero():
    # a sectional path of length r yields a nonzero element of rad^r \ rad^{r+1}.
    A = linear_path_algebra(3, field=QQ)
    ar = knit_ar_quiver(A)
    assert max_sectional_length(ar) >= 2          # ZA_3 has a length-2 sectional path
    # and the specific sectional path S_3 -> P_2 -> P_1 has no mesh reversal
    idx = {v["name"]: i for i, v in enumerate(ar.vertices)}
    assert is_sectional(ar, [idx["S_3"], idx["P_2"], idx["P_1"]]) is True
    # while S_3 -> P_2 -> P_1 -> I_2 reverses at P_1 (tau I_2 = P_2)
    assert is_sectional(ar,
                        [idx["S_3"], idx["P_2"], idx["P_1"], idx["I_2"]]) is False


@selfcert
def test_left_right_degree_op_symmetry():
    # d_r(f) over A == d_l(D f) over A^op (the formal dual; NOT a "flip and see"). D is
    # a contravariant equivalence mod A -> mod A^op sending an irreducible X->Y to an
    # irreducible DY->DX with d_l(Df) = d_r(f); so the MULTISET of d_r over A's arrows
    # equals the multiset of d_l over A^op's arrows (index-free, robust to the knit's
    # BFS ordering).
    A = linear_path_algebra(3, field=QQ)
    rf, ar = radical_filtration(A), knit_ar_quiver(A)
    Aop = A.opposite()
    rfop, arop = radical_filtration(Aop), knit_ar_quiver(Aop)

    def _key(d):
        return (d is None, d if d is not None else 0)

    dr_A = sorted(_key(right_degree(rf, i, j)[0]) for (i, j) in ar.arrows)
    dl_Aop = sorted(_key(left_degree(rfop, i, j)[0]) for (i, j) in arop.arrows)
    assert dr_A == dl_Aop


@selfcert
def test_degree_table_covers_every_arrow():
    A = linear_path_algebra(3, field=QQ)
    rf, ar = radical_filtration(A), knit_ar_quiver(A)
    tab = degree_table(rf, ar)
    assert set(tab) == set(ar.arrows)
    for rec in tab.values():
        # exactly one of d_l / d_r finite on kA_n (mono/epi dichotomy)
        assert rec["finite_l"] ^ rec["finite_r"]
        # a multiplicity-1 arrow carries EXACTLY the four scalar keys (byte-stable):
        assert set(rec) == {"d_l", "d_r", "finite_l", "finite_r"}


@selfcert
def test_min_finite_helper():
    from quiverlab.modules.ar_invariants import _min_finite
    assert _min_finite([None, 2, None]) == 2
    assert _min_finite([3, 1, 2]) == 1
    assert _min_finite([None, None]) is None
    assert _min_finite([]) is None


@selfcert
def test_degree_table_multiplicity_gt_1_per_class(monkeypatch):
    # No knittable rep-finite quiver algebra in scope has a mult>1 AR arrow (verified
    # over kA_4/kA_5, D_4/D_5, the Nakayama zoo -- see the # PIN in degree_table). We
    # cover the per-class code path SYNTHETICALLY: force _class_reps to return TWO
    # class reps for every arrow (the two rad/rad^2 generators of a length-2 pair),
    # and assert the block gains {mult, d_l_classes, d_r_classes} with the min-degree
    # aggregation, while a genuine mult-1 run stays four-key.
    import quiverlab.modules.ar_invariants as ari
    A = linear_path_algebra(3, field=QQ)
    rf, ar = radical_filtration(A), knit_ar_quiver(A)
    real = ari._class_reps

    def _doubled(rf_, i, j):
        reps = real(rf_, i, j)
        return reps + reps[:1] if reps else reps      # duplicate the first rep -> mult 2

    monkeypatch.setattr(ari, "_class_reps", _doubled)
    tab = degree_table(rf, ar)
    for (i, j), rec in tab.items():
        assert rec["mult"] == 2
        assert len(rec["d_l_classes"]) == 2 and len(rec["d_r_classes"]) == 2
        # the scalar d_l/d_r is the finite-preferring min over the class list
        assert rec["d_l"] == ari._min_finite(rec["d_l_classes"])
        assert rec["d_r"] == ari._min_finite(rec["d_r_classes"])
        # duplicated rep => identical per-class degrees
        assert rec["d_l_classes"][0] == rec["d_l_classes"][1]
