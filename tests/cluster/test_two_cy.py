# SPDX-License-Identifier: MIT
"""The 2-Calabi-Yau certificate on the module window (Plan 79 / R31, Task 4).

`Ext^1_C(X,Y) = Ext^1_A(X,Y) (+) D Ext^1_A(Y,X)` (BMRRT), so `dim Ext^1_C` is symmetric
BY CONSTRUCTION -- that half is bookkeeping. The computable content is the
Auslander-Reiten formula `Ext^1_A(X,Y) = D Hom-bar(Y, tau X)` (ASS Thm IV.2.13), checked
here on every ORDERED pair with ordinary Hom as the theorem-backed proxy for the
injectively-stable Hom-bar.
"""
import pytest

from quiverlab.cluster import ClusterCategory
from quiverlab.modules.ext import ext_dims

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert


@selfcert
@pytest.mark.parametrize("typ,pairs", [("A2", 9), ("A3", 36), ("A4", 100)])
def test_the_AR_formula_holds_on_every_ordered_pair(typ, pairs):
    # THE evidence for the ordinary-Hom proxy: a 0-mismatch sweep, reported rather than
    # assumed. A mismatch would be listed in `mismatches`, never smoothed away.
    c = ClusterCategory(typ).two_cy_certificate()
    assert c["pairs_checked"] == pairs
    assert c["ar_formula_holds"] is True
    assert c["mismatches"] == []


@selfcert
def test_the_AR_sweep_also_holds_on_D4_and_a_zigzag_A4():
    # The critic-expanded sweep from the plan: type D and a NON-linear orientation, so the
    # proxy is not only exercised on linearly-oriented type A.
    from quiverlab import Quiver
    from quiverlab.fields import QQ
    for C in (ClusterCategory("D4"),
              ClusterCategory(Quiver([1, 2, 3, 4],
                                     {"a": (1, 2), "b": (3, 2), "c": (3, 4)}),
                              field=QQ)):
        c = C.two_cy_certificate()
        assert c["ar_formula_holds"] is True and c["mismatches"] == []


@lit
def test_A3_has_exactly_five_nonzero_Ext1_pairs():
    # The plan's Spike-4 pin, reproduced: on kA3 exactly 5 ordered pairs have
    # dim Ext^1_A = 1 and the rest vanish.
    c = ClusterCategory("A3").two_cy_certificate()
    nonzero = c["nonzero_ext1_A_pairs"]
    assert len(nonzero) == 5
    assert all(c["ext1_A"][f"{i},{j}"] == 1 for i, j in nonzero)


@selfcert
@pytest.mark.parametrize("typ", ["A2", "A3", "D4"])
def test_ext1_C_is_symmetric_and_is_the_sum_of_the_two_directions(typ):
    c = ClusterCategory(typ).two_cy_certificate()
    assert c["ext1_C_symmetric"] is True
    n = c["modules_checked"]
    for i in range(n):
        for j in range(n):
            assert (c["ext1_C"][f"{i},{j}"]
                    == c["ext1_A"][f"{i},{j}"] + c["ext1_A"][f"{j},{i}"])
            assert c["ext1_C"][f"{i},{j}"] == c["ext1_C"][f"{j},{i}"]


@xeng
def test_the_ext1_table_agrees_with_an_INDEPENDENT_modules_ext_recompute():
    # The drift gate: recompute every entry straight from modules.ext, outside the
    # certificate's own loop.
    C = ClusterCategory("A3")
    c = C.two_cy_certificate()
    mods = [M for kind, M in C.indecomposables() if kind == "module"]
    assert len(mods) == c["modules_checked"]
    for i, X in enumerate(mods):
        for j, Y in enumerate(mods):
            assert ext_dims(C.algebra, X, Y, 1)[1] == c["ext1_A"][f"{i},{j}"]


@selfcert
def test_the_verdict_is_SCOPED_and_never_a_bald_true():
    # The MODERATE ruling: is_2_calabi_yau() must not claim the categorical statement --
    # the shifted P_v[1] pairs are cited (BMRRT), not computed.
    v = ClusterCategory("A3").is_2_calabi_yau()
    assert v["verdict"] is True
    assert v["scope"] == "module_window"
    assert v["shifted_by_citation"] is True
    assert "cited, not computed" in v["note"]
    assert set(v) >= {"verdict", "scope", "shifted_by_citation", "pairs_checked", "note"}


@selfcert
def test_the_payload_says_the_symmetry_is_by_construction():
    # Honesty about which half is content and which is bookkeeping: a+b = b+a is not a
    # verified theorem and the note must not let it read as one.
    c = ClusterCategory("A2").two_cy_certificate()
    assert "BY CONSTRUCTION" in c["note"]
    assert "cited" in c["note"] and "NOT computed" in c["note"]


@selfcert
def test_indecomposables_hands_back_MODULES_not_AR_quiver_records():
    # REGRESSION. ARQuiver.vertices is a list of RECORD DICTS {"name","dimvec","module"};
    # the first draft put those dicts straight into the fundamental domain, so every
    # "module" entry was a dict. The Task-1 tests missed it because they only counted
    # entries and read .name on the SHIFTED ones (which are real projectives). It
    # surfaced only when the 2-CY sweep tried to call .tau() on one.
    for kind, M in ClusterCategory("A3").indecomposables():
        assert hasattr(M, "tau") and hasattr(M, "dimension_vector"), (kind, type(M))
        assert not isinstance(M, dict)
