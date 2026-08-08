"""Plan 72 part (b): certified arrow removal / addition Hochschild reductions (R6).

CLMS 1812.07655 / Proc. AMS 148 (2020) 2421-2432. Oracles (Plan-32 markers):
- oracle_literature: the frozen HH_(A)=HH_(B) homology iso on P1-P4, the cohomology
  Ext-correction / center-disconnection deltas (P1/P2/P3).
- oracle_crossengine: HH_{>=2}(A) == HH_{>=2}(B) == direct; the T-flagship crosscheck.
- oracle_selfcert: inert-arrow witnesses (incl. the P4 binomial), dim drop, the
  addition round-trip, HH_0 provable invariance, loud refusals.
"""
import pytest

from quiverlab import GF, Quiver
from quiverlab.fields import QQ
from quiverlab.families import TrivialExtension
from quiverlab.hochschild.arrow_removal import (
    inert_arrows, remove_arrows, add_arrows, arrow_removal, _relation_to_grammar,
)


def build(verts, arrows, rels, field=QQ):
    return Quiver(verts, arrows).algebra(relations=rels, field=field)


def P1(field=QQ):   # loop x at 1 (x^2=0) + inert c:1->2
    return build([1, 2], {"x": (1, 1), "c": (1, 2)}, ["x^2"], field)


def P2(field=QQ):   # 1->2->3 (ab=0) + inert c:1->3
    return build([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (1, 3)}, ["a*b"], field)


def P3(field=QQ):   # loops x@1, y@2 (x^2=y^2=0) + inert c:1->2
    return build([1, 2], {"x": (1, 1), "y": (2, 2), "c": (1, 2)}, ["x^2", "y^2"], field)


def P4(field=QQ):   # commutative square a*c - b*d (binomial) + inert e:1->4
    return build([1, 2, 3, 4],
                 {"a": (1, 2), "b": (1, 3), "c": (2, 4), "d": (3, 4), "e": (1, 4)},
                 ["a*c - b*d"], field)


# =====================================================================
# Task B1 -- inert_arrows / remove_arrows / add_arrows
# =====================================================================
@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("A,inert,dimA,dimB", [
    (P1(), ["c"], 5, 3),
    (P2(), ["c"], 6, 5),
    (P3(), ["c"], 8, 4),
    (P4(), ["e"], 10, 9),        # the non-monomial (binomial) inert-detection PIN
])
def test_inert_arrows_and_dims(A, inert, dimA, dimB):
    assert inert_arrows(A) == inert
    assert A.dim == dimA
    B = remove_arrows(A, inert)
    assert B.dim == dimB


@pytest.mark.oracle_selfcert
def test_non_inert_arrows_are_detected():
    # x,a,b are in relations -> NOT inert; only the bridge arrows are.
    assert "x" not in inert_arrows(P1())
    assert set(inert_arrows(P2())) == {"c"} and "a" not in inert_arrows(P2())


@pytest.mark.oracle_selfcert
def test_remove_non_inert_refuses_with_witness():
    with pytest.raises(Exception) as exc:
        remove_arrows(P1(), ["x"])               # x appears in x^2
    assert "not inert" in str(exc.value) or "relation" in str(exc.value)


@pytest.mark.oracle_selfcert
def test_binomial_relation_scanned_correctly():
    # P4: a,b,c,d all appear across the two terms of a*c - b*d; only e is inert.
    A = P4()
    assert inert_arrows(A) == ["e"]
    assert _relation_to_grammar(A.relations[0]).replace(" ", "") in ("a*c-b*d", "a*c-1*b*d")


@pytest.mark.oracle_selfcert
def test_addition_round_trip():
    # add_arrows(remove_arrows(A, D), D_as_new) ~= A (dim + HH)
    for Af, arrow, ends in [(P1(), "c", (1, 2)), (P2(), "c", (1, 3))]:
        B = remove_arrows(Af, [arrow])
        A2 = add_arrows(B, {arrow: ends})
        assert A2.dim == Af.dim
        assert A2.hochschild_homology(4).dims == Af.hochschild_homology(4).dims


@pytest.mark.oracle_selfcert
def test_relative_cycle_addition_refuses():
    kx2 = build([1], {"x": (1, 1)}, ["x^2"])
    with pytest.raises(Exception) as exc:
        add_arrows(kx2, {"y": (1, 1)})           # a second loop -> relative cycle
    assert "relative cycle" in str(exc.value) or "infinite" in str(exc.value)


# =====================================================================
# Task B2 -- the certified HH reduction (Thm 3.2 homology / Thm 4.2 coh)
# =====================================================================
@pytest.mark.oracle_literature
@pytest.mark.oracle_crossengine
def test_P1_reduction():
    r = P1().arrow_removal(top=5)
    assert r.removed == ["c"]
    assert r.hom_A == r.hom_B == [3, 1, 1, 1, 1, 1]   # homology iso in ALL degrees
    assert r.hom_agrees is True and r.hom_low_agrees == [True, True]
    assert r.coh_A == [1, 1, 1, 1, 1, 1] and r.coh_B == [3, 1, 1, 1, 1, 1]
    assert r.coh_low_delta[0] == -2                  # center/disconnection at n=0
    assert r.coh_correction[2:] == [0, 0, 0, 0]      # n>=2 Ext correction is 0


@pytest.mark.oracle_literature
@pytest.mark.oracle_crossengine
def test_P2_reduction_asymmetry():
    r = P2().arrow_removal(top=5)
    assert r.hom_A == r.hom_B == [3, 0, 0, 0, 0, 0]
    assert r.coh_A == [1, 1, 1, 0, 0, 0] and r.coh_B == [1, 0, 0, 0, 0, 0]
    assert r.coh_low_delta[1] == 1                   # n=1 center effect
    assert r.coh_correction[2] == 1                  # n=2 Thm-4.2 Ext term, NONZERO


@pytest.mark.oracle_literature
@pytest.mark.oracle_crossengine
def test_P3_low_degree_boundary():
    r = P3().arrow_removal(top=5)
    assert r.hom_A == r.hom_B == [4, 2, 2, 2, 2, 2]  # HH_0 = 4 fixed despite disconnection
    assert r.coh_A == [1, 3, 2, 2, 2, 2] and r.coh_B == [4, 2, 2, 2, 2, 2]
    assert r.coh_low_delta == [-3, 1]                # disconnection at 0, bridge at 1
    assert r.coh_correction[2:] == [0, 0, 0, 0]      # n>=2 correction is 0


@pytest.mark.oracle_literature
def test_P4_homology():
    r = P4().arrow_removal(top=5)
    assert r.removed == ["e"]
    assert r.hom_A == r.hom_B == [4, 0, 0, 0, 0, 0]


@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("A", [P1(), P2(), P3(), P4()])
def test_hh0_provably_invariant(A):
    # inert removal never changes HH_0 (no cyclic path through an inert arrow).
    r = A.arrow_removal(top=3)
    assert r.hom_A[0] == r.hom_B[0]
    assert r.hom_low_agrees[0] is True


@pytest.mark.oracle_crossengine
def test_trivial_extension_crosscheck():
    # T(kA2) + an inert arrow to a fresh sink vertex reduces BACK to T(kA2):
    # HH_{>=2} preserved (ties R6 to the R5 flagship).
    B = Quiver([1, 2], {"a1": (1, 2)}).algebra(field=QQ)
    T = TrivialExtension(B)
    verts = list(T.quiver.vertices) + [3]
    arrows = dict(T.quiver.arrows)
    arrows["ee"] = (1, 3)                            # inert; vertex 3 is a sink -> finite
    rels = [_relation_to_grammar(r) for r in T.relations]
    A = Quiver(verts, arrows).algebra(relations=rels, field=QQ)
    assert inert_arrows(A) == ["ee"]
    r = A.arrow_removal(top=4)
    assert r.hom_agrees is True                      # HH_{>=2}(A) == HH_{>=2}(A\ee)
    assert r.hom_B[2:] == T.hochschild_homology(4).dims[2:]


@pytest.mark.oracle_crossengine
def test_gf32003_parity_P2():
    r = P2(GF(32003)).arrow_removal(top=4)
    assert r.hom_A == r.hom_B == [3, 0, 0, 0, 0]
    assert r.hom_agrees is True
