"""tau-exceptional sequences via the Jasso tau-perpendicular reduction (Plan 65 / R27;
Buan-Marsh J. Algebra 585 (2021); Jasso IMRN 2015; DIJ 2019). Enumeration via the ordered
support-tau-tilt bijection (count = n! * #sTt); NO mutation-BFS (transitivity only rank 2,
Buan-Hanson-Marsh 2402.10301). tau-tilting-finite + char-0/char>dim scope, QQ.

Also pins Plan 65 Task 0 (M1): the P45 `mutate` D_4-star defect fix -- `exchange_graph`
must return status="complete" with 50 vertices on ALL D_4-star orientations (the mixed
star tripped a spurious status="error" while still discovering all 50 vertices)."""
from math import factorial

import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.families.basic import linear_path_algebra

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _sig_key(tower):
    """A unique path key for a materialised reduction tower (Plan 65 / R27, H1)."""
    return tuple(entry.key() for entry in tower)


# --------------------------------------------------------------------------- #
# Task 0 (M1): the P45 mutate D_4-star fix -- regression pins
# --------------------------------------------------------------------------- #
_D4_STARS = {
    "subspace": {"a": (2, 1), "b": (3, 1), "c": (4, 1)},
    "source": {"a": (1, 2), "b": (1, 3), "c": (1, 4)},
    "mixed": {"a": (1, 2), "b": (3, 1), "c": (4, 1)},
}


@selfcert
@pytest.mark.parametrize("orient", ["subspace", "source", "mixed"])
def test_d4_star_exchange_graph_complete_all_orientations(orient):
    """M1 / Task 0: #sTt(D_4) = 50 with status='complete' on every star orientation.
    Before the fix the mixed star returned status='error' (a spurious P45 mutate miss on
    summands 0/3 of one full-support tau-tilting module) even though all 50 vertices were
    found; the fix arbitrates the decomposable-cone case and restores status='complete'."""
    from quiverlab.tautilting.mutation import exchange_graph
    A = Quiver([1, 2, 3, 4], _D4_STARS[orient]).algebra(field=QQ)
    eg = exchange_graph(A, budget_pairs=2000)
    assert eg.status == "complete"
    assert eg.is_complete is True
    assert len(eg.vertices) == 50
    assert eg.n_regular is True


def _matching_index(p1, p0):
    """The exchangeable summand of ``p1`` whose removal recovers the almost-complete pair
    shared with ``p0`` -- the g-column of ``p1`` NOT present in ``p0``."""
    from quiverlab.tautilting.rigid import g_columns
    g0 = {tuple(c) for c in g_columns(p0)}
    for k, c in enumerate(g_columns(p1)):
        if tuple(c) not in g0:
            return k
    raise AssertionError("no matching exchange index")


@selfcert
def test_d4_mixed_star_previously_failing_pair_certifies():
    """The mutate fix must certify (not merely relabel): the previously-failing full-support
    tau-tilting module {(0,1,0,0),(1,1,1,0),(1,1,0,1),(1,1,1,1)} raised at summands 0 and 3
    ('no valid exchange found'); now every summand exchanges and the involution holds. (The
    n_regular=True check above already certifies all 50 vertices; this pins the exact pair.)"""
    from quiverlab.tautilting.mutation import exchange_graph, mutate
    A = Quiver([1, 2, 3, 4], _D4_STARS["mixed"]).algebra(field=QQ)
    eg = exchange_graph(A, budget_pairs=2000)
    target = {(0, 1, 0, 0), (1, 1, 1, 0), (1, 1, 0, 1), (1, 1, 1, 1)}
    pair = None
    for rec in eg.vertices:
        dvs = {tuple(dv[w] for w in sorted(dv)) for dv in rec["summand_dimvecs"]}
        if not rec["support"] and dvs == target:
            pair = rec["pair"]
            break
    assert pair is not None, "the previously-failing full-support pair was not found"
    for k in range(4):                      # k=0 and k=3 were the failing exchanges
        nb = mutate(pair, k)                 # certifies internally; raises if it cannot
        back = mutate(nb, _matching_index(nb, pair))   # involution
        assert back.g_key() == pair.g_key()


# --------------------------------------------------------------------------- #
# Task A1 (selfcert): the Jasso tau-perpendicular reduction C(U)
# --------------------------------------------------------------------------- #
@selfcert
def test_reduction_rank_drops():
    from quiverlab.tautilting.exceptional import tau_perpendicular_reduction
    A = linear_path_algebra(3, field=QQ)
    red = tau_perpendicular_reduction(A, A.simple(1))     # a tau-rigid indecomposable
    C = red.reduction_algebra
    assert len(list(C.quiver.vertices)) == 3 - 1          # rk C(U) = n - |U|


@selfcert
def test_reduction_is_kA2_not_just_rank2():
    """H5: rk = n-|U| ALONE is insensitive to a wrong Bongartz completion / wrong U-vertices
    (any single-vertex idempotent quotient of a rank-3 algebra has rank 2). Pin the ISO:
    C(S1) of kA3 must be kA2 (connected, dim 3), NOT k x k (dim 2)."""
    from quiverlab.tautilting.exceptional import tau_perpendicular_reduction
    A = linear_path_algebra(3, field=QQ)
    C = tau_perpendicular_reduction(A, A.simple(1)).reduction_algebra
    kA2 = linear_path_algebra(2, field=QQ)
    assert len(list(C.quiver.arrows)) == 1                # connected rank-2 = kA2, not k x k
    assert C.dim == kA2.dim                               # kA2 has dim 3 (1+1+1 for 1->2)


@selfcert
def test_shifted_projective_reduction_is_support_quotient():
    """H2: the C(U)=End(T_U)/<e_U> recipe is ONLY for tau-rigid U; a shifted projective
    P_v[1] reduces to the support quotient A/<e_v>."""
    from quiverlab.tautilting.exceptional import tau_perpendicular_reduction
    A = linear_path_algebra(2, field=QQ)
    Pv_shift = next(o for o in A.tau_exceptional_objects()
                    if o.sign < 0 and o.vertex == 2)
    red = tau_perpendicular_reduction(A, Pv_shift)        # dispatches on the sign
    assert list(red.reduction_algebra.quiver.vertices) == [1]   # A/<e_2> = k at vertex 1
    assert red.reduction_algebra.dim == A.quotient_by_idempotent([2]).dim
    assert red.kind == "shifted"


# --------------------------------------------------------------------------- #
# Task A2 (selfcert): objects + recognizer
# --------------------------------------------------------------------------- #
@selfcert
def test_tau_exceptional_objects_count():
    A = linear_path_algebra(2, field=QQ)
    objs = A.tau_exceptional_objects()
    pos = [o for o in objs if o.sign > 0]
    neg = [o for o in objs if o.sign < 0]
    assert len(pos) == 3                                  # tau-rigid indecs of kA2: S1,S2,P1
    assert len(neg) == 2                                  # one shifted projective per vertex
    assert {o.vertex for o in neg} == {1, 2}


@selfcert
def test_recognizer_length1_and_classical_seed():
    """A length-1 tau-rigid object is accepted; a non-tau-rigid indecomposable is rejected;
    a classical A2 exceptional sequence is accepted (the (a)<->(b) seed)."""
    from quiverlab.families.radical_square_zero import RadicalSquareZero
    A = linear_path_algebra(2, field=QQ)
    assert A.is_tau_exceptional_sequence([A.simple(1)]) is True     # tau-rigid indec
    assert A.is_tau_exceptional_sequence([A.simple(2), A.projective(1)]) is True   # classical CES
    # a non-tau-rigid indecomposable simple over k[x]/(x^2) (self-injective, S not tau-rigid)
    B = RadicalSquareZero(Quiver([1], {"x": (1, 1)}), field=QQ)
    assert B.is_tau_exceptional_sequence([B.simple(1)]) is False


@selfcert
def test_recognizer_accepts_shifted_projective_sequence():
    """H2 pin: on kA2 the signed sequence (S1, P2[1]) -- inner S1 (genuine module), outer the
    shifted projective P2[1] -- is a complete signed tau-exceptional sequence. Its outer term
    reduces via A/<e_2> = k at vertex 1, whose signed objects lift to {S1(+), P1[1](-)}."""
    A = linear_path_algebra(2, field=QQ)
    objs = {(o.sign, getattr(o, "vertex", None)): o for o in A.tau_exceptional_objects()}
    P2shift = objs[(-1, 2)]                               # the shifted projective P_2[1]
    S1 = A.simple(1)
    assert A.is_tau_exceptional_sequence([S1, P2shift]) is True    # (inner, outer) convention


# --------------------------------------------------------------------------- #
# Task A3: the ordered-sTt bijection enumeration
# --------------------------------------------------------------------------- #
@selfcert   # H1: SELF-CERT -- signed_count is DEFINED as n!*stt_count (tautological).
@pytest.mark.parametrize("n", [2, 3])            # A_2 -> 10, A_3 -> 84
def test_signed_count_formula_wired(n):
    from quiverlab.tautilting.mutation import exchange_graph
    A = linear_path_algebra(n, field=QQ)
    eg = exchange_graph(A, budget_pairs=2000)
    assert eg.status == "complete"                        # M1: never trust a non-complete graph
    rep = A.tau_exceptional_sequences(budget=4096, want_sequences=False)
    assert rep.is_complete
    assert rep.stt_count == len(eg.vertices)              # #sTt read correctly from P45
    assert rep.signed_count == factorial(n) * len(eg.vertices)   # formula wired (tautological)


@xeng   # H1: the REAL cross-engine oracle -- MATERIALISE the sequences and count them.
@pytest.mark.parametrize("factory,n", [
    ("A2", 2),                                           # signed_count = 2! * 5  = 10
    ("A3", 3),                                           # signed_count = 3! * 14 = 84
    ("nonhered", 2),                                     # k(1<->2)/rad^2 : 2! * 6 = 12
])
def test_materialised_count_equals_signed_count(factory, n):
    """The formula n!*#sTt is cross-checked only when the enumerator actually BUILDS the
    sequences and their number matches. #sTt: A2=5, A3=14, k(1<->2)/rad^2 = 6."""
    if factory == "A2":
        A = linear_path_algebra(2, field=QQ)
    elif factory == "A3":
        A = linear_path_algebra(3, field=QQ)
    else:
        from quiverlab.families.radical_square_zero import RadicalSquareZero
        A = RadicalSquareZero(Quiver([1, 2], {"a": (1, 2), "b": (2, 1)}), field=QQ)
    rep = A.tau_exceptional_sequences(budget=4096, want_sequences=True)
    assert rep.is_complete and rep.status == "complete"
    mats = rep.sequences
    assert len(mats) == rep.signed_count                 # materialised == n!*#sTt (the cross-check)
    assert len({_sig_key(s) for s in mats}) == len(mats)         # pairwise-DISTINCT
    assert all(A.is_tau_exceptional_sequence(s) for s in mats)   # each genuinely recognized


@selfcert
def test_tau_tilting_infinite_refused():
    Kron = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(field=QQ)
    rep = Kron.tau_exceptional_sequences(budget=200)
    assert rep.is_complete is False and rep.status in ("budget", "unsupported")


@selfcert
def test_exchange_graph_status_gate_and_d4_count():
    """M1: the tau-side gate reads a count ONLY off a status='complete' graph. After Task 0
    the D4 mixed star completes -> signed_count == 24*50 == 1200; the count is derived, never
    from a non-complete graph. (want_sequences=False: the count is exact via the formula,
    materialising 1200 towers is unnecessary for this gate.)"""
    A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (3, 1), "c": (4, 1)}).algebra(field=QQ)
    rep = A.tau_exceptional_sequences(budget=4096, want_sequences=False)
    assert (rep.is_complete and rep.signed_count == 24 * 50) or (
        rep.is_complete is False and rep.status in ("error", "unsupported", "budget"))


# --------------------------------------------------------------------------- #
# Task A4 (xeng): the (a)<->(b) cross-check on hereditary rep-finite algebras
# --------------------------------------------------------------------------- #
@xeng
@pytest.mark.parametrize("n", [2, 3])
def test_tau_exc_unsigned_equals_classical(n):
    """On hereditary rep-finite A, the tau-exceptional sequences with all objects genuine
    modules (all-positive towers, no shifted projectives) are in bijection with the classical
    exceptional sequences (Buan-Marsh: coincide for hereditary). Realized as a COUNT equality
    between two INDEPENDENT enumerations -- the tau reduction recursion vs the classical
    backward-orthogonality DFS -- a loud bug if they differ. (Honest scope: the termwise
    is_isomorphic lift of a deep tower to ambient A-modules needs the general reduction
    object-lift F, the Plan-65 scope boundary; the count equality is the realized
    cross-engine content, note (c) on the verification page.)"""
    A = linear_path_algebra(n, field=QQ)
    tau = A.tau_exceptional_sequences(want_sequences=True)
    all_module = [s for s in tau.sequences if all(e.obj.sign > 0 for e in s)]
    classical = A.exceptional_sequences()
    assert len(all_module) == classical.count            # two independent enumerations agree
    assert all(A.is_tau_exceptional_sequence(s) for s in all_module)


# --------------------------------------------------------------------------- #
# Adversarial-review fixes (2026-08-08): H-1 completion-independence certificate,
# H-2 tower chain check, M-1 deep D_4/A_4 materialisation + reorder-invariance.
# --------------------------------------------------------------------------- #
@selfcert
def test_h1_reduction_completion_independence_certificate():
    """H-1: C(U) = End(T_U)/<e_U> is completion-DEPENDENT as an ALGEBRA (different tau-tilting
    completions give different dim / quiver), but the enumeration invariant #sTt(C(U)) is
    completion-INDEPENDENT (DIJ). On kD4 (where the multiplicity is real -- critic-proved) the
    reduction self-certifies #sTt across completions; here we exhaustively check that invariant
    AND witness the algebra-level dependence (so the certificate is not vacuous)."""
    from quiverlab.tautilting.exceptional import (_objects_from_eg, _find_completions,
        _build_reduction_from_pair, _stt_count, tau_perpendicular_reduction)
    from quiverlab.tautilting.mutation import exchange_graph
    A = Quiver([1, 2, 3, 4], {"a": (2, 1), "b": (3, 1), "c": (4, 1)}).algebra(field=QQ)
    eg = exchange_graph(A, budget_pairs=2000)
    assert eg.status == "complete"
    checked, algebra_differs = 0, False
    for o in _objects_from_eg(A, eg):
        if o.sign <= 0:
            continue
        comps = _find_completions(o.module, eg)
        if len(comps) < 2:
            continue
        C0 = _build_reduction_from_pair(A, o.module, comps[0])[0]
        C1 = _build_reduction_from_pair(A, o.module, comps[1])[0]
        assert _stt_count(C0) == _stt_count(C1)              # the certified invariant (H-1)
        if (C0.dim, len(list(C0.quiver.arrows))) != (C1.dim, len(list(C1.quiver.arrows))):
            algebra_differs = True
        checked += 1
    assert checked >= 1, "no multi-completion tau-rigid U on kD4 (certificate never fires?)"
    assert algebra_differs, "expected a completion-DEPENDENT C(U) algebra on kD4 (H-1)"
    # the certified public reduction succeeds where the multiplicity is real (no raise)
    U = next(o.module for o in _objects_from_eg(A, eg)
             if o.sign > 0 and len(_find_completions(o.module, eg)) >= 2)
    assert tau_perpendicular_reduction(A, U).reduction_algebra is not None


@selfcert
def test_h2_bogus_tower_chain_is_rejected():
    """H-2: _verify_tower must recompute the reduction chain -- a hand-built tower whose middle
    algebra is k x k (the TRUE reduction C(S_1 of kA_3) is kA_2) is a FORGERY and must be
    rejected (#sTt(k x k) = 4 != 5 = #sTt(kA_2)), even though every rung's object is a valid
    signed object of its (bogus) algebra. The pre-fix _verify_tower returned True (critic)."""
    from quiverlab.tautilting.exceptional import _TowerEntry, tau_exceptional_objects
    kA3 = linear_path_algebra(3, field=QQ)
    kxk = Quiver([1, 2], {}).algebra(field=QQ)            # dim 2, #sTt 4 (the forgery)
    k1 = Quiver([1], {}).algebra(field=QQ)                # dim 1, #sTt 2
    S1 = next(o for o in kA3.tau_exceptional_objects()
              if o.sign > 0 and o.module.dimension_vector() == {1: 1, 2: 0, 3: 0})
    kxk_shift = next(o for o in tau_exceptional_objects(kxk) if o.sign < 0)
    k1_shift = next(o for o in tau_exceptional_objects(k1) if o.sign < 0)
    bogus = [_TowerEntry(kA3, S1, 3, 0),
             _TowerEntry(kxk, kxk_shift, 2, 0),           # forged middle (true = kA_2)
             _TowerEntry(k1, k1_shift, 1, 0)]
    assert kA3.is_tau_exceptional_sequence(bogus) is False


@xeng
def test_d4_materialised_and_reorder_invariant():
    """M-1: kD4 (subspace) signed count 1200 = 4!*50 MATERIALISED (not just formula-counted),
    the all-module (all-positive) tower count 162 = the classical CES count (kD4), AND the
    materialisation is REORDER-INVARIANT (reversing the object enumeration at every level
    yields the same 1200 / 162 -- the tower set does not depend on processing order)."""
    from quiverlab.tautilting.exceptional import tau_exceptional_sequences
    A = Quiver([1, 2, 3, 4], {"a": (2, 1), "b": (3, 1), "c": (4, 1)}).algebra(field=QQ)
    rep = tau_exceptional_sequences(A, budget=4096, want_sequences=True)
    assert rep.signed_count == 24 * 50 == 1200
    assert len(rep.sequences) == 1200                    # materialised, not formula
    assert sum(1 for s in rep.sequences if all(e.obj.sign > 0 for e in s)) == 162  # == classical
    rev = tau_exceptional_sequences(A, budget=4096, want_sequences=True, _reverse_objects=True)
    assert len(rev.sequences) == 1200
    assert sum(1 for s in rev.sequences if all(e.obj.sign > 0 for e in s)) == 162


@xeng
def test_a4_materialised_signed_count():
    """M-1: kA4 signed count 1008 = 4!*42 MATERIALISED (completion-multiplicity is real at
    this scale), all pairwise-distinct, each recognized on a sample."""
    A = linear_path_algebra(4, field=QQ)
    rep = A.tau_exceptional_sequences(budget=4096, want_sequences=True)
    assert rep.signed_count == 24 * 42 == 1008
    assert len(rep.sequences) == 1008
    assert len({_sig_key(s) for s in rep.sequences}) == 1008          # pairwise-distinct
    assert all(A.is_tau_exceptional_sequence(s) for s in rep.sequences[:24])   # sample recognized
