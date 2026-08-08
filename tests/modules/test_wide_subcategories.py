"""The wide-subcategory poset via the core label order (Plan 64 / R26, Enomoto 2201.00595).
Literature: #wide(kA2) = 5 (poset M3: three incomparable atoms = the 3 bricks); #wide(kA3) =
14 (poset NC(A3)). Self-cert: the kappa order and the core label order coincide (Enomoto).
DISCRIMINATOR (ruling 1/2): #wide == #torsion iff A is representation-FINITE (Marks-Stovicek),
so on the rep-finite kA_n the COUNT never discriminates -- the poset STRUCTURE (M3, NC(A3))
does. The strict-inequality #wide < #torsion oracle needs a rep-INFINITE tau-tilting-finite
algebra and is DEFERRED to implementation (see the comment block at the end). QQ-scope."""
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.tautilting.congruence import torsion_lattice, wide_subcategories

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@lit
@pytest.mark.parametrize("n, nwide", [(2, 5), (3, 14)])   # kA_n rep-finite => #wide == #torsion
def test_wide_count_type_A(n, nwide):
    # NOTE: on kA_n the count is Catalan(n+1) AND equals #torsion, because kA_n is
    # representation-FINITE (Marks-Stovicek: #wide == #torsion iff rep-finite). The count is
    # therefore NON-discriminating here; the poset STRUCTURE (below) is the arbiter.
    W = wide_subcategories(linear_path_algebra(n, field=QQ))
    assert W.is_complete and W.size == nwide


@selfcert
def test_kappa_order_equals_core_label_order():
    # Enomoto: the two poset constructions coincide (and both = wide subcategories).
    W = wide_subcategories(linear_path_algebra(2, field=QQ))
    assert W.constructions_agree is True


# --- poset helpers: W.order[e] = frozenset of elements <= e (down-set incl. e), matching the
#     TorsionLattice.order convention (id -> frozenset(ids <= it)). --------------------------
def _strict_below(W, e):
    return set(W.order[e]) - {e}

def _bottom(W):
    b = [e for e in W.elements if _strict_below(W, e) == set()]
    assert len(b) == 1, ("unique bottom", b)
    return b[0]

def _top(W):
    t = [e for e in W.elements if set(W.order[e]) == set(W.elements)]
    assert len(t) == 1, ("unique top", t)
    return t[0]

def _atoms(W, bot):
    return [e for e in W.elements if _strict_below(W, e) == {bot}]

def _lower_covers(W, x):
    below = _strict_below(W, x)
    return [c for c in below if not any(c in _strict_below(W, d) for d in below if d != c)]

def _incomparable(W, a, b):
    return a not in W.order[b] and b not in W.order[a]


@lit
def test_kA2_wide_poset_is_M3():
    # M3 (== NC(A2)): bottom 0, THREE pairwise-incomparable atoms (add S1, add S2, add P1),
    # top mod A; each atom covers bottom and is a lower cover of top. Ruling 2: assert the
    # SHAPE, not just the count -- and actually USE the derived bottom/top/atoms.
    W = wide_subcategories(Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ))
    assert W.size == 5 and len(W.elements) == 5 and len(W.labels) == 5
    bot, top = _bottom(W), _top(W)
    atoms = _atoms(W, bot)
    assert len(atoms) == 3                                   # exactly three atoms
    for i in range(3):                                       # pairwise incomparable
        for j in range(i + 1, 3):
            assert _incomparable(W, atoms[i], atoms[j])
    for a in atoms:                                          # each atom covers bottom
        assert _strict_below(W, a) == {bot}
    assert set(_lower_covers(W, top)) == set(atoms)          # top covers EXACTLY the 3 atoms
    # rank (= |label set| = #simple objects of the wide subcat) histogram of M3 is 1,3,1
    assert sorted(len(lab) for lab in W.labels) == [0, 1, 1, 1, 2]
    assert sum(1 for lab in W.labels if len(lab) == 1) == 3  # 3 singleton-brick atoms


@lit
def test_kA3_wide_poset_is_NC_A3():
    # NC(A3) == NC of [4]: 14 elements; rank (= |label set|) Whitney numbers 1,6,6,1 (Narayana
    # N(4,k)); 6 atoms each covering bottom, 6 coatoms each covered by top; the poset is
    # self-dual. Ruling 2: pin the STRUCTURE beyond the count 14.
    from collections import Counter
    W = wide_subcategories(linear_path_algebra(3, field=QQ))
    assert W.is_complete and W.size == 14 and len(W.elements) == 14
    bot, top = _bottom(W), _top(W)
    assert dict(Counter(len(lab) for lab in W.labels)) == {0: 1, 1: 6, 2: 6, 3: 1}
    assert len(_atoms(W, bot)) == 6                          # 6 atoms (single-brick wides)
    assert len(_lower_covers(W, top)) == 6                   # 6 coatoms (top's lower covers)
    # self-duality of NC(4): down-set-size and up-set-size histograms coincide
    downs = Counter(len(W.order[e]) for e in W.elements)
    ups = Counter(sum(1 for f in W.elements if e in W.order[f]) for e in W.elements)
    assert downs == ups
    # RIGOROUS OPTION: an order-isomorphism check of W.order against the hardcoded NC(4)
    # covering table DERIVED in the plan appendix (14 partitions of {1,2,3,4}, 28 covers) is
    # the definitive arbiter; add it at implementation if this Whitney/atom/coatom/self-dual
    # battery is ever insufficient to distinguish NC(4) from another 14-element graded poset.


# --- DEFERRED (ruling 1, path b): the strict-inequality #wide < #torsion discriminator. -----
# A rep-FINITE algebra can NEVER witness #wide < #torsion (Marks-Stovicek: #wide == #torsion
# iff rep-finite; wide -> tors is always injective). So the prior draft's candidates
# (kZ2/rad^2, the rad^2 Nakayama 1->2->3) are GUARANTEED #wide == #torsion and were REMOVED --
# the assert `W.size < len(L.elements)` on them would ALWAYS FAIL. A genuine witness needs a
# representation-INFINITE but tau-tilting-finite algebra: the TAME preprojective algebras
# Pi(D4) (#torsion = |W(D4)| = 192) or Pi(A5) (720). At authoring the shipped P45 engine could
# NOT compute either cheaply -- Pi(D4)/QQ: exchange_graph did not reach 50 vertices in 120s;
# Pi(D4)/GF(31): mutate() raised QuiverlabError deep in the BFS (status="error" after 6
# vertices). So this oracle is DEFERRED to implementation:
#     A = PreprojectiveAlgebra("D4", field=QQ)          # or field=GF(p), p > dim = 28, for speed
#     L, W = torsion_lattice(A, budget=256), wide_subcategories(A, budget=256)   # budget >= 192
#     if L.is_complete and W.is_complete:
#         assert W.size < len(L.elements)               # 192 = #torsion; #wide < 192 is a THEOREM
# marked oracle_literature (the strict inequality is Marks-Stovicek, not an empirical guess). If
# the engine still refuses Pi(D4)/Pi(A5), this stays a documented honest-scope gap and the
# M3/NC structural arbiters above carry the wide-vs-torsion discrimination.
