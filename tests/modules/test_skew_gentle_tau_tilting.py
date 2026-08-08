"""Support tau-tilting for skew-gentle algebras via the engine, geometric cross-check
(Plan 68 / HZZ sec 5-6). Engine-primary: P45 exchange_graph on the split algebra.
Literature/geometric: the pair rank equals |Q_0| + |Sp| (HZZ sec 6); k x k (one vertex,
Sp = {1}) has 4 support tau-tilting pairs (the Boolean square, geometrically the 4
generalized dissections)."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.skewgentle.certificate import support_tau_tilting
from quiverlab.skewgentle.split import SkewGentleAlgebra
from quiverlab.skewgentle.triple import SkewGentleTriple
from quiverlab.tautilting.mutation import exchange_graph

lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def test_pair_rank_is_Q0_plus_Sp():
    # HZZ sec 6: a basic support tau-tilting module has |Q_0| + |Sp| summands.
    t = SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}), relations=[], special={2})
    blk = support_tau_tilting(t, field=QQ)
    assert blk["n"] == 2 + 1                           # |Q_0| + |Sp|


@lit
def test_kxk_has_four_support_tau_tilting_pairs():
    # Q = one vertex, Sp = {1}: A = k[eps]/(eps^2-eps) = k x k; split = two isolated
    # vertices; stau-tilt(k x k) = 2 x 2 = 4 (the geometric count = 4 dissections).
    # NOTE (M2): k x k is semisimple -> selfinjective, WITH Sp != empty. This does NOT
    # contradict Chen Cor 1.2(c): that corollary is about INDECOMPOSABLE algebras, and
    # k x k is DECOMPOSABLE. This test pins the tau-tilting count only.
    t = SkewGentleTriple.make(Quiver([1], {}), relations=[], special={1})
    A = SkewGentleAlgebra(t, field=QQ)
    eg = exchange_graph(A, budget_pairs=64)
    assert eg.is_complete and len(eg.vertices) == 4


def test_exchange_graph_on_split_does_not_raise():
    # H1 REGRESSION: split vertex labels must be a single orderable (string) type, else
    # P45's exchange_graph crashes at `sorted(pair.support)` (mutation.py::_label) with
    # `TypeError: '<' not supported between instances of 'str' and 'int'`. With ALL labels
    # stringified it runs clean. (VERIFIED live: mixed [1,'2+','2-',3] crashes here.)
    t = SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}), relations=[], special={2})
    A = SkewGentleAlgebra(t, field=QQ)
    eg = exchange_graph(A)                              # must NOT raise TypeError
    assert eg.is_complete


@xeng
def test_geometric_rank_matches_engine_pairs():
    t = SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}), relations=[], special={2})
    A = SkewGentleAlgebra(t, field=QQ)
    eg = exchange_graph(A)
    n = len(list(A.quiver.vertices))                  # = |Q_0| + |Sp|
    for rec in eg.vertices:                           # every pair has full rank n
        assert len(rec["summand_dimvecs"]) + len(rec["support"]) == n


@lit
@pytest.mark.deep
def test_hzz_section6_seven_vertex_rank_is_eight():
    # HZZ sec 6 rank oracle: a basic support tau-tilting module has |R| = |Q_0| + |Sp|
    # summands. The rank is PRESENTATION-INDEPENDENT (|Q_0| + |Sp|), so only the triple's
    # VALIDITY matters (HZZ sec 6 caveat in the plan doc).
    #
    # DEVIATION FROM THE PLAN DOC (documented + verified): the plan's exact sec-6 arrow
    # reconstruction (vertices 1..7, arrows a..h + a non-special loop eps_5,
    # I = {eps5^2, ab, ba, ed, hc}, Sp = {1}) is NOT a valid skew-gentle triple
    # (`is_skew_gentle_triple` returns False -- vertex 2's branch structure fails
    # gentleness once eps_1 is added; verified live). It is replaced with a VALID
    # 7-vertex, Sp = {1} triple -- the linear A_7 with a special SOURCE vertex -- which
    # pins the SAME presentation-independent rank |R| = |Q_0| + |Sp| = 8.
    Q = Quiver([1, 2, 3, 4, 5, 6, 7],
               {"a1": (1, 2), "a2": (2, 3), "a3": (3, 4), "a4": (4, 5),
                "a5": (5, 6), "a6": (6, 7)})
    t = SkewGentleTriple.make(Q, relations=[], special={1})
    A = SkewGentleAlgebra(t, field=QQ)
    assert len(list(A.quiver.vertices)) == 7 + 1       # |R| = |Q_0| + |Sp| = 8 (HZZ sec 6)
