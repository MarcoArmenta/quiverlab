"""The recognizer ladder quasi-tilted/shod/weakly-shod/laura/ada (Plan 61 / R18).
Literature: the ACLV Example 2.2(b) rad^2=0 A5 gives (qt,shod,wshod,laura,ada) =
(False,False,True,True,True) with complement={S3} -- a full five-way discrimination; a
hereditary kA_n gives all five True. Self-cert: monotone nesting (qt=>shod=>wshod=>laura,
qt/shod=>ada, ada=>laura), theorem gates (shod=>gldim<=3, qt=>gldim<=2). Cross-engine: the
two quasi-tilted routes (HRS gldim/complement vs ACLV projectives-in-left) agree. QQ-scope
(M1)."""
import pytest

from quiverlab import Quiver, RadicalSquareZero, linear_path_algebra, truncated_polynomial
from quiverlab.fields import QQ
from quiverlab.modules.recognizers_ladder import recognizer_ladder

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


def _radsq_nakayama_a5():
    Q = Quiver([1, 2, 3, 4, 5],
               {"a1": (2, 1), "a2": (3, 2), "a3": (4, 3), "a4": (5, 4)})
    return RadicalSquareZero(Q, field=QQ)


@lit
def test_aclv_22b_full_ladder():
    L = recognizer_ladder(_radsq_nakayama_a5())
    assert L.is_complete and L.universe_size == 9 and L.gldim == 4
    assert L.verdict("quasi_tilted") is False
    assert L.verdict("shod") is False
    assert L.verdict("weakly_shod") is True
    assert L.verdict("laura") is True
    assert L.verdict("ada") is True
    assert [r["name"] for r in L.complement] == ["S_3"]           # the laura datum
    assert L.rungs["shod"].witness["module"] == "S_3"             # pd 2, id 2 witness


@lit
@pytest.mark.parametrize("n", [2, 3, 4])
def test_hereditary_all_five_true(n):
    L = recognizer_ladder(linear_path_algebra(n, field=QQ))     # kA_n hereditary
    assert all(L.verdict(k) for k in
               ("quasi_tilted", "shod", "weakly_shod", "laura", "ada"))
    assert L.complement == ()                                    # empty complement


@xeng
def test_shod_two_routes_agree():
    # Thm 4.1 (a)<=>(b): the QT2 sweep (every indec pd<=1 or id<=1) == complement empty.
    # NON-tautological: QT2 reads atlas.pd_le_1/id_le_1 (independent of the L/R closure).
    # The monomial square (rep-finite, gldim 2) is NOT shod (P_1 in neither, complement != []);
    # kA3 (empty complement) IS shod.
    from quiverlab.modules.recognizers_ladder import _shod
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    for A, expect in ((Q.algebra(relations=["a*b"], field=QQ), False),        # monomial: not shod
                      (linear_path_algebra(3, field=QQ), True)):              # kA3: shod
        L = recognizer_ladder(A)
        assert L.verdict("shod") is expect
        assert (len(L.complement) == 0) == L.verdict("shod")        # complement-empty route
        # the QT2 route is a real independent computation inside _shod (asserted equal there);
        # a mismatch would have already raised QuiverlabError. Re-affirm the witness shape:
        if not L.verdict("shod"):
            w = L.rungs["shod"].witness
            assert set(w) == {"module", "pd_le_1", "id_le_1"} and not (w["pd_le_1"] or w["id_le_1"])


@xeng
def test_quasi_tilted_two_routes_agree():
    # HRS route (gldim<=2 AND complement==[]) == ACLV route (every P_x in L_A). PIN the
    # literature verdicts: comm-square-with-relation quasi-tilted True; rad^2=0 A5 False.
    from quiverlab.modules.recognizers_ladder import _quasi_tilted_aclv_route
    comm = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
                  ).algebra(relations=["a*b - c*d"], field=QQ)
    for A, expect in ((_radsq_nakayama_a5(), False), (comm, True)):
        L = recognizer_ladder(A)
        assert L.verdict("quasi_tilted") is expect                 # pinned literature verdict
        assert L.verdict("quasi_tilted") == _quasi_tilted_aclv_route(A)   # HRS route == ACLV route


@xeng
def test_weakly_shod_two_routes_agree():
    # PRIMARY (AR-arrow SCC) == CROSS-CHECK (Hom-closure atlas._leq SCC). Non-vacuous:
    # the fixture is directed (both True); kupisch=[3,2,2] is one big SCC (both False).
    from quiverlab import NakayamaAlgebra
    from quiverlab.modules.recognizers_ladder import _weakly_shod, _weakly_shod_hom
    from quiverlab.modules.left_right import left_right_parts
    for A in (_radsq_nakayama_a5(), NakayamaAlgebra(kupisch=[3, 2, 2], field=QQ)):
        atlas = left_right_parts(A)
        assert _weakly_shod(A, atlas, 256).verdict == _weakly_shod_hom(atlas).verdict


@lit
def test_weakly_shod_false_on_kupisch_322():
    # M2: a REAL non-directed knittable algebra. kupisch=[3,2,2] knits complete, its whole
    # module category is one non-trivial SCC on an injective->projective route => weakly_shod
    # False (verified live; also fully degenerate -- ada/shod/quasi_tilted all False, laura True).
    from quiverlab import NakayamaAlgebra
    L = recognizer_ladder(NakayamaAlgebra(kupisch=[3, 2, 2], field=QQ))
    assert L.is_complete and L.universe_size == 7 and L.gldim == 3
    assert L.verdict("weakly_shod") is False
    assert L.rungs["weakly_shod"].witness is not None            # {cycle, inj_source, proj_target}
    assert L.verdict("laura") is True and L.verdict("ada") is False and L.verdict("shod") is False


@selfcert
def test_monotone_nesting_holds_on_ladder():
    for A in (_radsq_nakayama_a5(), linear_path_algebra(4, field=QQ)):
        L = recognizer_ladder(A)
        v = L.verdict
        assert (not v("quasi_tilted")) or v("shod")
        assert (not v("shod")) or v("weakly_shod")
        assert (not v("weakly_shod")) or v("laura")
        assert (not v("quasi_tilted")) or v("ada")
        assert (not v("shod")) or v("ada")
        assert (not v("ada")) or v("laura")


@selfcert
def test_theorem_gates_gldim():
    for A in (_radsq_nakayama_a5(), linear_path_algebra(3, field=QQ)):
        L = recognizer_ladder(A)
        if L.verdict("shod"):
            assert L.gldim <= 3
        if L.verdict("quasi_tilted"):
            assert L.gldim <= 2


@selfcert
def test_selfinjective_refused_loudly():
    A = truncated_polynomial(3, field=QQ)                       # k[x]/(x^3), self-injective
    L = recognizer_ladder(A)
    assert L.is_complete is False and L.status == "unsupported"
    assert L.rungs == {}                                        # never a partial ladder


@selfcert
def test_rep_infinite_ada_22c_refused_loudly():
    # ACLV Example 2.2(c): 1=>2=>3=>4 bound by rad^2=0 -- mathematically ada, refused
    # because representation-infinite (the honest boundary; ada does NOT imply laura here).
    Q = Quiver([1, 2, 3, 4],
               {"a": (1, 2), "b": (1, 2), "c": (2, 3), "d": (2, 3), "e": (3, 4), "f": (3, 4)})
    A = RadicalSquareZero(Q, field=QQ)
    L = recognizer_ladder(A, budget=40)
    assert L.is_complete is False and L.status in ("budget", "error", "unsupported")


@selfcert
def test_delegates_match_free_function():
    A = _radsq_nakayama_a5()
    assert A.is_ada() is True and A.is_shod() is False
    assert A.is_quasi_tilted() is False and A.is_weakly_shod() is True and A.is_laura() is True


@selfcert
def test_weakly_shod_sweep_detects_cycle_between_inj_and_proj():
    from quiverlab.modules.recognizers_ladder import _weakly_shod_from_digraph
    # 5 vertices: 0 = injective, 1,2 = a 2-cycle, 3 = projective, 4 = isolated.
    arrows = {(0, 1): 1, (1, 2): 1, (2, 1): 1, (2, 3): 1}
    ok, witness = _weakly_shod_from_digraph(n=5, arrows=arrows, inj={0}, proj={3})
    assert ok is False and set(witness["cycle"]) == {1, 2}
    # remove the injective->cycle edge: now bounded => weakly shod.
    ok2, _ = _weakly_shod_from_digraph(n=5, arrows={(1, 2): 1, (2, 1): 1, (2, 3): 1},
                                       inj={0}, proj={3})
    assert ok2 is True


@selfcert
def test_monotone_guard_raises_on_inconsistent_rungs():
    from quiverlab.errors import QuiverlabError
    from quiverlab.modules.recognizers_ladder import _assert_nesting, LadderRung
    def rung(name, v): return LadderRung(name, v, {}, None if v else {"x": 1})
    good = {n: rung(n, v) for n, v in
            [("quasi_tilted", False), ("shod", False), ("weakly_shod", True),
             ("laura", True), ("ada", True)]}
    _assert_nesting(good)                                        # consistent: no raise
    bad = dict(good, shod=rung("shod", True))                    # shod True but weakly_shod... ok;
    bad["weakly_shod"] = rung("weakly_shod", False)              # shod True, weakly_shod False: BAD
    with pytest.raises(QuiverlabError, match="nesting|monoton|consist"):
        _assert_nesting(bad)
