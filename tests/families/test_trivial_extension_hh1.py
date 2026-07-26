"""Hochschild literature-oracle batteries around trivial extensions, incidence
algebras, and truncated algebras (Plan 29, Part 2 items 4-6). tests/families/
auto-assigns to the deep bucket (tests/conftest.py).

Sources (registry keys):

* ``HH^1(T(A)) != 0`` for every finite-dimensional A -- Cibils-Marcos-Redondo-
  Solotar, Glasgow Math. J. 45 (2003) 21-40 (``cmrs_split``): Z(A) is always a
  direct summand. The HH^1 decomposition and Example 2.20 -- Cibils-Redondo-
  Saorin, J. Algebra Appl. 3 (2004) 143-159 (``crs_trivial_ext_hh1``).
* ``HH^n(incidence algebra) = SH^n(order complex)`` -- Cibils, J. Pure Appl.
  Algebra 56 (1989) 221-232 (``cibils_incidence``); Redondo, J. London Math.
  Soc. 77 (2008) 465-480 (``redondo_incidence``).
* ``HH_n(A) = 0`` for ``n >= 1`` when the quiver is acyclic -- Cibils, LNM 1177
  (1986) 55-59 (``cibils_acyclic``).
* ``dim HH^*(kQ/R^N) < infinity  <=>  Q acyclic`` -- Xu-Han-Jiang, Sci. China
  Ser. A 50 (2007) 727-736 (``xhj_truncated``).

Provenance: ``docs/plans/2026-07-25-literature-oracles-deep-research.md``
(Cibils cluster); values re-verified live here.
"""
from quiverlab import CC, GF, IncidenceAlgebra, Quiver, TrivialExtension
from quiverlab.fields import QQ

F = GF(32003)
_DIAMOND = [("b", "x"), ("b", "y"), ("x", "t"), ("y", "t")]   # 2x2 grid = comm. square


def _base_algebras():
    # QQ (a presented-eligible field per Plan-31 D3: QQ/GF(p) coefficients are
    # string-representable for Quiver.algebra, so T(A) takes the genuine
    # presented route). Char 0, so the CMRS HH^1 values are unchanged from the
    # former CC bases; CC-algebraic bases fall back to the structure-constant
    # build (a distinct D3 branch, guarded in test_trivial_extension_presented.py).
    return {
        "kA2": Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=QQ),
        "kA3": Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=[], field=QQ),
        "2-Kronecker": Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=QQ),
        "comm-square": IncidenceAlgebra(_DIAMOND, field=QQ),
    }


# ---------------------------------------------- trivial-extension HH^1 (CMRS)

def test_trivial_extension_HH1_never_vanishes():
    """HH^1(T(A)) != 0 for every finite-dimensional A (CMRS 2003, ``cmrs_split``):
    Z(A) is always a summand. Since Plan 31 T(A) over QQ carries a genuine quiver
    presentation; over QQ the auto engine is still the normalized-bar oracle (QQ
    is not a prime field). The HH dims are iso-invariant (unchanged from the
    retained structure-constant build). Characteristic-independent; pinned here
    over char 0."""
    for name, A in _base_algebras().items():
        TA = TrivialExtension(A)
        assert TA.quiver is not None                   # presented route (Plan 31)
        hh1 = TA.hochschild_cohomology(1).dims[1]
        assert hh1 >= 1, f"{name}: HH^1(T(A)) must be nonzero"


def test_cs_now_computes_on_presented_trivial_extension():
    """Plan 31 STRENGTHENING (was: CS refuses the presentation-free T(A) with a
    ValueError). The presented T(kA_2) carries a quiver + relations -- exactly what
    the Chouhy-Solotar engine needs -- so ``engine="cs"`` now COMPUTES and agrees
    with the normalized-bar oracle. Over GF(32003) the auto engine is the fast
    bar-rank engine (same bar complex; identical dims): [3, 1, 1] at top 2.
    (``chouhy_solotar``; ``bar``; ``assem_book``.) The fuller cross-engine +
    iso-invariance battery is in tests/families/test_trivial_extension_presented.py."""
    TA = TrivialExtension(Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(32003)))
    assert TA.quiver is not None
    cs = TA.hochschild_cohomology(2, engine="cs").dims
    bar = TA.hochschild_cohomology(2).dims                           # auto = fast bar-rank over GF(p)
    assert cs == bar == [3, 1, 1]


def test_crs2004_example_2_20_cyclic_hh1_finding():
    """Cibils-Redondo-Saorin (2004), Example 2.20 (``crs_trivial_ext_hh1``): the
    Z_5 cycle 1->2->3->4->5->1 with the two length-3 monomial relations the paper
    writes RIGHT-TO-LEFT as a4 a3 a2 and a3 a2 a1. Translating to quiverlab's
    LEFT-TO-RIGHT composition: a4 a3 a2 = a2*a3*a4 and a3 a2 a1 = a1*a2*a3.

    FINDING (Plan 29): the secondary (ar5iv) transcription records HH^1(A) = 0 in
    characteristic 0, but quiverlab's validated bar oracle robustly gives
    dim HH^1(A) = 1 (over CC and GF(32003) alike) with dim Z(A) = 2. The extra
    central element is the SINGLE surviving oriented 5-cycle a3*a4*a5*a1*a2 -- the
    only length-5 cyclic path none of whose length-3 windows is a killed
    relation; the four other rotations each contain a1*a2*a3 or a2*a3*a4 and die.
    No 2-relation Z_5 construction (adjacent or gapped, either orientation) yields
    HH^1 = 0. The Cibils cluster report flags this family's per-degree values as
    transcription-lossy, so the paper's 0 is NOT frozen as a strict pin;
    quiverlab's computed value is pinned instead as a regression guard."""
    Q = Quiver([1, 2, 3, 4, 5],
               {"a1": (1, 2), "a2": (2, 3), "a3": (3, 4), "a4": (4, 5), "a5": (5, 1)})
    A = Q.algebra(relations=["a1*a2*a3", "a2*a3*a4"], field=CC)
    assert A.dim == 21
    assert A.center()[0] == 2                          # k.1 + the surviving 5-cycle
    assert A.hochschild_cohomology(1).dims == [2, 1]   # quiverlab's verified value (paper: 0)


# ------------------------------------------ incidence cohomology = nerve (SH)

def test_incidence_cohomology_crown_is_S1():
    """Cibils (1989) / Redondo (2008) (``cibils_incidence`` / ``redondo_incidence``):
    HH^n of an incidence algebra is the simplicial cohomology of the poset order
    complex. The crown (two minimal a,a'; two maximal b,b'; every a_i < b_j) has
    order complex = the 4-cycle a-b-a'-b'-a, homotopy-equivalent to S^1, so
    HH^* = [1,1,0,0] (derived from the theorem). It is hereditary -- no parallel
    paths, hence no commutativity relations -- so HH^{>=2} = 0 and top=2 exhibits
    both the S^1 pattern and its vanishing tail. (top=3 is infeasible: the bar
    differential d^3 is 19208 x 2744, past the guard.)"""
    crown = IncidenceAlgebra([("a", "b"), ("a", "bp"), ("ap", "b"), ("ap", "bp")], field=CC)
    assert crown.relations == []                       # no parallel paths -> hereditary
    assert crown.hochschild_cohomology(2).dims == [1, 1, 0]   # SH^0 = SH^1 = k (S^1)


def test_incidence_cohomology_contractible_chain_vanishes():
    """A chain poset 1<2<3 has a contractible order complex (a 2-simplex), so
    HH^{>=1} = 0: HH^* = [1,0,0,0] (Cibils 1989 / Redondo 2008). Contrast the
    crown's HH^1 = 1: the nerve topology is read off HH^1."""
    chain = IncidenceAlgebra([(1, 2), (2, 3)], field=CC)      # = kA_3
    assert chain.hochschild_cohomology(3).dims == [1, 0, 0, 0]


# ------------------------------------------- acyclic HH_* vanishing (Cibils 1986)

def test_acyclic_quiver_hochschild_homology_vanishes():
    """Cibils (1986) (``cibils_acyclic``): Q acyclic => HH_n(A) = 0 for n >= 1,
    with HH_0 concentrated in degree 0 (dim = #vertices here). Checked on kA_3 and
    the commutative square. (For the comm. square top=2: the homology boundary
    b_4 at top=3 is 4608 x 36864 = 170M cells, past the guard; HH_1 = HH_2 = 0
    already witnesses the vanishing.)"""
    A3 = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=[], field=F)
    assert A3.hochschild_homology(3).dims == [3, 0, 0, 0]
    csq = IncidenceAlgebra(_DIAMOND, field=F)
    assert csq.hochschild_homology(2).dims == [4, 0, 0]


# ---------------------------------- truncated finiteness boolean (Xu-Han-Jiang)

def test_truncated_finiteness_boolean_acyclic_vs_cyclic():
    """Xu-Han-Jiang (2007) Thm 3 (``xhj_truncated``): for a truncated algebra
    kQ/R^N, dim_k HH^*(A) < infinity  <=>  Q has no oriented cycle. A boolean
    oracle (no exact counts): an acyclic truncated algebra has a vanishing HH^*
    tail; a cyclic one never vanishes.

    The acyclic witness is linear A_3 with rad^2 = 0 (no oriented cycle),
    HH^{>=1} = 0. For the cyclic witness the plan named kZ_3/J^2, but over a field
    of characteristic != 3 that algebra's HH^* is sparse (HH^2 = HH^3 = HH^4 = 0;
    the next nonzero degrees sit near degree 6, past the bar guard). kZ_2/J^2 is
    the clean cyclic witness at feasible depth: HH^n = 1 for EVERY n, so it never
    vanishes -- exhibited here at degrees 3 and 4."""
    # acyclic: kQ/R^2 with no oriented cycle -> finite (in fact trivial) total HH^*
    acyc = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(relations=["a*b"], field=F)
    hh = acyc.hochschild_cohomology(3).dims
    assert hh[2] == 0 and hh[3] == 0                   # vanishing tail (gl.dim 2)
    # cyclic: kZ_2/J^2, an oriented 2-cycle with rad^2 = 0 -> HH^* never vanishes
    cyc = Quiver([0, 1], {"a": (0, 1), "b": (1, 0)}).algebra(relations=["a*b", "b*a"], field=F)
    hh = cyc.hochschild_cohomology(4).dims
    assert hh[3] > 0 and hh[4] > 0                     # nonzero past degree 2 (infinite total HH^*)
