"""Plan 29 (Part 3) literature-oracle battery -- CS-engine HH/HC value pins.

Companion to ``test_battery_literature.py``: pins ``cs_(co)homology_dims`` and the
Plan-20/21 native cup against values that EXIST OUTSIDE the library (published
Hochschild computations), so a regression corrupting BOTH CS and a second engine
identically would still be caught.  Every pin carries its provenance + the
verification-status label from ``docs/plans/2026-07-25-literature-oracles-deep-
research.md`` (the labels there are BINDING).

Batteries here (plan Part 3, items 1/2/6):

  1a. Bergh-Erdmann general quantum complete intersection -- CC ONLY (the GF(p)
      root-of-unity trap is real, see the module note below).
  1b. Redondo-Roman triangular string A_n (monomial, char-independent); the A_3
      degree-3 revival is cross-checked in-test against an INDEPENDENT deep engine
      (the minimal A^e engine, Plan 16) BEFORE the literature vector is pinned.
  1c. Redondo-Roman Thm 5.2 cup-triviality on A_2 (one decisive instance; the
      RR-2018 cup-nonvanishing case is deferred, see the honest-scope skip).

===========================================================================
CHAR CAVEAT (Bergh-Erdmann, verbatim from the research doc, Solotar cluster):
  A = k<x,y>/(x^a, y^b, yx - q*xy) is a codimension-2 QCI stated for "q not a root
  of unity", which needs an INFINITE field -> testable ONLY over CC / char 0.  Over
  any GF(p), q is always a root of unity (a different regime, with infinite-
  dimensional HH^*), so the [2,2,1,0,...] / [a+b-1, a+b-2, ...] vectors are NOT run
  over GF(p).  The char-p homology branches (Thm 3.1) need an infinite field of
  char p, which quiverlab lacks -- not testable, noted honestly here and on
  docs/verification.md.

DEVIATION (A_3 cross-check engine): the research doc recommends a "bar cross-check"
of the A_3/A_5 revival at pin time.  For A_3 = tri_string(3) (dim 16) the normalized
bar / fast-rank cochain complex exceeds max_cells at degree 3 (C^3 ~ 54000 wide),
so the *bar* oracle cannot reach the revival.  The independent cross-check is
therefore the minimal A^e engine (Plan 16, ``minimal_cohomology_dims`` -- iterated
A^e syzygies over GF(p), disjoint from both CS and the bar), which computes A_3 =
[1,3,0,2,0] in ~1.5s.  This honors the plan's intent (an independent engine anchors
the revival before the literature vector is pinned).
===========================================================================
"""
import pytest

from quiverlab import Quiver, CC, GF
from quiverlab.resolutions_cs.homology import cs_cohomology_dims, cs_homology_dims
from quiverlab.resolutions_cs.comparison import Comparison
from quiverlab.engine.adapter import to_engine
from quiverlab.engine.resolutions_minimal import minimal_cohomology_dims

pytest.importorskip("quiverlab.groebner")

pytestmark = [pytest.mark.oracle_literature]


# --------------------------------------------------------------------------- #
# builders                                                                     #
# --------------------------------------------------------------------------- #
def _qci(a, b, field):
    """Bergh-Erdmann codim-2 QCI  k<x,y>/(x^a, y^b, yx - 2*xy)  (q = 2, not a root
    of unity over CC).  dim A = a*b.  Mirrors the existing a=b=2 pin
    (``test_bgms_quantum_ci_homology``)."""
    return Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=[f"x^{a}", f"y^{b}", "y*x - 2*x*y"], field=field)


def _tri_string(n, field):
    """Redondo-Roman (2014) Example 3: A_n = kQ/I with vertices 0..n, two parallel
    arrows a_i, b_i : (i-1) -> i, and I = <a_i a_{i+1}, b_i b_{i+1}> (length-2
    monomial; mixed products survive).  Monomial, gentle, triangular."""
    arrows = {}
    for i in range(1, n + 1):
        arrows[f"a{i}"] = (i - 1, i)
        arrows[f"b{i}"] = (i - 1, i)
    rels = [f"a{i}*a{i+1}" for i in range(1, n)] + [f"b{i}*b{i+1}" for i in range(1, n)]
    return Quiver(list(range(n + 1)), arrows).algebra(relations=rels, field=field)


# =========================================================================== #
# 1a.  Bergh-Erdmann general (a,b) quantum complete intersection  (CC ONLY)     #
# =========================================================================== #
# Citation: P. A. Bergh, K. Erdmann, "Homology and cohomology of quantum complete
# intersections," Algebra & Number Theory 2 (2008), no. 5, 501-522.
# DOI 10.2140/ant.2008.2.501; arXiv:0709.3029.  Thm 3.2 (cohomology), Thm 3.1
# (homology).  Verification status: fetched-source (Thms 3.1/3.2 verbatim).
@pytest.mark.parametrize("a,b", [(2, 3), (3, 3)])
def test_bergh_erdmann_qci_general_ab(a, b):
    """Thm 3.2: dim HH^* = [2,2,1,0,...] for EVERY a,b>1 (independent of a,b and of
    char).  Thm 3.1: dim HH_0 = a+b-1 and dim HH_n = a+b-2 (n>=1) when char divides
    neither a nor b -- so over CC (char 0), HH_* = [a+b-1] + [a+b-2]*(...).  Extends
    the existing a=b=2 pin ([3,2,2,...]) to the family.  CC only (q=2 not a root of
    unity -- see the module CHAR CAVEAT)."""
    A = _qci(a, b, CC)
    assert cs_cohomology_dims(A, 6).dims == [2, 2, 1, 0, 0, 0, 0]
    assert cs_homology_dims(A, 6).dims == [a + b - 1] + [a + b - 2] * 6


# --- Plan 33: the (a,b) QCI oracle pushed to degree 8 at higher dimension ----
@pytest.mark.parametrize("a,b", [(2, 4), (3, 4), (4, 4), (2, 5), (5, 5)])
def test_bergh_erdmann_qci_scale_deg8(a, b):
    """Plan 33 SCALE of ``test_bergh_erdmann_qci_general_ab``: the same Bergh-Erdmann
    Thm 3.2 dim HH^* = [2,2,1,0,...] and Thm 3.1 dim HH_* = [a+b-1] + [a+b-2]*(...)
    pushed to DEGREE 8 over the enlarged family (a,b) in {(2,4),(3,4),(4,4),(2,5),
    (5,5)}, dim A = a*b up to 25 at (5,5).  The [2,2,1] cohomology plateau is
    independent of (a,b) and of degree; the homology plateau tracks a+b-2 -- the
    Bergh-Erdmann signature at scale (``qci_hh_oracle``, DOI 10.2140/ant.2008.2.501,
    Thms 3.1/3.2).  CC ONLY (q = 2 is not a root of unity -- see the module CHAR
    CAVEAT).  Recomputed live (Plan 33) before pinning; the CS engine handles the
    dim-25 (5,5) case in ~13 s/direction (the deep-bucket ceiling here)."""
    A = _qci(a, b, CC)
    assert cs_cohomology_dims(A, 8).dims == [2, 2, 1, 0, 0, 0, 0, 0, 0]
    assert cs_homology_dims(A, 8).dims == [a + b - 1] + [a + b - 2] * 8


# =========================================================================== #
# 1b.  Redondo-Roman triangular string A_n family  (all char)                   #
# =========================================================================== #
# Citation: M. J. Redondo, L. Roman, "Hochschild cohomology of triangular string
# algebras and its ring structure," J. Pure Appl. Algebra 218 (2014), no. 5,
# 925-936.  arXiv:1301.0516.  Example 3 (explicit A_n dims), Thm 4.3 (general
# combinatorial count => char-independent), Thm 5.2 (cup).  Verification status:
# fetched-source (Thm 4.3, Ex. 3, Thm 5.2 verbatim; A_1 Euler-verified).
@pytest.mark.parametrize("n,top,expected", [
    (1, 3, [1, 3, 0, 0]),      # 2-Kronecker (no relations); Euler chi = 2 - 4 = -2
    (2, 3, [1, 2, 0, 0]),      # A_{2m}: [1, 2m, 0, ...]
    (4, 4, [1, 4, 0, 0, 0]),   # A_{2m}: [1, 2m, 0, ...]
])
@pytest.mark.parametrize("fld", [GF(2), GF(3), CC], ids=["GF2", "GF3", "CC"])
def test_tri_string_even_and_kronecker(n, top, expected, fld):
    """A_1 [1,3,0,0], A_2 [1,2,0,0], A_4 [1,4,0,0,0] over GF(2), GF(3), CC (Thm 4.3
    is char-independent)."""
    A = _tri_string(n, fld)
    assert cs_cohomology_dims(A, top).dims == expected


@pytest.mark.oracle_crossengine
@pytest.mark.parametrize("p", [2, 3])
def test_tri_string_A3_revival_with_minimal_crosscheck(p):
    """A_3 = [1, 3, 0, 2, 0] over GF(p): HH^2 vanishes then HH^3 = 2 REVIVES (a
    non-monotone discriminator; A_{2m+1} has HH^{2m+1}=2).  The revival value is
    cross-checked in-test against the minimal A^e engine (Plan 16 -- iterated A^e
    syzygies, INDEPENDENT of both CS and the bar) BEFORE the literature vector is
    pinned, per the research doc's anchoring recommendation.  (The bar oracle
    exceeds max_cells at degree 3 for this dim-16 algebra -- see the module
    DEVIATION note.)"""
    A = _tri_string(3, GF(p))
    cs = cs_cohomology_dims(A, 4).dims
    minimal = minimal_cohomology_dims(to_engine(A), 4, primes=(p,))[p]
    assert minimal == cs, f"CS vs minimal A^e disagree on A_3 over GF({p}): {cs} vs {minimal}"
    assert cs == [1, 3, 0, 2, 0]


def test_tri_string_A3_revival_cc():
    """A_3 = [1, 3, 0, 2, 0] over CC.  Char-independence of the revival is anchored by
    the GF(p) minimal-engine cross-check above (Thm 4.3 is a combinatorial count)."""
    assert cs_cohomology_dims(_tri_string(3, CC), 4).dims == [1, 3, 0, 2, 0]


# --- Plan 33: the revival at scale -- A_5 (the deepest revival reached here) --
@pytest.mark.parametrize("fld", [GF(2), GF(3), CC], ids=["GF2", "GF3", "CC"])
def test_tri_string_A5_revival_scale(fld):
    """Plan 33 SCALE of the A_3 revival: A_5 = A_{2m+1} (m = 2, dim 36) has
    HH^* = [1, 5, 0, 0, 0, 2, 0] -- HH^{2,3,4} vanish then HH^5 = 2 REVIVES (the
    A_{2m+1} family has HH^{2m+1} = 2; this is the deepest revival reached in the
    battery).  Pinned over GF(2), GF(3) AND CC: char-independence is Redondo-Roman
    Thm 4.3 (the dim is a combinatorial count), the very theorem that lets the A_3
    GF(p) minimal-A^e anchor above carry to every characteristic
    (``redondo_roman_2014``, Ex 3 / Thm 4.3).

    DEVIATION from the A_3 pattern: no in-test minimal-engine cross-check.  A_5 has
    dim 36, where the minimal A^e engine's iterated syzygies exceed the deep budget
    (minutes vs CS's ~0.15 s to degree 6); the revival MECHANISM is already anchored
    by the dim-16 A_3 minimal cross-check, and Thm 4.3's char-independence (verified
    here across three characteristics) carries it two degrees deeper."""
    assert cs_cohomology_dims(_tri_string(5, fld), 6).dims == [1, 5, 0, 0, 0, 2, 0]


# =========================================================================== #
# 1c.  Redondo-Roman (2014) Thm 5.2 -- cup-triviality in positive degrees       #
# =========================================================================== #
def test_cup_triviality_triangular_string_rr2014():
    """RR 2014 Thm 5.2: on ANY triangular string algebra the Hochschild ring is
    TRIVIAL in positive degrees -- HH^n cup HH^m = 0 for all n,m>0.  Decisive
    instance: A_2 (NOT hereditary -- a genuine string algebra of infinite global
    dimension, so HH^2 = 0 is a real computation) has HH^* = [1,2,0,0], i.e. a
    2-DIMENSIONAL HH^1 -- so we cup GENUINELY NONZERO classes and every product must
    land in HH^2 (which is 0-dimensional, making Thm 5.2's "=0" the specialization
    of the theorem to this member).  Both cup routes are exercised: the Plan-20
    native CS cup and the transported (in-window) cup.

    Citation: Redondo-Roman (2014), arXiv:1301.0516, Thm 5.2.  Verification status:
    fetched-source (Thm 5.2 verbatim)."""
    A = _tri_string(2, GF(3))
    assert cs_cohomology_dims(A, 3).dims == [1, 2, 0, 0]     # HH^1 = k^2, HH^2 = 0

    comp = Comparison(A)
    comp._ensure(3)                                          # native cup needs S up to 3
    assert comp.window >= 2, f"expected an in-window degree-2 cup, got window {comp.window}"
    reps1 = comp.cs_cohomology_basis(1)
    assert len(reps1) == 2, "A_2 must have a 2-dimensional HH^1 (real classes to cup)"
    zero2 = [0] * len(comp._res._basis(2, "coh"))

    for engine in ("native", "transport"):
        for i in range(2):
            for j in range(2):
                u, v = comp.hh_class_cs(1, i), comp.hh_class_cs(1, j)
                cup = comp.cup_of_cs_classes(u, v, engine=engine)
                assert comp.same_cohomology_class(cup, zero2, degree=2), \
                    f"HH^1 cup HH^1 is not a coboundary ({engine}, reps {i},{j}) -- " \
                    f"contradicts RR-2014 Thm 5.2"


@pytest.mark.skip(reason="RR-2018 (arXiv:1504.02495) Ex 3.1 cup-NONvanishing (Thms "
                         "4.8/4.9): the paper presents dim HH^n via combinatorial "
                         "sets, not assembled integer vectors, so the exact nonzero "
                         "products are convention-risky without a clean bar anchor. "
                         "Deferred per Plan 29 honest-scope (NOT pinned).")
def test_cup_nonvanishing_rr2018_deferred():
    """Placeholder documenting the deferred RR-2018 cup-nonvanishing oracle.  Would
    pin a NONZERO HH cup on the quadratic string algebra k Q / <a1 a2, a2 a1, b1 a2>
    (vertices 1,2; arrows a1,b1:1->2, a2:2->1; dim 6, infinite gl.dim) once an exact
    nonzero product is anchored through the bar/transported route."""
    pass


# =========================================================================== #
# 3.  Multi-vertex ZOO depth to degree 24 (Plan 33 -- CS is ~free this deep)     #
# =========================================================================== #
# The CS engine reaches degree 24 in milliseconds on the standing multi-vertex zoo
# records (Plan 18), far past the bar oracle's degree-5 ceiling on such algebras.
from quiverlab.families.zoo import load_catalog, build_from_record  # noqa: E402


def _zoo(name, field):
    rec = next(r for r in load_catalog() if r.get("name") == name)
    return build_from_record(rec, field=field)


@pytest.mark.parametrize("name,nverts", [("comm_square", 4), ("line_abc_cde", 6)])
def test_zoo_acyclic_depth_cs_deg24(name, nverts):
    """The ACYCLIC multi-vertex zoo records (Plan 18): comm_square (dim 9) and
    line_abc_cde (dim 16) have no oriented cycles, so HH_n = 0 for n >= 1 to any
    depth (Cibils' acyclic vanishing, ``cibils_acyclic``) and HH^* = [1, 0, ...]
    (connected).  Pinned to DEGREE 24 via CS (~0.01 s) -- a depth showcase where the
    bar oracle is hopeless.  HH_0 = #vertices (the separable degree-0 part)."""
    A = _zoo(name, GF(32003))
    assert cs_cohomology_dims(A, 24).dims == [1] + [0] * 24
    assert cs_homology_dims(A, 24).dims == [nverts] + [0] * 24


@pytest.mark.oracle_crossengine
def test_zoo_cyclic_nakayama_depth_cs_vs_bardzell_deg24():
    """The CYCLIC zoo record cn_3_2 = kZ_3/rad^2 (dim 6, self-injective, radical
    square zero): its HH_* is PERIODIC (period 6) and nonzero infinitely often
    (complexity 1) -- a Cibils rad^2 = 0 computation (``cibils_radsq``).  Cross-engine
    at DEPTH: the CS engine (general reduction-system route) and the Bardzell engine
    (monomial minimal bimodule resolution) agree degreewise to degree 24 --
    HH_* = [3, 0, 1, 1, 0, 0, ...] repeating (``bardzell``; ``chouhy_solotar``).  The
    bar oracle cannot reach this depth on a self-injective algebra."""
    A = _zoo("cn_3_2", GF(32003))
    cs = cs_homology_dims(A, 24).dims
    from quiverlab.engine.hh_engine import hochschild_homology_dims
    from quiverlab.engine.resolutions_bardzell import BardzellResolution, MonomialPresentation
    from quiverlab.engine.coxeter2 import cyclic_nakayama
    Ae, _ = cyclic_nakayama(3, 2)
    res = BardzellResolution(MonomialPresentation.cyclic_nakayama(3, 2))
    bard = hochschild_homology_dims(Ae, 24, resolution=res)[32003]
    assert cs == bard, f"CS vs Bardzell disagree on cn_3_2 to deg 24: {cs} vs {bard}"
    assert cs[:8] == [3, 0, 1, 1, 0, 0, 0, 0]      # period-6 revival pattern
