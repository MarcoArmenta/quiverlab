"""Plan 33 -- nontrivial family examples AT SCALE (batteries).

Companion to the Plan-29 ``test_battery_literature_p29.py`` engine/CS batteries; this
file collects the FAMILY-constructed nontrivial examples the plan catalog lists
(``docs/plans/2026-07-26-plan-33-nontrivial-examples.md``): preprojective algebras
of Dynkin type, exterior algebras, commutative complete-intersection tensor
products, the Boolean-lattice incidence algebra, trivial extensions at scale, and
the generalized quantum complete intersection surface.

Oracle-class markers are PER TEST (Plan 32), because the honest-scope labels in the
plan are BINDING and split classes within the file:

* ``oracle_literature`` -- a paper-pinned dim / closed form / theorem identity
  (preprojective dims C(n+2,3), self-injectivity + Loewy = h-1, exterior Koszulity,
  the symmetric ``HH^n = HH_n`` identity, incidence acyclicity, the Cartan identity
  ``C_{T(A)} = C_A + C_A^T``, ``HH^1(T(A)) != 0``, the Bergh-Erdmann QCI values).
* ``oracle_crossengine`` -- two INDEPENDENT engines agree LIVE (CS vs bar, CS vs
  minimal A^e, the presented-vs-TensorProduct Kunneth build).  The preprojective
  and exterior HH VALUES are cross-engine ONLY -- NO literature pin (binding
  honest-scope: no published HH^*/HH_* table was consulted for those two families).
* ``oracle_selfcert`` -- an internal certificate (``dim T(A) = 2 dim A``).

``tests/families/`` auto-assigns to the deep CI bucket (``tests/conftest.py``).
Every value below was recomputed live (Plan 33) before pinning.
"""
import pytest

from quiverlab import Quiver, CC, GF
from quiverlab.families import (
    PreprojectiveAlgebra, ExteriorAlgebra, IncidenceAlgebra, TensorProduct,
    TrivialExtension, QuantumCI, truncated_polynomial, linear_path_algebra,
    dynkin_quiver,
)
from quiverlab.resolutions_cs.homology import cs_cohomology_dims, cs_homology_dims
from quiverlab.modules.koszul import g_quadratic_certificate
from quiverlab.engine.adapter import to_engine
from quiverlab.engine.resolutions_minimal import minimal_cohomology_dims

pytest.importorskip("quiverlab.groebner")

F = GF(32003)              # characteristic-0 proxy (fast int64 rank engine)


# =========================================================================== #
# 1.  Preprojective algebras of Dynkin type  Pi(A4/A5/D4/D5)  (Plan 33)         #
# =========================================================================== #
# Citation: the preprojective algebra Pi(Q) of a Dynkin quiver (``preprojective``,
# ``assem_book``; Erdmann-Snashall for self-injectivity + Loewy length).  Pi(A_n) is
# finite-dimensional of dim C(n+2, 3); Pi(D_n) of dim n(n-1)(2n-1)/3; it is
# self-injective with Loewy length = h - 1 (Coxeter number h).  Coxeter numbers:
# A4 -> 5, A5 -> 6, D4 -> 6, D5 -> 8, so Loewy = 4, 5, 5, 7.  HH values are cross-
# engine ONLY (binding honest-scope: no published table consulted).
_PREPROJ = [("A4", 20, 4), ("A5", 35, 5), ("D4", 28, 5), ("D5", 60, 7)]


@pytest.mark.oracle_literature
@pytest.mark.parametrize("typ,dim,loewy", _PREPROJ)
def test_preprojective_structural(typ, dim, loewy):
    """Pi(A4/A5/D4/D5): builder DIM (20/35/28/60), SELF-INJECTIVE, and Loewy length
    = h - 1 (4/5/5/7 for Coxeter numbers 5/6/6/8) (``preprojective``, ``assem_book``;
    Erdmann-Snashall).  Built with the Plan-33 auto degree bound (no kwarg)."""
    A = PreprojectiveAlgebra(typ, field=F)
    assert A.dim == dim
    assert A.is_selfinjective() is True
    assert A.loewy_length() == loewy


@pytest.mark.oracle_literature
@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_preprojective_dim_closed_form_An(n):
    """dim Pi(A_n) = C(n+2, 3): 4, 10, 20, 35 for n = 2..5 (the closed form verified
    on these members; ``preprojective``, ``assem_book``)."""
    from math import comb
    assert PreprojectiveAlgebra(f"A{n}", field=F).dim == comb(n + 2, 3)


@pytest.mark.oracle_literature
def test_preprojective_auto_bound_matches_explicit():
    """Plan-33 D4 (src): PreprojectiveAlgebra builds Pi(A5) and Pi(D5) with NO
    degree_bound kwarg (auto per-type table) -- the auto build is byte-identical to
    the explicit-bound build (same dim + Cartan) and reproduces the closed-form dims
    35 / 60 (``preprojective``, ``assem_book``)."""
    for typ, dim, bound in [("A5", 35, 12), ("D5", 60, 16)]:
        auto = PreprojectiveAlgebra(typ, field=F)
        explicit = PreprojectiveAlgebra(typ, field=F, degree_bound=bound)
        assert auto.dim == explicit.dim == dim
        assert auto.cartan_matrix() == explicit.cartan_matrix()


@pytest.mark.oracle_crossengine
def test_preprojective_hh_crossengine_pi_a4():
    """Pi(A4) (dim 20) HH^* cross-engine, shallow, NO literature pin (binding honest-
    scope): the CS engine agrees with the fast bar-rank engine at degree 1 ([2, 2])
    and with the INDEPENDENT minimal A^e engine (Plan 16) at degree 2 ([2, 2, 2]) --
    three disjoint resolutions of the same bimodule.  Only Pi(A4) is cross-checkable
    within budget: the fast engine's bar basis blows up past deg 1 at this dimension,
    the minimal A^e engine costs ~24 s at deg 2 (and hangs beyond dim ~30, so A5/D5
    are out), and CS's own degree bound blocks the non-quadratic-tip completion of
    A5/D4/D5."""
    A = PreprojectiveAlgebra("A4", field=F)
    cs2 = cs_cohomology_dims(A, 2).dims
    assert cs2 == [2, 2, 2]
    fast1 = A.hochschild_cohomology(1, engine="fast").dims
    assert fast1 == cs2[:2] == [2, 2]                      # CS vs fast bar-rank, deg 1
    minimal2 = minimal_cohomology_dims(to_engine(A), 2, primes=(32003,))[32003]
    assert minimal2 == cs2 == [2, 2, 2]                    # CS vs minimal A^e, deg 2


# =========================================================================== #
# 2.  Exterior algebras  Lambda(k^3), Lambda(k^4)  (Plan 33)                     #
# =========================================================================== #
# Citation: the exterior algebra Lambda(k^n) is Koszul (``priddy`` -- G-quadratic
# PBW certificate; ``froberg_koszul``) and self-injective with Loewy length n + 1.
# HH values are cross-engine ONLY (no published table consulted).
@pytest.mark.oracle_literature
@pytest.mark.parametrize("n,loewy", [(3, 4), (4, 5)])
def test_exterior_koszul_selfinjective_loewy(n, loewy):
    """Lambda(k^n) (dim 2^n = 8, 16): KOSZUL via the G-quadratic (Priddy PBW)
    certificate (``priddy``; NOT ext_algebra().koszul, which is a >120 s trap at this
    scale), SELF-INJECTIVE, Loewy length n + 1 (4, 5) (``froberg_koszul``,
    ``assem_book``).  Char-qualified over GF(32003)."""
    A = ExteriorAlgebra(n, field=F)
    assert A.dim == 2 ** n
    assert g_quadratic_certificate(A) is True
    assert A.is_selfinjective() is True
    assert A.loewy_length() == loewy


@pytest.mark.oracle_crossengine
@pytest.mark.parametrize("n,deg", [(3, 2), (4, 1)])
def test_exterior_hh_cs_vs_bar(n, deg):
    """Lambda(k^n) HH cross-engine, CS =(deg)= bar, NO literature pin (binding honest-
    scope -- no published HH table).  In the bar window (deg 2 for dim-8 Lambda(k^3);
    only deg 1 for dim-16 Lambda(k^4), where the bar basis blows up at deg 2) the CS
    engine and the exponential bar oracle agree degreewise, HH^* and HH_*
    (``priddy``/``chouhy_solotar`` engines; ``bar`` oracle).  Char-qualified over
    GF(32003).  CS reaches deeper alone (Lambda(k^3) HH^* = [5,12,24,40] to deg 3),
    but only the in-window bar agreement is a second-engine oracle."""
    A = ExteriorAlgebra(n, field=F)
    cs_co = cs_cohomology_dims(A, deg).dims
    cs_ho = cs_homology_dims(A, deg).dims
    bar_co = A.hochschild_cohomology(deg, engine="bar").dims
    bar_ho = A.hochschild_homology(deg, engine="bar").dims
    assert cs_co == bar_co, f"Lambda(k^{n}) coh CS != bar at deg {deg}: {cs_co} vs {bar_co}"
    assert cs_ho == bar_ho, f"Lambda(k^{n}) hom CS != bar at deg {deg}: {cs_ho} vs {bar_ho}"


# =========================================================================== #
# 3.  Commutative CI tensor products  (Kunneth + symmetric HH^n = HH_n)          #
# =========================================================================== #
# Citation: for a symmetric algebra HH^n = HH_n (the trace-form duality), and HH is
# multiplicative on tensor factors (Kunneth, ``tensor_product``).  k[x]/(x^3) (x)
# k[x]/(x^3) (dim 9) and k[x,y,z]/(x^2,y^2,z^2) = (k[x]/x^2)^{(x)3} (dim 8) are
# commutative Frobenius => symmetric.
def _comm_x3x3(field):
    """k<x,y>/(x^3, y^3, xy - yx) = k[x]/(x^3) (x) k[x]/(x^3), presented so CS runs
    past the bar oracle's dim-9 degree-3 blow-up."""
    return Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x^3", "y^3", "x*y - y*x"], field=field)


def _comm_triple(field):
    """k[x,y,z]/(x^2, y^2, z^2) = (k[x]/x^2)^{(x)3} (the COMMUTATIVE triple, dim 8)."""
    return Quiver([1], {"x": (1, 1), "y": (1, 1), "z": (1, 1)}).algebra(
        relations=["x^2", "y^2", "z^2", "x*y - y*x", "x*z - z*x", "y*z - z*y"], field=field)


@pytest.mark.oracle_literature
@pytest.mark.parametrize("builder,dim,expected", [
    (_comm_x3x3, 9, [9, 12, 16, 20, 24, 28, 32]),
    (_comm_triple, 8, [8, 12, 18, 25, 33, 42, 52]),
])
def test_tensor_ci_symmetric_hh(builder, dim, expected):
    """Commutative CI tensor products are SYMMETRIC, so HH^n = HH_n degreewise
    (``tensor_product`` Kunneth + the symmetric trace-form duality): to degree 6,
    k[x]/(x^3)^{(x)2} (dim 9) has HH^* = HH_* = [9,12,16,20,24,28,32], and
    k[x,y,z]/(x^2,y^2,z^2) (dim 8) has [8,12,18,25,33,42,52].  Computed via CS on the
    presented form (the bar oracle blows up at degree 3 for these dims)."""
    A = builder(F)
    assert A.dim == dim
    co = cs_cohomology_dims(A, 6).dims
    ho = cs_homology_dims(A, 6).dims
    assert co == ho == expected                            # symmetric: HH^n = HH_n


@pytest.mark.oracle_crossengine
def test_tensor_ci_kunneth_crosscheck():
    """Kunneth cross-engine: the explicit structure-constant TensorProduct build
    k[x]/(x^3) (x) k[x]/(x^3) (bar engine) agrees with the presented CS build at the
    bar window (degree 2): HH^* = HH_* = [9, 12, 16] both ways -- the Kunneth
    isomorphism realized as agreement between two INDEPENDENT constructions
    (structure-constant tensor vs quiver presentation) and two engines
    (``tensor_product``; ``bar``/``chouhy_solotar``)."""
    T = TensorProduct(truncated_polynomial(3, field=F), truncated_polynomial(3, field=F))
    pres = _comm_x3x3(F)
    assert T.hochschild_cohomology(2).dims == cs_cohomology_dims(pres, 2).dims == [9, 12, 16]
    assert T.hochschild_homology(2).dims == cs_homology_dims(pres, 2).dims == [9, 12, 16]


# =========================================================================== #
# 4.  Incidence algebra of the Boolean lattice B_3  (Plan 33)                    #
# =========================================================================== #
# Citation: HH^n of an incidence algebra = simplicial cohomology of the poset order
# complex (Gerstenhaber-Schack / Cibils, ``cibils_incidence``, ``redondo_incidence``).
# B_3 has a 0-hat and 1-hat, so its order complex is CONTRACTIBLE => HH^0 = 1,
# HH^{>=1} = 0 at any depth.
def _boolean_b3_covers():
    """Cover relations of the Boolean lattice 2^{a,b,c}: x < y iff y = x with one more
    bit set (subsets encoded as bitmasks 0..7)."""
    return [(x, x | (1 << b)) for x in range(8) for b in range(3) if not (x >> b) & 1]


@pytest.mark.oracle_literature
@pytest.mark.parametrize("field", [CC, F], ids=["CC", "GF32003"])
def test_incidence_boolean_b3(field):
    """Incidence algebra of the Boolean lattice B_3 (dim 27, 8 vertices, 12 Hasse
    covers): the order complex is CONTRACTIBLE (B_3 has 0-hat and 1-hat), so
    HH^* = [1, 0, 0, ...] to degree 10 (Cibils/Gerstenhaber-Schack incidence-vs-nerve
    identity, ``cibils_incidence``, ``redondo_incidence``).  Non-monomial admissible
    => CS; acyclic => cheap even to degree 10."""
    A = IncidenceAlgebra(_boolean_b3_covers(), field=field)
    assert A.dim == 27
    assert cs_cohomology_dims(A, 10).dims == [1] + [0] * 10


# =========================================================================== #
# 5.  Trivial extensions AT SCALE  T(kD4), T(kA5), T(kA6)  (Plan 33 / Plan 31)   #
# =========================================================================== #
# Citation: T(A) = A |x D(A) is symmetric for every f.d. A (``happel_trivial_extension``,
# ``skowronski_yamagata``); C_{T(A)} = C_A + C_A^T (``assem_book``); HH^1(T(A)) != 0
# always (Cibils-Marcos-Redondo-Solotar, ``cmrs_split``).  Extends the Plan-31
# presented battery (kA2/3/4) to kD4 and the larger kA5/kA6 (dim 42).
_TE_BASES = {
    "kD4": lambda: dynkin_quiver("D4", "linear").algebra(relations=[], field=F),
    "kA5": lambda: linear_path_algebra(5, field=F),
    "kA6": lambda: linear_path_algebra(6, field=F),
}


def _plus_transpose(C):
    n = len(C)
    return [[C[i][j] + C[j][i] for j in range(n)] for i in range(n)]


@pytest.mark.oracle_literature
@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("name", ["kD4", "kA5", "kA6"])
def test_trivial_extension_scale(name):
    """T(kD4) (dim 18), T(kA5) (dim 30), T(kA6) (dim 42): the ``dim T = 2 dim A``
    CERTIFICATE (``oracle_selfcert``); SYMMETRIC / weakly symmetric / Frobenius /
    self-injective (``happel_trivial_extension``, ``skowronski_yamagata``); the Cartan
    identity C_{T(A)} = C_A + C_A^T (``assem_book``); and HH^1(T(A)) != 0
    (``cmrs_split`` -- Z(A) is always a summand).  HH^1 via bar for the dim-18 T(kD4)
    ([5, 1]) and via CS for the dim-30/42 T(kA5)/T(kA6) ([6, 1] / [7, 1]) where the
    bar basis blows up at degree 1."""
    A = _TE_BASES[name]()
    T = TrivialExtension(A)
    assert T.quiver is not None                            # presented route (not the |x fallback)
    assert T.dim == 2 * A.dim                              # dim certificate (selfcert)
    assert T.is_symmetric() is True
    assert T.is_weakly_symmetric() is True
    assert T.is_frobenius() is True
    assert T.is_selfinjective() is True
    assert T.cartan_matrix() == _plus_transpose(A.cartan_matrix())
    if name == "kD4":
        hh1 = T.hochschild_cohomology(1, engine="bar").dims      # dim 18: bar feasible
    else:
        hh1 = cs_cohomology_dims(T, 1).dims                      # dim 30/42: CS
    assert hh1[1] >= 1, f"HH^1(T({name})) must be nonzero (cmrs_split): {hh1}"


# =========================================================================== #
# 6.  Generalized quantum CI surface  QuantumCI(q, a, b)  (Plan 33 D3)           #
# =========================================================================== #
# The Plan-33 src generalization QuantumCI(q, a, b) builds k<x,y>/(x^a, y^b,
# xy + q yx) (quiverlab's established xy + q yx convention).  Its HH matches the raw-
# Quiver Bergh-Erdmann pins in test_battery_literature_p29.py::
# test_bergh_erdmann_qci_scale_deg8 (values are convention-independent at generic q).
@pytest.mark.oracle_literature
@pytest.mark.parametrize("a,b", [(2, 4), (3, 4), (4, 4), (2, 5), (5, 5)])
def test_quantum_ci_family_surface_deg4(a, b):
    """The generalized ``QuantumCI(q=2, a, b)`` FAMILY surface (Plan 33 D3): dim = a*b
    (8/12/16/10/25) and HH^* = [2,2,1,0,0], HH_* = [a+b-1] + [a+b-2]*(...) to degree 4
    over CC -- the Bergh-Erdmann signature reproduced through the no-code family
    builder (``qci_hh_oracle``).  q = 2 is not a root of unity (CC only) and avoids
    the convention-divergent q in {1, -1, E(3)}."""
    A = QuantumCI(2, a, b, field=CC)
    assert A.dim == a * b
    assert cs_cohomology_dims(A, 4).dims == [2, 2, 1, 0, 0]
    assert cs_homology_dims(A, 4).dims == [a + b - 1] + [a + b - 2] * 4
