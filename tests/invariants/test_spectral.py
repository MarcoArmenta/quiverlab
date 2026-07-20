"""Exact spectral radius / Mahler measure (spec section 5 component 8), SOUND for complex
off-circle roots. The bank floats survive only as test oracles here; src is float-free
(asserted by tests/test_no_floats.py)."""
import sympy as sp
from quiverlab import Quiver, CC
from quiverlab.invariants.spectral import spectral_radius, mahler_measure
from quiverlab.engine.coxeter_spectrum import is_cyclotomic_product

t = sp.symbols("t")
LEHMER = t**10 + t**9 - t**7 - t**6 - t**5 - t**4 - t**3 + t + 1


def _rad2_algebra(M, field=CC):
    """The radical-square-zero algebra whose arrow-count matrix is M (Cartan = I + M):
    M[i][j] parallel arrows i+1 -> j+1, all length-2 paths killed."""
    n = len(M)
    arrows = {}
    for i in range(n):
        for j in range(n):
            for c in range(M[i][j]):
                arrows[f"x{i+1}_{j+1}_{c}"] = (i + 1, j + 1)
    Q = Quiver(list(range(1, n + 1)), arrows)
    rels = [f"{a}*{b}" for a, (_sa, ta) in arrows.items()
            for b, (sb, _tb) in arrows.items() if ta == sb]
    return Q.algebra(relations=rels, field=field)


def test_returns_exact_types_and_none():
    assert spectral_radius(None) is None
    assert mahler_measure(None) is None
    rho = spectral_radius(t**2 - 3 * t + 1)
    assert not isinstance(rho, float)              # exact sympy object, never a float


def test_constant_and_zero_poly_are_guarded_not_bogus_one():
    # Fix 3: a constant / zero polynomial has degree < 1 -> no reciprocal-root
    # structure, so both entry points must return the guarded value (None, matching
    # the poly-is-None convention), NEVER the silent bogus 1 they used to yield.
    for p in (sp.Integer(5), sp.Integer(1), sp.Integer(0), sp.Integer(-7), sp.Integer(2) * t**0):
        assert spectral_radius(p) is None
        assert mahler_measure(p) is None
    # existing valid-degree cases are unchanged
    assert spectral_radius(t**3 - 1) == 1
    assert mahler_measure((t + 1)**2 * (t**2 - t + 1)) == sp.Integer(1)


def test_cyclotomic_short_circuit_is_exact_one():
    assert spectral_radius((t + 1)**3 * (t - 1)**2) == sp.Integer(1)
    assert mahler_measure((t + 1)**2 * (t**2 - t + 1)) == sp.Integer(1)
    assert spectral_radius(t**3 - 1) == 1
    assert is_cyclotomic_product(LEHMER) is False   # Lehmer is NOT cyclotomic


def test_quadratic_wild_radius_is_golden_squared():
    rho = spectral_radius(t**2 - 3 * t + 1)         # Branch A (real off-circle roots)
    assert sp.simplify(rho - (3 + sp.sqrt(5)) / 2) == 0      # exact
    assert abs(float(rho.evalf(30)) - (3 + 5**0.5) / 2) < 1e-9   # bank float oracle
    assert sp.simplify(mahler_measure(t**2 - 3 * t + 1) - rho) == 0


def test_lehmer_spectral_radius_is_lehmer_number():
    rho = spectral_radius(LEHMER)                   # Branch A (Salem: one real root > 1)
    assert abs(float(rho.evalf(30)) - 1.17628081825991) < 1e-9   # Lehmer's number
    assert abs(float(mahler_measure(LEHMER).evalf(30)) - 1.17628081825991) < 1e-9


def test_mahler_of_mixed_product():
    m = mahler_measure((t + 1)**4 * (t**2 - 3 * t + 1))   # cyclotomic (t+1)^4 -> 1
    assert abs(float(m.evalf(30)) - (3 + 5**0.5) / 2) < 1e-9


def test_complex_off_circle_roots_counterexample_polynomial():
    # SOUNDNESS: (t+1)(t^4-7t^3+16t^2-7t+1) has a COMPLEX conjugate pair of modulus
    # 3.54645... that real_roots-only would drop (returning a wrong 1). Branch B.
    poly = (t + 1) * (t**4 - 7 * t**3 + 16 * t**2 - 7 * t + 1)
    assert is_cyclotomic_product(poly) is False
    rho = spectral_radius(poly)
    assert abs(float(rho.evalf(30)) - 3.54645544468500) < 1e-9
    assert (rho - 3).is_positive                    # exactly > 3, not a silent 1
    m = mahler_measure(poly)
    assert abs(float(m.evalf(30)) - 12.5773462211358) < 1e-9


def test_complex_off_circle_roots_end_to_end_rad2_algebra():
    # SOUNDNESS end-to-end: a genuine dim-22 rad^2=0 algebra whose Coxeter polynomial
    # t^4+4t^3+14t^2+4t+1 has complex modulus-3.58474... roots. real_roots-only returned 1.
    M = [[0, 0, 2, 0], [1, 0, 1, 1], [3, 3, 0, 1], [2, 2, 2, 0]]
    A = _rad2_algebra(M)
    assert A.dim == 22
    p = A.coxeter_polynomial()
    assert sp.expand(p.as_expr() - (t**4 + 4 * t**3 + 14 * t**2 + 4 * t + 1)) == 0
    rho = spectral_radius(p.as_expr())
    assert rho != 1 and (rho - sp.Rational(35, 10)).is_positive    # NOT a silent 1
    assert abs(float(rho.evalf(30)) - 3.58474330285921) < 1e-9


def test_spectral_of_cartan_of_quiver_first_ever_coverage():
    # first-ever assertion coverage of star_quiver + cartan_of_quiver (ledger)
    from quiverlab.engine.coxeter_spectrum import star_quiver, cartan_of_quiver
    from quiverlab.engine.coxeter2 import coxeter_polynomial_from_cartan
    n, arrows = star_quiver([1, 1, 1])         # D_4 star, 4 vertices
    assert n == 4
    C, _ = cartan_of_quiver(n, arrows)
    poly, _Phi = coxeter_polynomial_from_cartan(C)
    assert spectral_radius(poly) == 1          # Dynkin D_4 is cyclotomic -> radius 1
