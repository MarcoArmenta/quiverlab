"""Certified algebraic ρ/M, Φ_n labelling, off-circle count -- self-certifying
(minpoly annihilates, interval brackets, count sound) + the record's exact
literature values."""
import sympy as sp
import pytest

from quiverlab.errors import QuiverlabError
from quiverlab.invariants.spectral import (
    certify_real_algebraic, cyclotomic_factorization, off_circle_root_count,
    spectral_radius, mahler_measure, _noncyclotomic_part, _real_roots_suffice)

t, x = sp.Symbol("t"), sp.Symbol("x")
LEHMER = t**10 + t**9 - t**7 - t**6 - t**5 - t**4 - t**3 + t + 1


@pytest.mark.oracle_selfcert
def test_certificate_is_self_consistent():
    for alpha in (spectral_radius(t**2 - 7*t + 1), spectral_radius(LEHMER)):
        c = certify_real_algebraic(alpha)
        m = sp.Poly(c["minpoly"], x)
        # certify's minpoly is EXACTLY sympy's minimal polynomial of alpha -- a
        # COEFFICIENT ROUND-TRIP, not an independent annihilation check. (A direct
        # m(alpha)==0 via Poly.eval would spuriously FAIL on the degree-10 Lehmer
        # CRootOf, whose powers do not auto-reduce; this compares against the same
        # sympy.minimal_polynomial the certificate is built from.) The INDEPENDENT
        # anchoring -- that this is the RIGHT polynomial -- is carried by the two
        # literature pins below (3-Kronecker [1,-7,1] and Lehmer's degree-10
        # polynomial), whose coefficients are fixed by the record, not by sympy.
        assert m.as_expr() == sp.minimal_polynomial(alpha, x)      # round-trips sympy's minpoly
        a, b = sp.Rational(c["interval"][0]), sp.Rational(c["interval"][1])
        assert (alpha - a).is_nonnegative and (b - alpha).is_nonnegative  # brackets
        ivs = m.intervals()
        assert (a, b) == ivs[c["root_index"]][0]                   # index consistent
        assert not any(isinstance(z, float) for z in c["minpoly"])  # never a float


@pytest.mark.oracle_literature
def test_three_kronecker_certified_value():
    c = certify_real_algebraic(spectral_radius(t**2 - 7*t + 1))
    assert c["minpoly"] == [1, -7, 1] and c["degree"] == 2
    assert (sp.Rational(c["interval"][0]), sp.Rational(c["interval"][1])) == (6, 7)
    assert sp.simplify(sp.sympify(c["value"]) - (7 + 3*sp.sqrt(5))/2) == 0


@pytest.mark.oracle_literature
def test_lehmer_minpoly_and_interval():
    c = certify_real_algebraic(spectral_radius(LEHMER))
    assert sp.Poly(c["minpoly"], x) == sp.Poly(LEHMER.subs(t, x), x)
    assert (sp.Rational(c["interval"][0]), sp.Rational(c["interval"][1])) == (1, 2)


@pytest.mark.oracle_literature
def test_cyclotomic_labelling():
    facs = cyclotomic_factorization(t**2 + t + 1)                  # Φ_3
    assert facs == [{"factor": "t**2 + t + 1", "latex": facs[0]["latex"],
                     "multiplicity": 1, "cyclotomic_index": 3}]
    d4 = cyclotomic_factorization((t + 1)**2 * (t**2 - t + 1))     # Φ_2^2 · Φ_6
    idx = sorted((f["cyclotomic_index"], f["multiplicity"]) for f in d4)
    assert idx == [(2, 2), (6, 1)]
    mixed = cyclotomic_factorization(LEHMER)                       # irreducible, none
    assert [f["cyclotomic_index"] for f in mixed] == [None]


@pytest.mark.oracle_selfcert
def test_product_mahler_exceeds_spectral_radius():
    """M != ρ when TWO roots lie outside the circle -- the product-Mahler branch the
    single-outside-root battery never exercised. χ = (t²−7t+1)(t²−14t+1) has outside
    roots (7+3√5)/2 and 7+4√3, so the Mahler measure M = their product (a degree-4
    algebraic number, minpoly [1,−98,243,−98,1], interval (95,96)) STRICTLY exceeds the
    spectral radius ρ = 7+4√3 (minpoly [1,−14,1], interval (13,14)). Both certify exactly
    and quickly (the critic measured ~0.01 s) -- the multi-outside-root M path is exact."""
    chi = (t**2 - 7*t + 1) * (t**2 - 14*t + 1)
    assert off_circle_root_count(chi) == 2
    rho, M = spectral_radius(chi), mahler_measure(chi)
    assert (M - rho).is_positive                                  # M strictly exceeds ρ
    cR = certify_real_algebraic(rho)
    assert cR["minpoly"] == [1, -14, 1] and cR["degree"] == 2
    assert (sp.Rational(cR["interval"][0]), sp.Rational(cR["interval"][1])) == (13, 14)
    cM = certify_real_algebraic(M)
    assert cM["minpoly"] == [1, -98, 243, -98, 1] and cM["degree"] == 4
    assert (sp.Rational(cM["interval"][0]), sp.Rational(cM["interval"][1])) == (95, 96)
    # the degree-4 certificate self-checks: minpoly == sympy's, Sturm-unique bracket
    m = sp.Poly(cM["minpoly"], x)
    assert m.as_expr() == sp.minimal_polynomial(M, x)
    assert (sp.Rational(cM["interval"][0]), sp.Rational(cM["interval"][1])) \
        == m.intervals()[cM["root_index"]][0]


@pytest.mark.oracle_selfcert
def test_rational_and_count():
    c = certify_real_algebraic(sp.Integer(1))
    assert c["is_rational"] and c["degree"] == 1 and c["minpoly"] == [1, -1]
    assert off_circle_root_count(t**3 - 1) == 0                    # cyclotomic
    assert off_circle_root_count(t**2 - 7*t + 1) == 1              # one Salem/Pisot root
    assert off_circle_root_count(LEHMER) == 1
    # The complex-dominant gate (§D3): on the dim-5 non-hereditary rad²=0 example the
    # dominant root is a |complex root| and minimal_polynomial HANGS (~121 s, measured)
    # -- it does NOT raise. The deterministic shipped predicate must decide
    # complex-dominance up front so the assembler refuses ρ/M without ever calling it,
    # while STILL reporting the outside-circle count (a conjugate pair).
    chi_cx = (t + 1) * (t**4 - 7*t**3 + 16*t**2 - 7*t + 1)
    assert _real_roots_suffice(_noncyclotomic_part(chi_cx)) is False
    assert off_circle_root_count(chi_cx) == 2
