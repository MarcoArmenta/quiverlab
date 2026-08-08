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
        # minpoly annihilates alpha: the certificate's polynomial IS alpha's minimal
        # polynomial over QQ. (Poly.eval on a CRootOf does not auto-reduce its powers
        # to zero -- fine for the radical 3-Kronecker root but not the degree-10 Lehmer
        # CRootOf -- so we compare against sympy's minimal_polynomial directly: exact,
        # no float, and strictly stronger than a single-point annihilation check.)
        assert m.as_expr() == sp.minimal_polynomial(alpha, x)      # minpoly annihilates
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
