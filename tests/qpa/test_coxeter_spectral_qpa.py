"""QPA (GAP) crosscheck for the Plan-58 Coxeter spectral surface.

QPA 1.37 binds ``CoxeterPolynomial(A)`` on a quiver algebra and it agrees with
quiverlab's ``coxeter_polynomial`` = charpoly(-C^{-T}C) EXACTLY (variable rename
x_1 -> t; probe 2026-08-07: kA2 -> x_1^2 + x_1 + 1). That anchors the FOUNDATION of
the whole spectral surface cross-engine.

QPA has NO Mahler / spectral-radius / cyclotomic-type recognizer surface, so ρ, M,
the cyclotomic verdict, and the Lehmer-class note are covered by literature pins
(dlPena 2014/2013, Lehmer/E10), cross-engine (the shipped spectral_radius /
mahler_measure primitives), and self-cert (minpoly annihilates + Sturm interval +
ρ = M when one root is outside). The honest-scope guard below FAILS (never skips) if
QPA ever grows such a surface -- the standing signal to wire a real crosscheck.

qpa-marked: skips locally without GAP, mandatory under QUIVERLAB_REQUIRE_QPA=1.
"""
import sympy as sp
import pytest

from quiverlab import Quiver, linear_path_algebra
from quiverlab.families import dynkin_quiver
from quiverlab.fields import QQ
from quiverlab.qpa import session
from quiverlab.qpa.scripts import quiver_and_algebra_script

pytestmark = [pytest.mark.qpa,
              pytest.mark.skipif(session.should_skip_qpa(),
                                 reason="[qpa] backend not installed")]

t = sp.Symbol("t")


def _kron(m, field=QQ):
    return Quiver([1, 2], {f"a{i}": (1, 2) for i in range(m)}).algebra(field=field)


def _mine_low_to_high(A):
    """quiverlab's Coxeter polynomial coefficients, low->high degree (QPA order)."""
    poly = sp.Poly(A.coxeter_polynomial().as_expr(), t)
    return [int(c) for c in reversed(poly.all_coeffs())]


def _qpa_low_to_high(A):
    """QPA's CoxeterPolynomial(A) coefficients, low->high degree."""
    lg = session.libgap_handle()
    session.run(quiver_and_algebra_script(A))                 # binds global A
    coeffs = lg.eval("CoefficientsOfUnivariatePolynomial(CoxeterPolynomial(A))")
    return [int(c) for c in coeffs]


def test_coxeter_polynomial_matches_qpa():
    """kA2, kA4, 3-Kronecker, D4, E6 over QQ: CoxeterPolynomial agrees exactly."""
    algebras = {
        "kA2": linear_path_algebra(2, field=QQ),
        "kA4": linear_path_algebra(4, field=QQ),
        "3-Kronecker": _kron(3),
        "D4": dynkin_quiver("D4").algebra(relations=[], field=QQ),
        "E6": dynkin_quiver("E6").algebra(relations=[], field=QQ),
    }
    for name, A in algebras.items():
        assert _qpa_low_to_high(A) == _mine_low_to_high(A), f"{name} Coxeter mismatch"


# Names QPA might plausibly use for a Mahler / spectral-radius / Lehmer /
# cyclotomic-type ALGEBRA recognizer surface (Cyclotomic* in GAP are polynomial
# CONSTRUCTORS -- CyclotomicPolynomial / Cyclotomics -- not predicates on an algebra).
_SPECTRAL_NAMES = (
    "MahlerMeasure",
    "SpectralRadius",
    "LehmerPolynomial",
    "IsOfCyclotomicType",
    "CyclotomicTypeOfAlgebra",
    "SpectralRadiusOfCoxeterTransformation",
)


def test_qpa_has_no_spectral_surface():
    lg = session.libgap_handle()
    # (1) Scan the global name table FIRST (an IsBoundGlobal query REGISTERS the name
    #     into NamesGVars(), so scanning after would echo the probes back).
    gvar_names = [str(n) for n in lg.eval("NamesGVars()")]
    spectral_like = sorted(
        n for n in gvar_names
        if any(w in n.lower() for w in ("mahler", "lehmer"))
        or ("spectral" in n.lower() and "radius" in n.lower()))
    # (2) None of the named spectral entry points are bound (callable).
    bound = {name: bool(lg.eval(f'IsBoundGlobal("{name}")'))
             for name in _SPECTRAL_NAMES}
    present = [name for name, ok in bound.items() if ok]

    if not present and not spectral_like:
        pytest.skip(
            "QPA 1.37 exposes no Mahler / spectral-radius / cyclotomic-type recognizer "
            "surface (CoxeterPolynomial is present and IS crosschecked; Cyclotomic* are "
            "GAP polynomial constructors, not algebra predicates). Covering oracles for "
            "ρ/M/cyclotomic-verdict/Lehmer-class: literature pins (dlPena 2014/2013, "
            "Lehmer/E10) + cross-engine (spectral_radius/mahler_measure) + self-cert "
            "(minpoly + Sturm interval).")

    pytest.fail(
        "QPA now exposes a spectral surface -- wire a real crosscheck against the "
        f"Plan-58 rho/M/cyclotomic tables. Found: present={present}, "
        f"spectral_like_names={spectral_like}.")
