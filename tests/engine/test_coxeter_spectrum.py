"""Guards for hanlab/coxeter_spectrum.py (Search VII).

Pins: the exact cyclotomicity test, the spectral radius / Mahler measure layer, the
Lehmer star T(2,3,7), the Q21 reconciliation (t^3-1 vs t^3+1 as the suspension sign),
the K0-blindness collapse on trivial extensions / local algebras, the paracyclic
column degeneration (monodromy induces the identity on the twisted homological
columns -- ledger Q3), its NEGATIVE control (the Nakayama action on the untwisted
column is genuinely nontrivial), and the duality reflection (cochain slot map on
column 1 = Theta on the Han column)."""

import numpy as np
import sympy as sp
import pytest

from quiverlab.invariants.spectral import spectral_radius, mahler_measure
from quiverlab.engine.coxeter_spectrum import (is_cyclotomic_product,
                              cartan_of_quiver, trivial_extension_cartan, star_quiver,
                              nakayama_charpoly_hh, column_degeneration,
                              dual_column_action, quantum_ci_algebra)
from quiverlab.engine.coxeter2 import (coxeter_polynomial_from_cartan, cyclic_nakayama,
                      charpoly_of_induced)
from quiverlab.engine.coxeter import nakayama_automorphism, induced_on_HH_homology, is_identity

t = sp.symbols('t')
P = 32003

LEHMER = t**10 + t**9 - t**7 - t**6 - t**5 - t**4 - t**3 + t + 1


# --------------------------------------------------------------------
# spectral layer
# --------------------------------------------------------------------
def test_cyclotomic_detection():
    assert is_cyclotomic_product(t**3 - 1) is True
    assert is_cyclotomic_product(t**3 + 1) is True
    assert is_cyclotomic_product((t + 1)**2 * (t**2 - t + 1)) is True
    assert is_cyclotomic_product(t**2 - 3*t + 1) is False
    assert is_cyclotomic_product((t + 1)**4 * (t**2 - 3*t + 1)) is False
    assert is_cyclotomic_product(LEHMER) is False
    assert is_cyclotomic_product(None) is None


def test_spectral_radius_and_mahler():
    # restored (was deferred to Plan 05); now against the exact reimplementation
    rho = spectral_radius(t**2 - 3 * t + 1)
    assert abs(float(rho.evalf(30)) - (3 + 5**0.5) / 2) < 1e-9
    assert sp.simplify(mahler_measure(t**2 - 3 * t + 1) - rho) == 0
    assert spectral_radius((t + 1)**3 * (t - 1)**2) == 1        # exact, cyclotomic
    m = mahler_measure((t + 1)**4 * (t**2 - 3 * t + 1))
    assert abs(float(m.evalf(30)) - (3 + 5**0.5) / 2) < 1e-9


def test_lehmer_star_237():
    # T(2,3,7) = star_quiver([1,2,6]); Coxeter poly is Lehmer's polynomial; spectral
    # radius is Lehmer's number 1.176280818259917...  (ledger: exact half restored)
    n, arrows = star_quiver([1, 2, 6])
    assert n == 10
    C, _alg = cartan_of_quiver(n, arrows)
    poly, _Phi = coxeter_polynomial_from_cartan(C)
    assert sp.expand(poly - LEHMER) == 0
    assert is_cyclotomic_product(poly) is False
    assert abs(float(spectral_radius(poly).evalf(30)) - 1.17628081825991) < 1e-9


# --------------------------------------------------------------------
# Q21: the two Coxeter polynomials of kZ_3/rad^2 reconciled
# --------------------------------------------------------------------
def test_q21_reconciliation_kZ3():
    alg, _ = cyclic_nakayama(3, 2)
    Pn = np.zeros((3, 3), dtype=np.int64)
    for i in range(3):
        Pn[(i + 1) % 3, i] = 1
    C = np.eye(3, dtype=np.int64) + Pn
    polyF, _Phi = coxeter_polynomial_from_cartan(C)
    polyT, d0 = nakayama_charpoly_hh(alg, P, 0)
    assert d0 == 3
    assert sp.expand(polyF - (t**3 + 1)) == 0
    assert sp.expand(polyT - (t**3 - 1)) == 0
    # the suspension-sign dictionary:  p_formal(t) = (-1)^3 * p_Theta(-t)
    assert sp.expand(polyF - (-1)**3 * polyT.subs(t, -t)) == 0
    # the twisted Cartan symmetry  C^T = P_nu^{-1} C
    assert np.array_equal(C.T, (np.linalg.matrix_power(Pn, 2) @ C))


# --------------------------------------------------------------------
# K0-blindness
# --------------------------------------------------------------------
def test_local_algebras_constant_coxeter():
    for a in (2, 3, 4):
        poly, Phi = coxeter_polynomial_from_cartan(np.array([[a]], dtype=np.int64))
        assert sp.expand(poly - (t + 1)) == 0
        assert Phi == -sp.eye(1)


def test_trivial_extension_collapse():
    # rep-finite input vs wild input: K0 spectra of T(A) are identically (t+1)^v
    C2, _ = cartan_of_quiver(2, [(0, 1)])                # kA_2 (rep-finite)
    nW, aW = star_quiver([1, 1, 1, 1, 1])                # 5-subspace (wild)
    CW, _ = cartan_of_quiver(nW, aW)
    pA, _ = coxeter_polynomial_from_cartan(C2)
    pW, _ = coxeter_polynomial_from_cartan(CW)
    assert spectral_radius(pA) == 1
    assert (spectral_radius(pW) - sp.Rational(26, 10)).is_positive       # > 2.6, exact
    for C, v in ((C2, 2), (CW, nW)):
        CT = trivial_extension_cartan(C)
        polyT, PhiT = coxeter_polynomial_from_cartan(CT)
        assert PhiT == -sp.eye(v)
        assert sp.expand(polyT - (t + 1)**v) == 0
    # Euclidean D~4: symmetrized Cartan is singular -> Coxeter undefined
    nE, aE = star_quiver([1, 1, 1, 1])
    CE, _ = cartan_of_quiver(nE, aE)
    polyE, PhiE = coxeter_polynomial_from_cartan(trivial_extension_cartan(CE))
    assert polyE is None and PhiE is None


def test_nakayama_singular_cartan_theta_still_defined():
    # kZ_4/rad^2: det C = 0 (formal Coxeter undefined) but Theta is fine
    alg, _ = cyclic_nakayama(4, 2)
    Pn = np.zeros((4, 4), dtype=np.int64)
    for i in range(4):
        Pn[(i + 1) % 4, i] = 1
    C = np.eye(4, dtype=np.int64) + Pn
    polyF, _ = coxeter_polynomial_from_cartan(C)
    assert polyF is None
    polyT, d0 = nakayama_charpoly_hh(alg, P, 0)
    assert d0 == 4
    assert sp.expand(polyT - (t**4 - 1)) == 0


# --------------------------------------------------------------------
# column degeneration (monodromy = Coxeter; ledger Q3) + negative control
# --------------------------------------------------------------------
def test_column_degeneration_kZ3():
    alg, _ = cyclic_nakayama(3, 2)
    rows = column_degeneration(alg, P, N=2, powers=(1,))
    assert all(r["hom_is_id"] for r in rows)
    assert any(r["dim_hom"] > 0 for r in rows)            # non-vacuous


def test_column_degeneration_quantum_ci():
    alg = quantum_ci_algebra(2)
    rows = column_degeneration(alg, P, N=2, powers=(1, 2))
    assert all(r["hom_is_id"] for r in rows)
    # the nu-twisted column of Lambda_q is genuinely nonzero in low degrees
    assert any(r["dim_hom"] >= 2 for r in rows if r["m"] == 1)


def test_degeneration_negative_control():
    # the SAME slot map on the UNTWISTED column is the Nakayama action -- NOT id.
    alg, _ = cyclic_nakayama(3, 2)
    S, _si = nakayama_automorphism(alg, P)
    M0, d0 = induced_on_HH_homology(alg, 0, S, P)
    assert d0 == 3 and not is_identity(M0, P)
    assert sp.expand(charpoly_of_induced(M0, P) - (t**3 - 1)) == 0


def test_duality_reflection():
    # cochain slot map on the nu-twisted cohomological column == Theta on HH_*
    alg, _ = cyclic_nakayama(3, 2)
    rows = dual_column_action(alg, P, N=2, mpow=1)
    assert all(r["matches_theta_on_HH"] for r in rows)
    alg = quantum_ci_algebra(2)
    rows = dual_column_action(alg, P, N=2, mpow=1)
    assert all(r["matches_theta_on_HH"] for r in rows)


def test_lambda_q_chi_eigenvalues_q_power_s():
    # chi_n(Lambda_q) has roots q^s, q^{-s} with s = 2*floor(n/2)+1  (CLMS pattern)
    q = 2
    alg = quantum_ci_algebra(q)
    S, _si = nakayama_automorphism(alg, P)
    for n in (1, 2, 3):
        Mn, dn = induced_on_HH_homology(alg, n, S, P)
        assert dn == 2
        chi = sp.Poly(charpoly_of_induced(Mn, P), t)
        s = 2 * (n // 2) + 1
        for root in (pow(q, s, P), pow(pow(q, P - 2, P), s, P)):
            assert int(chi.eval(root)) % P == 0
