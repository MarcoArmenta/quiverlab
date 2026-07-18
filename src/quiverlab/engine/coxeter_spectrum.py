# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
"""
coxeter_spectrum -- Search VII: K0-level Coxeter spectra vs operator-level spectra.

The mathematical background (paper: "Coxeter spectra at the operator level"; see also
notes/theorem_B_convention.md, which this module upgrades from note to instrument):

  * K0 level.  The classical Coxeter transformation Phi = -C^{-T}C of an elementary
    algebra is *blind* off the hereditary locus:
      - local algebras:      C = [dim A]  =>  Phi = [-1], char poly t+1, constant;
      - self-injective:      C^{-T}C = P_nu (Nakayama permutation), so Phi = -P_nu and
                             the char poly is a signed product of cyclotomics -- spectral
                             radius 1 and Mahler measure 1 *regardless* of representation
                             type or HH growth;
      - symmetric:           C = C^T  =>  Phi = -I, char poly (t+1)^v, constant.
    On the hereditary locus the spectrum is the classical tame/wild dictionary
    (rho = 1 iff the underlying graph is Dynkin/Euclidean); the tables below exhibit
    both regimes side by side.

  * Operator level.  On the self-injective locus the Serre-twist (higher Coxeter)
    action sigma_bullet is the Nakayama action on HH_bullet; its char polys chi_n are
    the higher Coxeter polynomials.  The *paracyclic monodromy* T = sigma^{(n+1)} of
    the m-th Serre-power column acts on the nu^m-twisted complex, and the para-mixed
    identity  b_sigma B + B b_sigma = id - T  forces T to induce the *identity* on the
    homology of each twisted column ("diagonal degeneration").  `column_degeneration`
    verifies this exactly, in F_p, for both homology and cohomology of the twisted
    columns -- the computational check of the monodromy = Coxeter identification
    (ledger Q3) and of the Q21 reconciliation.

Everything is exact: sympy over Q for the spectral layer, int64 ranks mod p for the
homological layer.  Driver:  python -m hanlab.coxeter_spectrum report

spectral_radius/mahler_measure: not ported; exact reimplementations arrive with the
invariants plan (Plan 05).
"""

import itertools
import json
import sys

import numpy as np
import sympy as sp

from quiverlab.engine.hh_engine import Algebra, hochschild_homology_dims

PRIME = 32003  # default large-prime char-0 proxy (mirrors hanlab.PRIME)
from quiverlab.engine.coxeter import (nakayama_automorphism, induced_on_HH_homology,
                     quotient_induced, sigma_chain_matrix, sigma_cochain_matrix,
                     nullspace_mod_p, colspace_basis_mod_p, is_identity,
                     inverse_mod_p)
from quiverlab.engine.coxeter2 import (charpoly_of_induced, coxeter_polynomial_from_cartan,
                      cartan_from_raw, quiver_path_algebra, dynkin_quiver,
                      cyclic_nakayama)
from quiverlab.engine.bimodule import (twisted_bimodule, cn_basis_coeff, differential_matrix_coeff,
                      cochain_basis_coeff, coboundary_matrix_coeff)

_t = sp.symbols('t')


# ====================================================================
# spectral layer on a (Coxeter) polynomial -- exact where it matters
# ====================================================================
def is_cyclotomic_product(poly):
    """True iff `poly` (sympy expr or Poly in t, integer coefficients) is, up to sign,
    a product of cyclotomic polynomials -- equivalently, all roots are roots of unity
    (Kronecker).  Exact: factor over Q and compare each irreducible factor against the
    finitely many cyclotomic polynomials of its degree (phi(n) = d  =>  n <= 2 d^2)."""
    if poly is None:
        return None
    P = sp.Poly(sp.expand(poly), _t)
    if P.degree() == 0:
        return bool(abs(P.LC()) == 1)
    const, factors = sp.factor_list(P.as_expr(), _t)
    if abs(const) != 1:
        return False
    for f, _mult in factors:
        fp = sp.Poly(f, _t)
        d = fp.degree()
        if d == 0:
            continue
        if fp.LC() < 0:
            fp = sp.Poly(-fp.as_expr(), _t)
        if fp == sp.Poly(_t, _t):          # root 0 is not a root of unity
            return False
        found = False
        for nn in range(1, 2 * d * d + 3):
            if sp.totient(nn) == d and fp == sp.Poly(sp.cyclotomic_poly(nn, _t), _t):
                found = True
                break
        if not found:
            return False
    return True


# _roots_abs, spectral_radius, mahler_measure: not ported; the float root-finding
# (mpmath nroots) layer and its callers move to the invariants plan (Plan 05), which
# supplies exact reimplementations.


# ====================================================================
# Cartan helpers (path algebras, trivial extensions)
# ====================================================================
def cartan_of_quiver(n, arrows, name="kQ"):
    """(Cartan matrix, Algebra) of the path algebra of an acyclic tree quiver."""
    alg, paths = quiver_path_algebra(n, arrows, name=name)
    C = cartan_from_raw(n, None, [(pth[0], pth[-1]) for pth in paths])
    return C, alg


def trivial_extension_cartan(C):
    """Cartan matrix of the trivial extension T(A) = A x DA:  C_T = C + C^T
    (dim e_y DA e_x = dim e_x A e_y).  Always symmetric, so the formal Coxeter
    matrix of T(A) is -I whenever det(C + C^T) != 0 -- the K0-blindness collapse."""
    C = np.asarray(C, dtype=np.int64)
    return C + C.T


def star_quiver(arm_lengths):
    """(n, arrows) for the star tree with the given arm lengths (edges per arm),
    all arms oriented toward the center (vertex 0).  T(2,3,7) = star_quiver([1,2,6])
    carries Lehmer's polynomial as its Coxeter polynomial."""
    n = 1 + sum(arm_lengths)
    # each arm is a chain v_L -> ... -> v_1 -> center(0)
    arrows = []
    v = 1
    for L in arm_lengths:
        chain = list(range(v, v + L))
        for i, w in enumerate(chain):
            arrows.append((w, 0 if i == 0 else chain[i - 1]))
        v += L
    return n, arrows


# ====================================================================
# Theta / Nakayama layer (operator level, K0 shadow)
# ====================================================================
def nakayama_charpoly_hh(alg, p=PRIME, n=0):
    """(char poly of the Nakayama action Theta on HH_n, dim HH_n).  For n = 0 and a
    self-injective Nakayama algebra this is det(t I - P_nu) -- the *unsigned*
    normalization of the Coxeter polynomial (Q21: t^3 - 1, not t^3 + 1)."""
    S, _Sinv = nakayama_automorphism(alg, p)
    M, d = induced_on_HH_homology(alg, n, S, p)
    return charpoly_of_induced(M, p), d


# ====================================================================
# column degeneration: the paracyclic monodromy on the m-th Serre power
# ====================================================================
def _mat_power_mod(S, k, p):
    out = np.eye(S.shape[0], dtype=np.int64)
    for _ in range(k):
        out = (out @ S) % p
    return out


def induced_on_twisted_homology(alg, Sm, p, n):
    """Induced action of the slot map (apply nu^m to the coefficient and to every
    reduced tensor slot -- the paracyclic monodromy T of column m, up to inverse)
    on HH_n(A, {}_1 A_{nu^m}).  Returns (matrix, dim)."""
    M = twisted_bimodule(alg, psi=Sm)
    basis_n = cn_basis_coeff(alg, M, n)
    idx_n = {g: i for i, g in enumerate(basis_n)}
    if n >= 1:
        basis_nm1 = cn_basis_coeff(alg, M, n - 1)
        idx_nm1 = {g: i for i, g in enumerate(basis_nm1)}
        bn = differential_matrix_coeff(alg, M, n, basis_n, idx_nm1) % p
        ker = nullspace_mod_p(bn, p)
    else:
        ker = np.eye(len(basis_n), dtype=np.int64)
    basis_np1 = cn_basis_coeff(alg, M, n + 1)
    bnp1 = differential_matrix_coeff(alg, M, n + 1, basis_np1, idx_n) % p
    im = colspace_basis_mod_p(bnp1, p)
    Sigma, _, _ = sigma_chain_matrix(alg, n, Sm, p)
    return quotient_induced(Sigma, im, ker, p)


def induced_on_twisted_cohomology(alg, Sm, Sminv, p, n):
    """Cohomological twin: the slot map on HH^n(A, {}_1 A_{nu^m})."""
    M = twisted_bimodule(alg, psi=Sm)
    basis_n = cochain_basis_coeff(alg, M, n)
    idx_n = {g: i for i, g in enumerate(basis_n)}
    basis_np1 = cochain_basis_coeff(alg, M, n + 1)
    idx_np1 = {g: i for i, g in enumerate(basis_np1)}
    dn = coboundary_matrix_coeff(alg, M, n, basis_n, idx_np1) % p
    if n >= 1:
        basis_nm1 = cochain_basis_coeff(alg, M, n - 1)
        idx_nm1 = {g: i for i, g in enumerate(basis_nm1)}
        dnm1 = coboundary_matrix_coeff(alg, M, n - 1, basis_nm1, idx_n) % p
        im = colspace_basis_mod_p(dnm1, p)
    else:
        im = np.zeros((len(basis_n), 0), dtype=np.int64)
    ker = nullspace_mod_p(dn, p)
    Sigma, _, _ = sigma_cochain_matrix(alg, n, Sm, Sminv, p)
    return quotient_induced(Sigma, im, ker, p)


def column_degeneration(alg, p=PRIME, N=3, powers=(1,)):
    """Verify, exactly over F_p, that the paracyclic monodromy slot map induces the
    IDENTITY on the *homological* nu^m-twisted column HH_n(A, {}_1 A_{nu^m}), for
    n <= N and each m in `powers` -- the para-mixed degeneration
    b_sigma B + B b_sigma = id - T, the computational face of the monodromy = Coxeter
    identification (ledger Q3).

    The degeneration is a statement about the homological columns: the *cochain* slot
    map on HH^n(A, {}_1 A_{nu^m}) is NOT a monodromy -- by the duality reflection
    D({}_1 A_{sigma}) = {}_1 A_{sigma^{-1} nu} it computes the Serre-twist action on
    the dual of the (1-m)-th homological column (see `dual_column_action`).

    Returns a list of dicts {m, n, dim_hom, hom_is_id}.  `quotient_induced` raises if
    the slot map is not a chain map, so a sign/convention error cannot silently pass."""
    S, _ = nakayama_automorphism(alg, p)
    out = []
    for mpow in powers:
        Sm = _mat_power_mod(S, mpow, p)
        for n in range(N + 1):
            Mh, dh = induced_on_twisted_homology(alg, Sm, p, n)
            out.append(dict(
                m=mpow, n=n,
                dim_hom=int(dh), hom_is_id=bool(is_identity(Mh, p)) if dh else True))
    return out


def dual_column_action(alg, p=PRIME, N=3, mpow=1):
    """The cochain slot map on the nu^m-twisted cohomological column
    HH^n(A, {}_1 A_{nu^m}).  For m = 1 this column is the dual of the *untwisted* Han
    column (D({}_1 A_nu) = A), and the slot map computes the Nakayama/Serre-twist
    action Theta there: its char polys must equal the Theta char polys chi_n on
    HH_n(A) degree by degree -- the duality-reflection cross-check.

    Returns rows {n, dim, charpoly} plus, when m = 1, {matches_theta_on_HH}."""
    S, _ = nakayama_automorphism(alg, p)
    Sm = _mat_power_mod(S, mpow, p)
    Sminv = inverse_mod_p(Sm, p)
    rows = []
    for n in range(N + 1):
        Mc, dc = induced_on_twisted_cohomology(alg, Sm, Sminv, p, n)
        row = dict(n=n, dim=int(dc), charpoly=str(charpoly_of_induced(Mc, p)))
        if mpow == 1:
            Mh, dh = induced_on_HH_homology(alg, n, S, p)
            row["matches_theta_on_HH"] = (int(dh) == int(dc) and
                str(charpoly_of_induced(Mh, p)) == row["charpoly"])
        rows.append(row)
    return rows


# ====================================================================
# quantum complete intersection builder (operator-level test fixture)
# ====================================================================
# The report driver (report / _hereditary_row) and its __main__ are not ported: they
# depend on the deleted spectral_radius / mahler_measure numeric layer (Plan 05).
def quantum_ci_algebra(q):
    """Lambda_q = k<x,y>/(x^2, y^2, yx - q xy): the quantum complete intersection, in
    the repo convention of `scan3.quantum_ci` (nu = diag(1, q^{-1}, q, 1), not inner
    for q != 1).  Local (single vertex), self-injective."""
    from quiverlab.engine.scan3 import quantum_ci
    return quantum_ci(q)
