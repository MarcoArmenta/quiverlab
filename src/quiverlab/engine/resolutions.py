# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
"""Pluggable projective-resolution backends for the Hochschild (co)homology engine.

A `Resolution` supplies the six resolution-specific primitives the rest of hanlab
depends on; everything downstream (rank_mod_p, the rank identities, quotient_induced,
the induced-action descent, the scans) is resolution-agnostic. `BarResolution` is the
reference backend -- the normalized bar complex implemented in hh_engine / scan3 /
coxeter -- and the ground-truth oracle against which faster/smaller backends
(Bardzell, Chouhy-Solotar; see notes/B1_chouhy_solotar_plan.md) must be cross-checked.

Contract every backend must honour:
  * matrices are integer (int64), NEVER pre-reduced mod p (small primes carry signal);
  * differential_matrix returns shape (dim_{n-1}, dim_n); coboundary_matrix returns
    (dim_{n+1}, dim_n) -- so the rank identities in the consumers are unchanged;
  * term_basis / cochain_basis return ordered lists of hashable generators, and the
    differential matrices are expressed in exactly that ordering.

Imports of the bar-complex functions are done lazily inside the methods so that this
module can be imported in any order relative to hh_engine / scan3 / coxeter.
"""

import numpy as np


class Resolution:
    """Interface for a projective-resolution backend (see module docstring)."""

    def term_basis(self, alg, n):
        raise NotImplementedError

    def differential_matrix(self, alg, n, basis_n, index_nm1):
        raise NotImplementedError

    def cochain_basis(self, alg, n):
        raise NotImplementedError

    def coboundary_matrix(self, alg, n, basis_n, index_np1):
        raise NotImplementedError

    def sigma_chain_matrix(self, alg, n, S, p):
        raise NotImplementedError

    def sigma_cochain_matrix(self, alg, n, S, Sinv, p):
        raise NotImplementedError


class BarResolution(Resolution):
    """The normalized bar complex C_n = A (x) (A/k1)^{(x)n}: the reference backend."""

    def term_basis(self, alg, n):
        from quiverlab.engine.hh_engine import cn_basis
        return cn_basis(alg, n)

    def differential_matrix(self, alg, n, basis_n, index_nm1):
        from quiverlab.engine.hh_engine import differential_matrix
        return differential_matrix(alg, n, basis_n, index_nm1)

    def cochain_basis(self, alg, n):
        from quiverlab.engine.scan3 import cochain_basis
        return cochain_basis(alg, n)

    def coboundary_matrix(self, alg, n, basis_n, index_np1):
        from quiverlab.engine.scan3 import coboundary_matrix
        return coboundary_matrix(alg, n, basis_n, index_np1)

    def sigma_chain_matrix(self, alg, n, S, p):
        from quiverlab.engine.coxeter import sigma_chain_matrix
        return sigma_chain_matrix(alg, n, S, p)

    def sigma_cochain_matrix(self, alg, n, S, Sinv, p):
        from quiverlab.engine.coxeter import sigma_cochain_matrix
        return sigma_cochain_matrix(alg, n, S, Sinv, p)


class TruncatedPolynomialResolution(Resolution):
    """Periodic minimal bimodule resolution of A = k[x]/(x^a) -- the first small
    backend (B1 Increment 2), proving the depth unlock.

    The minimal A^e-resolution of k[x]/(x^a) is periodic with a single free generator
    per degree: P_n = A^e for all n, with differentials alternating multiplication by
    u = x(x)1 - 1(x)x and by v = sum_i x^i (x) x^{a-1-i}. After A (x)_{A^e} - the
    complex becomes A in every degree with induced maps
        d_n = 0                        (n odd; u acts as x*(-) - (-)*x = 0 on a
                                        commutative algebra),
        d_n = multiply by a*x^{a-1}    (n even; v acts as sum_i x^i(-)x^{a-1-i}
                                        = a*x^{a-1}*(-)).
    So dim P_n = a for all n (vs the bar complex's a*(a-1)^n), and HH_n is computable
    to arbitrary depth. Integer entries are kept un-reduced, so a prime p | a correctly
    collapses every map (HH_n = a for all n in characteristic p).

    Homology only (term_basis + differential_matrix); validated against the bar engine
    in tests/test_truncated_resolution.py. Intended for algebras built by
    truncated_polynomial(a); reads a = alg.m.
    """

    def term_basis(self, alg, n):
        return list(range(alg.m))                 # monomials 1, x, ..., x^{a-1}

    def differential_matrix(self, alg, n, basis_n, index_nm1):
        a = alg.m
        M = np.zeros((a, a), dtype=np.int64)       # shape (dim_{n-1}, dim_n) = (a, a)
        if n % 2 == 0:                             # d_n = multiply by a*x^{a-1}
            M[a - 1, 0] = a                        # x^j -> a*x^{a-1+j}: nonzero only j=0
        return M                                   # n odd: d_n = 0


def _default(resolution):
    """Return `resolution`, or a fresh BarResolution if it is None."""
    return BarResolution() if resolution is None else resolution
