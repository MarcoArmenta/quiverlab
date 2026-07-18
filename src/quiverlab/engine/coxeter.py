# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
"""
Search IV -- the COXETER / NAKAYAMA route.

Marco Armenta's theory: the automorphism Theta of the Tamarkin-Tsygan calculus induced by the
Auslander-Reiten / Serre (= Nakayama, in the self-injective case) bimodule acts as
   * the IDENTITY on Hochschild cohomology HH^*        (Thm C; recovers Suarez-Alvarez),
   * the COXETER transformation on Hochschild homology HH_*  (Thm B; char poly = Coxeter poly).

This is exactly an "asymmetric ingredient": something homology carries that cohomology does not.
Proof-relevant consequence: if Theta has an eigenvalue != 1 on HH_n, then HH_n != 0, and this
certificate is INVISIBLE to cohomology (where Theta = id).  A nontrivial Coxeter eigenvalue that
persists across degrees would force non-vanishing in infinitely many degrees -- Han's conjecture --
by an argument no cohomological method can make.

This module computes the induced action of an algebra automorphism sigma on HH_n(A) and HH^n(A)
(as explicit matrices on the (co)homology spaces, over F_p) and compares to the identity.

We focus on LOCAL Frobenius algebras (one idempotent = the unit), where sigma preserves the radical
= reduced space cleanly:
   * symmetric (k[x]/x^a, k[x,y]/(x^2,y^2)):     Nakayama automorphism nu = id,
   * quantum complete intersection Lambda_q:      nu = diag(1, q^{-1}, q, 1) on {1,x,y,xy}, nontrivial.
"""

import json
import sys
import numpy as np
from quiverlab.engine.hh_engine import (Algebra, differential_matrix, cn_basis,
                       hochschild_homology_dims, truncated_polynomial, two_gen_local)
from quiverlab.engine.scan3 import (cochain_basis, coboundary_matrix, hochschild_cohomology_dims)


# ====================================================================
# exact linear algebra over F_p
# ====================================================================
def rref_mod_p(M, p):
    A = (np.asarray(M, dtype=np.int64) % p).copy()
    rows, cols = A.shape
    piv = []
    r = 0
    for c in range(cols):
        pr = -1
        for i in range(r, rows):
            if A[i, c] % p != 0:
                pr = i
                break
        if pr == -1:
            continue
        A[[r, pr]] = A[[pr, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        for i in range(rows):
            if i != r and A[i, c] % p != 0:
                A[i] = (A[i] - A[i, c] * A[r]) % p
        piv.append(c)
        r += 1
        if r == rows:
            break
    return A, piv


def nullspace_mod_p(A, p):
    A = np.asarray(A, dtype=np.int64)
    if A.size == 0:
        return np.zeros((A.shape[1], 0), dtype=np.int64)
    R, piv = rref_mod_p(A, p)
    rows, cols = A.shape
    pivset = set(piv)
    free = [c for c in range(cols) if c not in pivset]
    cols_out = []
    for f in free:
        v = np.zeros(cols, dtype=np.int64)
        v[f] = 1
        for i, pc in enumerate(piv):
            v[pc] = (-R[i, f]) % p
        cols_out.append(v % p)
    if not cols_out:
        return np.zeros((cols, 0), dtype=np.int64)
    return np.stack(cols_out, axis=1)


def colspace_basis_mod_p(A, p):
    A = np.asarray(A, dtype=np.int64)
    if A.size == 0:
        return np.zeros((A.shape[0], 0), dtype=np.int64)
    R, piv = rref_mod_p(A, p)
    if not piv:
        return np.zeros((A.shape[0], 0), dtype=np.int64)
    return (A[:, piv] % p)


def solve_mod_p(A, y, p):
    """solve A x = y over F_p; A full column rank, y in colspan. returns x or None."""
    A = np.asarray(A, dtype=np.int64) % p
    y = (np.asarray(y, dtype=np.int64) % p).reshape(-1, 1)
    N, c = A.shape
    M = np.concatenate([A, y], axis=1)
    R, piv = rref_mod_p(M, p)
    if c in piv:        # pivot in augmented column => inconsistent
        return None
    x = np.zeros(c, dtype=np.int64)
    for i, j in enumerate(piv):
        if j < c:
            x[j] = R[i, -1] % p
    return x


def quotient_induced(Sigma, B_im, B_ker, p):
    """
    Matrix of endomorphism Sigma on the subquotient colspan(B_ker)/colspan(B_im),
    where colspan(B_im) <= colspan(B_ker) and Sigma preserves both.
    Columns of B_im, B_ker are vectors in the ambient space.
    """
    Sigma = np.asarray(Sigma, dtype=np.int64) % p
    nim = B_im.shape[1] if B_im.size else 0
    if B_ker.size == 0 or B_ker.shape[1] == 0:
        return np.zeros((0, 0), dtype=np.int64), 0
    M = np.concatenate([B_im, B_ker], axis=1) if nim else B_ker.copy()
    R, piv = rref_mod_p(M, p)
    reps_idx = [c - nim for c in piv if c >= nim]
    reps = B_ker[:, reps_idx]
    d = reps.shape[1]
    if d == 0:
        return np.zeros((0, 0), dtype=np.int64), 0
    Bbasis = np.concatenate([B_im, reps], axis=1) if nim else reps.copy()
    Mind = np.zeros((d, d), dtype=np.int64)
    for j in range(d):
        yv = (Sigma @ reps[:, j]) % p
        x = solve_mod_p(Bbasis, yv, p)
        if x is None:
            raise RuntimeError("Sigma does not preserve the subspace -- bug or sigma not a chain map")
        Mind[:, j] = x[nim:nim + d] % p
    return Mind, d


# ====================================================================
# sigma-action matrices on chains C_n and cochains C^n
# ====================================================================
def sigma_chain_matrix(alg, n, S, p):
    """matrix of sigma on C_n = A (x) \\bar A^{(x)n}; S = matrix of sigma on A (columns=images)."""
    m = alg.m
    t = alg.t
    R = alg.R
    basis = cn_basis(alg, n)
    index = {g: i for i, g in enumerate(basis)}
    dim = len(basis)
    Mat = np.zeros((dim, dim), dtype=np.int64)
    Sp = S % p
    # reduced action Sbar = sigma descended to Abar = A/k1: keep j in R, dropping the
    # unit coordinate. Well-defined for ANY algebra automorphism (sigma(1)=1); the
    # "sigma preserves the radical" reading is only literally true for local algebras.
    for col, gen in enumerate(basis):
        a0 = gen[0]
        rs = gen[1:]
        # image of e_{a0}: sum_i Sp[i,a0] e_i  (full)
        i_nz = np.nonzero(Sp[:, a0])[0]
        # images of reduced factors: list of (index_in_R_value, coeff) per position
        factor_options = []
        ok = True
        for r in rs:
            js = [(j, int(Sp[j, r])) for j in R if Sp[j, r] % p != 0]
            if not js:
                ok = False
                break
            factor_options.append(js)
        if not ok:
            continue
        for i in i_nz:
            ci = int(Sp[i, a0])
            # iterate cartesian product of factor options
            stack = [(0, (int(i),), ci)]
            # do it iteratively
            partials = [((int(i),), ci)]
            for opts in factor_options:
                newp = []
                for (word, coeff) in partials:
                    for (j, cj) in opts:
                        newp.append((word + (j,), (coeff * cj) % p))
                partials = newp
            for (word, coeff) in partials:
                Mat[index[word], col] = (Mat[index[word], col] + coeff) % p
    return Mat, basis, index


def sigma_cochain_matrix(alg, n, S, Sinv, p):
    """matrix of sigma on C^n = Hom(\\bar A^{(x)n}, A): (s.f)(u)=s(f(s^{-1}u))."""
    m = alg.m
    R = alg.R
    basis = cochain_basis(alg, n)
    index = {g: i for i, g in enumerate(basis)}
    dim = len(basis)
    Mat = np.zeros((dim, dim), dtype=np.int64)
    Sp = S % p
    Sinvp = Sinv % p
    import itertools
    u_words = [()] if n == 0 else list(itertools.product(R, repeat=n))
    for col, (w, jj) in enumerate(basis):
        # image: sum_u sum_i [ prod_k Sinv[w_k, u_k] * S[i, jj] ] (u, i)
        i_nz = [(i, int(Sp[i, jj])) for i in range(m) if Sp[i, jj] % p != 0]
        for u in u_words:
            coeff_w = 1
            for k in range(n):
                coeff_w = (coeff_w * int(Sinvp[w[k], u[k]])) % p
            if coeff_w == 0:
                continue
            for (i, sij) in i_nz:
                val = (coeff_w * sij) % p
                if val:
                    Mat[index[(u, i)], col] = (Mat[index[(u, i)], col] + val) % p
    return Mat, basis, index


# ====================================================================
# induced action on HH_n (homology) and HH^n (cohomology)
# ====================================================================
def induced_on_HH_homology(alg, n, S, p, resolution=None):
    """returns (matrix of sigma on HH_n, dim HH_n)."""
    from quiverlab.engine.resolutions import _default
    resolution = _default(resolution)
    basis_n = resolution.term_basis(alg, n)
    idx_n = {g: i for i, g in enumerate(basis_n)}
    basis_nm1 = resolution.term_basis(alg, n - 1) if n >= 1 else None
    idx_nm1 = {g: i for i, g in enumerate(basis_nm1)} if n >= 1 else None
    basis_np1 = resolution.term_basis(alg, n + 1)
    idx_np1 = {g: i for i, g in enumerate(basis_np1)}
    bn = resolution.differential_matrix(alg, n, basis_n, idx_nm1) if n >= 1 else np.zeros((1, len(basis_n)), dtype=np.int64)
    bnp1 = resolution.differential_matrix(alg, n + 1, basis_np1, idx_n)
    ker = nullspace_mod_p(bn, p) if n >= 1 else np.eye(len(basis_n), dtype=np.int64)
    im = colspace_basis_mod_p(bnp1, p)
    Sigma, _, _ = resolution.sigma_chain_matrix(alg, n, S, p)
    return quotient_induced(Sigma, im, ker, p)


def induced_on_HH_cohomology(alg, n, S, Sinv, p, resolution=None):
    """returns (matrix of sigma on HH^n, dim HH^n)."""
    from quiverlab.engine.resolutions import _default
    resolution = _default(resolution)
    basis_n = resolution.cochain_basis(alg, n)
    idx_n = {g: i for i, g in enumerate(basis_n)}
    basis_np1 = resolution.cochain_basis(alg, n + 1)
    idx_np1 = {g: i for i, g in enumerate(basis_np1)}
    dn = resolution.coboundary_matrix(alg, n, basis_n, idx_np1)
    if n >= 1:
        basis_nm1 = resolution.cochain_basis(alg, n - 1)
        idx_nm1 = {g: i for i, g in enumerate(basis_nm1)}
        dnm1 = resolution.coboundary_matrix(alg, n - 1, basis_nm1, idx_n)
        im = colspace_basis_mod_p(dnm1, p)
    else:
        im = np.zeros((len(basis_n), 0), dtype=np.int64)
    ker = nullspace_mod_p(dn, p)
    Sigma, _, _ = resolution.sigma_cochain_matrix(alg, n, S, Sinv, p)
    return quotient_induced(Sigma, im, ker, p)


def is_identity(M, p):
    """True if the square integer matrix M reduces to the identity modulo p.

    The empty (0x0) matrix counts as the identity -- the induced action on a
    zero-dimensional HH group is vacuously trivial.
    """
    if M.shape[0] == 0:
        return True
    return np.array_equal(M % p, np.eye(M.shape[0], dtype=np.int64))


# ====================================================================
# automorphisms
# ====================================================================
def nakayama_quantum_ci(q, p):
    """nu for Lambda_q = k<x,y>/(x^2,y^2,yx-q xy): diag(1, q^{-1}, q, 1) on {1,x,y,xy}."""
    qinv = pow(q % p, p - 2, p)
    S = np.diag([1, qinv % p, q % p, 1]).astype(np.int64)
    Sinv = np.diag([1, q % p, qinv % p, 1]).astype(np.int64)
    return S, Sinv


def inverse_mod_p(M, p):
    """Inverse of a square matrix over F_p (raises if singular).  Deterministic
    Gauss-Jordan via the same `rref_mod_p` used everywhere else."""
    M = np.asarray(M, dtype=np.int64) % p
    n = M.shape[0]
    aug = np.concatenate([M, np.eye(n, dtype=np.int64)], axis=1)
    R, piv = rref_mod_p(aug, p)
    if piv[:n] != list(range(n)):
        raise ValueError("matrix is singular mod %d" % p)
    return R[:, n:] % p


def frobenius_form(alg, p, tries=200):
    """A Frobenius form for `alg` over F_p, or None if the algebra is not Frobenius.

    A finite-dimensional algebra A is *Frobenius* iff there is a linear functional
    lambda: A -> k whose associated bilinear form (a,b) |-> lambda(a*b) is
    non-degenerate.  This returns the *covector* lambda (length-m int64, f-basis
    coords) together with its Gram matrix G_{ij} = lambda(e_i e_j); `None` if no
    non-degenerate lambda is found (so the algebra is not Frobenius / not
    self-injective).  Candidates are tried deterministically: the coordinate
    functionals e_c^* first (these hit the socle for the local validation algebras),
    then the all-ones functional, then seeded pseudo-random integer covectors."""
    m = alg.m

    def gram(lam):
        G = np.zeros((m, m), dtype=np.int64)
        for i in range(m):
            for j in range(m):
                G[i, j] = int(np.dot(lam, alg.mult_full(i, j))) % p
        return G % p

    cands = [np.eye(m, dtype=np.int64)[c] for c in range(m)]
    cands.append(np.ones(m, dtype=np.int64))
    rng = np.random.default_rng(0)
    for _ in range(tries):
        cands.append(rng.integers(1, p, size=m).astype(np.int64))
    for lam in cands:
        G = gram(lam)
        _, piv = rref_mod_p(G, p)
        if len(piv) == m:
            return lam % p, G
    return None


def is_frobenius(alg, p):
    """True iff `alg` admits a non-degenerate Frobenius form over F_p (== Frobenius
    == self-injective, for a finite-dimensional algebra)."""
    return frobenius_form(alg, p) is not None


def nakayama_automorphism(alg, p):
    """The Nakayama automorphism nu of a Frobenius algebra, computed from its
    structure, as `(S, S_inv)` of integer matrices over F_p (f-basis, columns =
    images) -- the drop-in structural replacement for `nakayama_quantum_ci`.

    For a Frobenius algebra with form lambda and Gram matrix G_{ij} = lambda(e_i e_j),
    the Nakayama automorphism is defined by  lambda(ab) = lambda(b nu(a)); writing it
    as a matrix N (columns = images) this forces  N = G^{-1} G^T  (so N = id exactly
    when G is symmetric, i.e. A is *symmetric*).  The DA-bimodule of A is {}_1 A_nu,
    so this is the self-injective case of the dualizing bimodule (PLAN item A) -- the
    induced Theta on HH_* is `induced_on_HH_homology(alg, n, S, p)`.

    The Nakayama automorphism is canonical only up to an inner automorphism (a
    different Frobenius form gives an inner-conjugate nu), but inner automorphisms act
    trivially on Hochschild (co)homology, so the induced Theta -- and its Coxeter
    characteristic polynomial -- is independent of the choice.

    Raises ValueError if `alg` is not Frobenius (no non-degenerate form exists)."""
    found = frobenius_form(alg, p)
    if found is None:
        raise ValueError(
            "algebra %r is not Frobenius over F_%d (no non-degenerate Frobenius "
            "form), so it has no Nakayama automorphism -- the genuine "
            "non-self-injective Theta is the open PLAN item A frontier."
            % (getattr(alg, "name", "<algebra>"), p))
    _, G = found
    Ginv = inverse_mod_p(G, p)
    S = (Ginv @ (G.T % p)) % p
    Sinv = inverse_mod_p(S, p)
    return S, Sinv


def identity_auto(m, p):
    """The identity automorphism of an m-dimensional algebra as `(S, S_inv)`.

    Returns the pair of m x m identity matrices expected by the induced-action
    engine; used as the trivial-automorphism baseline in the validation checks
    (sigma = id must induce the identity on every HH_n and HH^n).
    """
    return np.eye(m, dtype=np.int64), np.eye(m, dtype=np.int64)
