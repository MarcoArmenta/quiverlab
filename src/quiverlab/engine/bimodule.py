# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
"""Hochschild (co)homology with coefficients in an A-bimodule (PLAN item A step 2).

The genuine non-self-injective Auslander-Reiten action Theta is built from the
**dualizing bimodule** `DA = Hom_k(A,k)`, not an algebra automorphism.  This module
adds the substrate that machinery needs: Hochschild homology `HH_*(A, M)` and
cohomology `HH^*(A, M)` with coefficients in an arbitrary A-bimodule `M`, via the
normalized bar complex

    C_n(A, M) = M (x) (A/k1)^{(x)n},     C^n(A, M) = Hom_k((A/k1)^{(x)n}, M),

with the standard differentials (the wrap terms now use the left/right action of A
on M).  For `M = A` (the regular bimodule) this reduces *exactly* to the homology /
cohomology engines (`hh_engine`, `scan3`) -- the reproduction oracle.

The headline bimodule is `DA = Hom_k(A,k)` with `(a·f·b)(x) = f(b x a)`.  It carries
two classical dualities that this module verifies:

    dim HH_n(A, DA) = dim HH^n(A)        and        dim HH^n(A, DA) = dim HH_n(A),

the bar complex for one being the k-linear dual of the bar complex for the other.
These are the bridge the genuine `DA`-twist / cap-product `Theta` is built on; for a
self-injective algebra `DA ≅ {}_1A_ν` recovers the Nakayama picture
(`nakayama_automorphism`).

HARD CONTRACT (as elsewhere): integer (`int64`) matrices, NEVER reduced mod p; the
engine applies `dim HH_n = dim C_n - rank b_n - rank b_{n+1}` (and the cohomology
dual) with `rank_mod_p` per prime.
"""

import itertools
import numpy as np

from quiverlab.engine.hh_engine import rank_mod_p


# ---------------------------------------------------------------------------
# Bimodules
# ---------------------------------------------------------------------------
class Bimodule:
    """A finite-dimensional A-bimodule M (dim `mu`), over the algebra `alg`.

    The actions are stored as int64 tensors of shape (alg.m, mu, mu):
        e_j · m  =  Lact[j] @ m        (left action by the j-th basis vector of A)
        m · e_k  =  Ract[k] @ m        (right action)
    i.e. column `s` of `Lact[j]` is the M-vector `e_j · m_s`.  Built un-reduced;
    `check_bimodule` verifies the module axioms over F_p.
    """

    def __init__(self, alg, mu, Lact, Ract, name=""):
        self.alg = alg
        self.mu = mu
        self.Lact = np.asarray(Lact, dtype=np.int64)
        self.Ract = np.asarray(Ract, dtype=np.int64)
        self.name = name

    def left(self, j, m):
        return self.Lact[j] @ m

    def right(self, k, m):
        return self.Ract[k] @ m


def regular_bimodule(alg):
    """M = A with its own multiplication: e_j·a = e_j a, a·e_k = a e_k.  HH_*(A, A)
    and HH^*(A, A) over this bimodule reproduce the standard engines exactly."""
    m = alg.m
    L = np.zeros((m, m, m), dtype=np.int64)
    R = np.zeros((m, m, m), dtype=np.int64)
    for j in range(m):
        for s in range(m):
            L[j][:, s] = alg.mult_full(j, s)
    for k in range(m):
        for s in range(m):
            R[k][:, s] = alg.mult_full(s, k)
    return Bimodule(alg, m, L, R, name="A")


def dual_bimodule(alg):
    """DA = Hom_k(A,k), the dualizing bimodule, with `(a·f·b)(x) = f(b x a)`.

    In the dual basis {phi_i}, phi_i(e_c) = delta_{ic}:
        (e_j · phi_i)(e_c) = phi_i(e_c e_j)  =>  (e_j·phi_i)[c] = (e_c e_j)[i],
        (phi_i · e_k)(e_c) = phi_i(e_k e_c)  =>  (phi_i·e_k)[c] = (e_k e_c)[i].
    """
    m = alg.m
    L = np.zeros((m, m, m), dtype=np.int64)
    R = np.zeros((m, m, m), dtype=np.int64)
    for j in range(m):
        for i in range(m):
            for c in range(m):
                L[j][c, i] = alg.mult_full(c, j)[i]
    for k in range(m):
        for i in range(m):
            for c in range(m):
                R[k][c, i] = alg.mult_full(k, c)[i]
    return Bimodule(alg, m, L, R, name="DA")


def twisted_bimodule(alg, phi=None, psi=None, name=None):
    """The twisted bimodule {}_phi A_psi: a, b act on x in A by  a·x·b = phi(a) x psi(b).

    `phi`, `psi` are algebra automorphisms of A as m x m matrices in the f-basis
    (columns = images of the basis vectors); `None` means the identity.  This is the
    self-injective case of the dualizing bimodule: for a Frobenius algebra
    `DA ≅ {}_1 A_ν` with ν the Nakayama automorphism (`nakayama_automorphism`), so
    `twisted_bimodule(alg, psi=nu)` realises `DA` and the ν-twisted Hochschild homology
    `HH_*(A, {}_1 A_ν)` whose induced automorphism is the Nakayama action Θ
    (`induced_on_HH_homology`).  The genuine **non**-self-injective Θ (where no such ν
    exists) is the open PLAN item A frontier.
    """
    m = alg.m
    I = np.eye(m, dtype=np.int64)
    phi = I if phi is None else np.asarray(phi, dtype=np.int64)
    psi = I if psi is None else np.asarray(psi, dtype=np.int64)
    L = np.zeros((m, m, m), dtype=np.int64)
    R = np.zeros((m, m, m), dtype=np.int64)
    for j in range(m):                      # e_j · x = phi(e_j) * x
        pj = phi[:, j]
        for s in range(m):
            acc = np.zeros(m, dtype=np.int64)
            for c in np.nonzero(pj)[0]:
                acc += int(pj[c]) * alg.mult_full(int(c), s)
            L[j][:, s] = acc
    for k in range(m):                      # x · e_k = x * psi(e_k)
        pk = psi[:, k]
        for s in range(m):
            acc = np.zeros(m, dtype=np.int64)
            for c in np.nonzero(pk)[0]:
                acc += int(pk[c]) * alg.mult_full(s, int(c))
            R[k][:, s] = acc
    if name is None:
        name = "_{phi}A_{psi}"
    return Bimodule(alg, m, L, R, name=name)


def check_bimodule(M, p=32003):
    """Verify the A-bimodule axioms over F_p: both actions are associative, they
    commute, and the unit acts as the identity.  Returns True/False."""
    alg = M.alg
    m = alg.m
    mu = M.mu
    t = alg.t
    Id = np.eye(mu, dtype=np.int64)
    # unit acts as identity
    if not (np.array_equal(M.Lact[t] % p, Id) and np.array_equal(M.Ract[t] % p, Id)):
        return False
    for i in range(m):
        for j in range(m):
            eij = alg.mult_full(i, j)
            # left: e_i·(e_j·m) = (e_i e_j)·m
            lhs = (M.Lact[i] @ M.Lact[j]) % p
            rhs = sum(int(eij[c]) * M.Lact[c] for c in np.nonzero(eij)[0]) % p \
                if np.any(eij) else np.zeros((mu, mu), dtype=np.int64)
            if not np.array_equal(lhs, rhs % p):
                return False
            # right: (m·e_i)·e_j = m·(e_i e_j)
            lhs = (M.Ract[j] @ M.Ract[i]) % p
            rhs = sum(int(eij[c]) * M.Ract[c] for c in np.nonzero(eij)[0]) % p \
                if np.any(eij) else np.zeros((mu, mu), dtype=np.int64)
            if not np.array_equal(lhs, rhs % p):
                return False
            # commute: e_i·(m·e_j) = (e_i·m)·e_j
            if not np.array_equal((M.Lact[i] @ M.Ract[j]) % p,
                                  (M.Ract[j] @ M.Lact[i]) % p):
                return False
    return True


# ---------------------------------------------------------------------------
# Homology  HH_*(A, M):  C_n = M (x) (A/k1)^{(x)n}
# ---------------------------------------------------------------------------
def cn_basis_coeff(alg, M, n):
    """Basis of C_n(A, M): (s, r_1,...,r_n) with s in 0..mu-1, r_i in R."""
    if n == 0:
        return [(s,) for s in range(M.mu)]
    Rt = alg.R
    return [(s,) + tail for s in range(M.mu)
            for tail in itertools.product(Rt, repeat=n)]


def differential_matrix_coeff(alg, M, n, basis_n, index_nm1):
    """b_n : C_n(A,M) -> C_{n-1}(A,M), shape (dim_{n-1}, dim_n), int64.

    b(m (x) a_1 (x) ... (x) a_n) = (m·a_1) (x) a_2 ... a_n
        + sum_{i=1}^{n-1} (-1)^i m (x) ... a_i a_{i+1} ...
        + (-1)^n (a_n·m) (x) a_1 ... a_{n-1}.
    """
    t = alg.t
    Out = np.zeros((len(index_nm1), len(basis_n)), dtype=np.int64)
    if n == 0:
        return Out
    for col, gen in enumerate(basis_n):
        s = gen[0]
        rs = gen[1:]
        ms = np.zeros(M.mu, dtype=np.int64); ms[s] = 1
        # term 0: (m · a_1) (x) a_2..a_n
        v = M.right(rs[0], ms)
        for c in np.nonzero(v)[0]:
            Out[index_nm1[(int(c),) + rs[1:]], col] += int(v[c])
        # interior terms
        for i in range(1, n):
            prod = alg.mult_full(rs[i - 1], rs[i])
            sign = -1 if (i % 2 == 1) else 1
            for idx in np.nonzero(prod)[0]:
                if idx == t:
                    continue
                merged = rs[:i - 1] + (int(idx),) + rs[i + 1:]
                Out[index_nm1[(s,) + merged], col] += sign * int(prod[idx])
        # term n: (a_n · m) (x) a_1..a_{n-1}
        v = M.left(rs[n - 1], ms)
        sign = -1 if (n % 2 == 1) else 1
        for c in np.nonzero(v)[0]:
            Out[index_nm1[(int(c),) + rs[:n - 1]], col] += sign * int(v[c])
    return Out


def hochschild_homology_with_coefficients(alg, M, N, primes=(32003, 2, 3)):
    """dict p -> [dim HH_0(A,M), ..., dim HH_N(A,M)]."""
    bases = {n: cn_basis_coeff(alg, M, n) for n in range(N + 2)}
    index = {n: {g: i for i, g in enumerate(bases[n])} for n in range(N + 2)}
    Bm = {n: differential_matrix_coeff(alg, M, n, bases[n], index[n - 1])
          for n in range(1, N + 2)}
    out = {}
    for p in primes:
        dims = []
        for n in range(N + 1):
            rn = rank_mod_p(Bm[n], p) if n >= 1 else 0
            rnp1 = rank_mod_p(Bm[n + 1], p)
            dims.append(len(bases[n]) - rn - rnp1)
        out[p] = dims
    return out


# ---------------------------------------------------------------------------
# Cohomology  HH^*(A, M):  C^n = Hom((A/k1)^{(x)n}, M)
# ---------------------------------------------------------------------------
def cochain_basis_coeff(alg, M, n):
    """Basis of C^n(A, M): (word, s), word in R^n, s in 0..mu-1 (the functional
    sending `word` to the s-th basis vector of M)."""
    words = [()] if n == 0 else list(itertools.product(alg.R, repeat=n))
    return [(w, s) for w in words for s in range(M.mu)]


def coboundary_matrix_coeff(alg, M, n, basis_n, index_np1):
    """d^n : C^n(A,M) -> C^{n+1}(A,M), shape (dim_{n+1}, dim_n), int64.

    (d f)(a_1,...,a_{n+1}) = a_1·f(a_2..a_{n+1})
        + sum_{i=1}^n (-1)^i f(a_1.. a_i a_{i+1} ..a_{n+1})
        + (-1)^{n+1} f(a_1..a_n)·a_{n+1}.
    """
    t = alg.t
    Out = np.zeros((len(index_np1), len(basis_n)), dtype=np.int64)
    u_words = list(itertools.product(alg.R, repeat=n + 1))
    for col, (w, s) in enumerate(basis_n):
        ms = np.zeros(M.mu, dtype=np.int64); ms[s] = 1
        for u in u_words:
            out = np.zeros(M.mu, dtype=np.int64)
            # T0: a_1 · f(a_2..a_{n+1})
            if u[1:] == w:
                out = out + M.left(u[0], ms)
            # interior
            for i in range(1, n + 1):
                P = alg.mult_full(u[i - 1], u[i])
                sign = -1 if (i % 2 == 1) else 1
                for c in np.nonzero(P)[0]:
                    if c == t:
                        continue
                    if u[:i - 1] + (int(c),) + u[i + 1:] == w:
                        out = out + sign * int(P[c]) * ms
            # T_{n+1}: f(a_1..a_n) · a_{n+1}
            if u[:n] == w:
                sign = -1 if ((n + 1) % 2 == 1) else 1
                out = out + sign * M.right(u[n], ms)
            for c in np.nonzero(out)[0]:
                Out[index_np1[(u, int(c))], col] += int(out[c])
    return Out


def hochschild_cohomology_with_coefficients(alg, M, N, primes=(32003, 2, 3)):
    """dict p -> [dim HH^0(A,M), ..., dim HH^N(A,M)]."""
    bases = {n: cochain_basis_coeff(alg, M, n) for n in range(N + 2)}
    index = {n: {g: i for i, g in enumerate(bases[n])} for n in range(N + 2)}
    Dm = {n: coboundary_matrix_coeff(alg, M, n, bases[n], index[n + 1])
          for n in range(N + 1)}
    out = {}
    for p in primes:
        dims = []
        for n in range(N + 1):
            rn = rank_mod_p(Dm[n], p)
            rnm1 = rank_mod_p(Dm[n - 1], p) if n >= 1 else 0
            dims.append(len(bases[n]) - rn - rnm1)
        out[p] = dims
    return out
