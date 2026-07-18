# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
"""The Tamarkin-Tsygan calculus substrate: cup product on HH^* and the cap-product
module action HH^* (x) HH_* -> HH_* (PLAN item A step 3 / item D).

PLAN item A's headline is the *genuine* Auslander-Reiten action Theta on
non-self-injective algebras.  Following Armenta, *On the cap product in Hochschild
theory* (arXiv:1908.02255), that action is realised by **capping with the class of
the inverse dualizing complex** -- so the cup ring and the cap module action are the
substrate the genuine Theta is built on.  This module supplies them on the normalized
bar complex (the reference oracle), with the products computed exactly over the
integer structure constants and reduced mod p only at the linear-algebra step.

Conventions (matching `hh_engine` homology and `scan3` cohomology).

  * Cochains.  C^p = Hom(\\bar A^{(x)p}, A); the engine basis element (w, j) with
    w in R^p, j in 0..m-1, is the functional sending the reduced word w to the
    basis vector e_j (and every other word to 0).  A cochain is thus a map
    w |-> f(w) in A.
  * Cup product.  (f ⌣ g)(a_1,...,a_{p+q}) = f(a_1,...,a_p) · g(a_{p+1},...,a_{p+q}),
    the A-product of the two values.  Associative on the nose; the unit is the
    0-cochain 1_A in C^0 = A (= e_t in the unit-adapted basis).  Leibniz w.r.t. the
    `scan3` coboundary:  δ(f ⌣ g) = δf ⌣ g + (-1)^p f ⌣ δg, so ⌣ descends to a
    graded-commutative product on HH^*.
  * Chains.  C_n = A (x) \\bar A^{(x)n}; basis (a_0, r_1,...,r_n), a_0 in 0..m-1,
    r_i in R.
  * Cap product.  C^p (x) C_n -> C_{n-p},
        f ∩ (a_0 (x) a_1 (x) ... (x) a_n) = a_0·f(a_1,...,a_p) (x) a_{p+1} (x) ... (x) a_n.
    The unit 1_A in C^0 acts as the identity; and (with this single-sided
    convention) the module law is  f ∩ (g ∩ z) = (g ⌣ f) ∩ z  (no sign), making
    HH_* a graded module over the cup ring HH^*.

Every identity above is asserted directly on basis (co)chains in
tests/test_tt_calculus.py -- those nose-level identities are the rigorous oracle that
the products descend to (co)homology; a handful of induced structure constants are
frozen there as regression pins.
"""

import itertools

import numpy as np

from quiverlab.engine.hh_engine import cn_basis, differential_matrix
from quiverlab.engine.scan3 import cochain_basis, coboundary_matrix
from quiverlab.engine.coxeter import (
    nullspace_mod_p,
    colspace_basis_mod_p,
    solve_mod_p,
    rref_mod_p,
)


# ---------------------------------------------------------------------------
# A-arithmetic and cochain <-> dict conversion
# ---------------------------------------------------------------------------
def _mul(alg, u, v):
    """Product u·v in A of two full length-m int64 vectors (exact, un-reduced)."""
    out = np.zeros(alg.m, dtype=np.int64)
    for a in np.nonzero(u)[0]:
        ua = u[a]
        for b in np.nonzero(v)[0]:
            out += ua * v[b] * alg.mult_full(int(a), int(b))
    return out


def cochain_to_dict(alg, p, vec):
    """Cochain vector (over `cochain_basis(alg, p)`) -> dict {word: A-vector}."""
    m = alg.m
    d = {}
    for c, (w, j) in zip(vec, cochain_basis(alg, p)):
        if c == 0:
            continue
        if w not in d:
            d[w] = np.zeros(m, dtype=np.int64)
        d[w][j] += int(c)
    return d


def dict_to_cochain(alg, p, d):
    """dict {word: A-vector} -> cochain vector over `cochain_basis(alg, p)`."""
    basis = cochain_basis(alg, p)
    idx = {gen: i for i, gen in enumerate(basis)}
    v = np.zeros(len(basis), dtype=np.int64)
    for w, av in d.items():
        for j in np.nonzero(av)[0]:
            v[idx[(w, int(j))]] += int(av[j])
    return v


# ---------------------------------------------------------------------------
# Cup product on cochains:  C^p (x) C^q -> C^{p+q}
# ---------------------------------------------------------------------------
def cup_cochain(alg, p, q, f, g):
    """Vector of f ⌣ g over `cochain_basis(alg, p+q)`.  Integer (un-reduced)."""
    fd = cochain_to_dict(alg, p, f)
    gd = cochain_to_dict(alg, q, g)
    basis_pq = cochain_basis(alg, p + q)
    idx = {gen: i for i, gen in enumerate(basis_pq)}
    out = np.zeros(len(basis_pq), dtype=np.int64)
    for wf, fv in fd.items():
        for wg, gv in gd.items():
            prod = _mul(alg, fv, gv)
            w = wf + wg
            for jj in np.nonzero(prod)[0]:
                out[idx[(w, int(jj))]] += int(prod[jj])
    return out


def unit_cochain(alg):
    """The cup unit 1_A in C^0 = A (= e_t in the unit-adapted basis)."""
    basis0 = cochain_basis(alg, 0)            # [((), j) for j in 0..m-1]
    v = np.zeros(len(basis0), dtype=np.int64)
    v[basis0.index(((), alg.t))] = 1
    return v


# ---------------------------------------------------------------------------
# Cap product:  C^p (x) C_n -> C_{n-p}
# ---------------------------------------------------------------------------
def cap_cochain(alg, p, n, f, z):
    """Vector of f ∩ z over `cn_basis(alg, n-p)`.  f in C^p, z in C_n.  Integer."""
    if p > n:
        raise ValueError("cap product needs p <= n (C^%d ∩ C_%d)" % (p, n))
    fd = cochain_to_dict(alg, p, f)
    basis_n = cn_basis(alg, n)
    basis_out = cn_basis(alg, n - p)
    idx = {gen: i for i, gen in enumerate(basis_out)}
    out = np.zeros(len(basis_out), dtype=np.int64)
    for c, gen in zip(z, basis_n):
        if c == 0:
            continue
        a0 = gen[0]
        rs = gen[1:]
        wp = rs[:p]
        rest = rs[p:]
        fv = fd.get(wp)
        if fv is None:
            continue
        # A-slot = e_{a0} · f(wp)
        slot = np.zeros(alg.m, dtype=np.int64)
        for b in np.nonzero(fv)[0]:
            slot += int(fv[b]) * alg.mult_full(int(a0), int(b))
        for a0p in np.nonzero(slot)[0]:
            key = (int(a0p),) + rest
            out[idx[key]] += int(c) * int(slot[a0p])
    return out


# ---------------------------------------------------------------------------
# Representatives of HH^n and HH_n over F_p (reusing the coxeter primitives)
# ---------------------------------------------------------------------------
class _Quotient:
    """A basis of representatives for a subquotient colspan(ker)/colspan(im) over F_p,
    with a method to read off the class coordinates of any element of `ker`."""

    def __init__(self, im, ker, p):
        self.p = p
        self.nim = im.shape[1] if im.size else 0
        if ker.size == 0 or ker.shape[1] == 0:
            self.reps = np.zeros((ker.shape[0], 0), dtype=np.int64)
            self.dim = 0
            self.basis = im.copy() if self.nim else np.zeros((ker.shape[0], 0), dtype=np.int64)
            return
        M = np.concatenate([im, ker], axis=1) if self.nim else ker.copy()
        _, piv = rref_mod_p(M, p)
        reps_idx = [c - self.nim for c in piv if c >= self.nim]
        self.reps = ker[:, reps_idx]
        self.dim = self.reps.shape[1]
        self.basis = (np.concatenate([im, self.reps], axis=1)
                      if self.nim else self.reps.copy())

    def coords(self, v):
        """Class coordinates (length self.dim) of v, which must lie in colspan(ker)."""
        if self.dim == 0:
            return np.zeros(0, dtype=np.int64)
        x = solve_mod_p(self.basis, np.asarray(v, dtype=np.int64) % self.p, self.p)
        if x is None:
            raise RuntimeError("element is not a (co)cycle representative -- "
                               "the product failed to descend to (co)homology")
        return x[self.nim:self.nim + self.dim] % self.p


def cohomology_classes(alg, n, p):
    """`_Quotient` for HH^n(A) over F_p (cocycles mod coboundaries, bar complex)."""
    bn = cochain_basis(alg, n)
    idx_np1 = {g: i for i, g in enumerate(cochain_basis(alg, n + 1))}
    dn = coboundary_matrix(alg, n, bn, idx_np1)
    ker = nullspace_mod_p(dn, p)
    if n >= 1:
        bnm1 = cochain_basis(alg, n - 1)
        idx_n = {g: i for i, g in enumerate(bn)}
        dnm1 = coboundary_matrix(alg, n - 1, bnm1, idx_n)
        im = colspace_basis_mod_p(dnm1, p)
    else:
        im = np.zeros((len(bn), 0), dtype=np.int64)
    return _Quotient(im, ker, p)


def homology_classes(alg, n, p):
    """`_Quotient` for HH_n(A) over F_p (cycles mod boundaries, bar complex)."""
    bn = cn_basis(alg, n)
    if n >= 1:
        idx_nm1 = {g: i for i, g in enumerate(cn_basis(alg, n - 1))}
        b_n = differential_matrix(alg, n, bn, idx_nm1)
        ker = nullspace_mod_p(b_n, p)
    else:
        ker = np.eye(len(bn), dtype=np.int64)
    idx_n = {g: i for i, g in enumerate(bn)}
    b_np1 = differential_matrix(alg, n + 1, cn_basis(alg, n + 1), idx_n)
    im = colspace_basis_mod_p(b_np1, p)
    return _Quotient(im, ker, p)


# ---------------------------------------------------------------------------
# Induced products on (co)homology, as structure-constant tensors over F_p
# ---------------------------------------------------------------------------
def cup_product_matrix(alg, p, q, prime):
    """Structure constants of the cup product HH^p (x) HH^q -> HH^{p+q} over F_p.

    Returns (C, dp, dq, dpq) where C has shape (dpq, dp, dq): the class of
    rep_i ⌣ rep_j in HH^{p+q} is sum_k C[k,i,j]·basis_k.  Asserts each cup of
    cocycles is again a cocycle (the descent check)."""
    Hp = cohomology_classes(alg, p, prime)
    Hq = cohomology_classes(alg, q, prime)
    Hpq = cohomology_classes(alg, p + q, prime)
    C = np.zeros((Hpq.dim, Hp.dim, Hq.dim), dtype=np.int64)
    for i in range(Hp.dim):
        for j in range(Hq.dim):
            h = cup_cochain(alg, p, q, Hp.reps[:, i], Hq.reps[:, j])
            C[:, i, j] = Hpq.coords(h)
    return C, Hp.dim, Hq.dim, Hpq.dim


def cap_product_matrix(alg, p, n, prime):
    """Structure constants of the cap action HH^p (x) HH_n -> HH_{n-p} over F_p.

    Returns (C, dp, dn, dnp) where C has shape (dnp, dp, dn): the class of
    rep_i ∩ rep_j in HH_{n-p} is sum_k C[k,i,j]·basis_k.  Asserts each cap of a
    cocycle with a cycle is again a cycle (the descent check)."""
    Hp = cohomology_classes(alg, p, prime)
    Hn = homology_classes(alg, n, prime)
    Hnp = homology_classes(alg, n - p, prime)
    C = np.zeros((Hnp.dim, Hp.dim, Hn.dim), dtype=np.int64)
    for i in range(Hp.dim):
        for j in range(Hn.dim):
            w = cap_cochain(alg, p, n, Hp.reps[:, i], Hn.reps[:, j])
            C[:, i, j] = Hnp.coords(w)
    return C, Hp.dim, Hn.dim, Hnp.dim


# ---------------------------------------------------------------------------
# Gerstenhaber bracket on HH^* (the Lie side of the TT calculus, PLAN item D)
# ---------------------------------------------------------------------------
#
# Circle product  f bar{o} g : C^p (x) C^q -> C^{p+q-1},
#     (f o g)(a_1,...,a_{p+q-1})
#         = sum_{i=0}^{p-1} (-1)^{(q-1) i}
#               f(a_1,...,a_i, g(a_{i+1},...,a_{i+q}), a_{i+q+1},...,a_{p+q-1}),
# where the value g(...) in A is reduced mod k.1 before it is fed to f (f is a
# normalized cochain, so it kills any unit component anyway).  The Gerstenhaber
# bracket is the graded commutator of the circle product:
#     [f, g] = f o g - (-1)^{(p-1)(q-1)} g o f.
#
# It is graded-antisymmetric, descends to a graded Lie bracket on HH^* (degree -1),
# restricts on HH^1 to the commutator of (outer) derivations, and the multiplication
# 2-cochain m satisfies [m, m] = 0 (<=> associativity).  Every one of these is pinned
# in tests/test_gerstenhaber.py; descent is checked there directly (the bracket of
# cocycles is a cocycle; the bracket of a coboundary with a cocycle is a coboundary).
def circle_cochain(alg, p, q, f, g):
    """Vector of the Gerstenhaber circle product f o g over `cochain_basis(p+q-1)`.

    Integer (un-reduced).  Requires p >= 1; q >= 0."""
    if p < 1:
        raise ValueError("circle product needs p >= 1")
    t = alg.t
    m = alg.m
    fd = cochain_to_dict(alg, p, f)
    gd = cochain_to_dict(alg, q, g)
    out = {}
    for u in itertools.product(alg.R, repeat=p + q - 1):
        acc = np.zeros(m, dtype=np.int64)
        for i in range(0, p):
            gval = gd.get(u[i:i + q])
            if gval is None:
                continue
            sign = 1 if ((q - 1) * i) % 2 == 0 else -1
            for c in np.nonzero(gval)[0]:
                if c == t:            # reduce g's value mod k.1 (f is normalized)
                    continue
                fval = fd.get(u[:i] + (int(c),) + u[i + q:])
                if fval is None:
                    continue
                acc = acc + sign * int(gval[c]) * fval
        if np.any(acc):
            out[u] = acc
    return dict_to_cochain(alg, p + q - 1, out)


def gerstenhaber_bracket_cochain(alg, p, q, f, g):
    """Vector of the Gerstenhaber bracket [f, g] over `cochain_basis(p+q-1)`.

    [f, g] = f o g - (-1)^{(p-1)(q-1)} g o f.  Integer (un-reduced)."""
    fg = circle_cochain(alg, p, q, f, g)
    gf = circle_cochain(alg, q, p, g, f)
    sign = 1 if ((p - 1) * (q - 1)) % 2 == 0 else -1
    return fg - sign * gf


def multiplication_cochain(alg):
    """The multiplication 2-cochain m in C^2: m(r_1, r_2) = r_1 . r_2 (in A)."""
    out = {}
    for (r1, r2) in itertools.product(alg.R, repeat=2):
        prod = alg.mult_full(r1, r2)
        if np.any(prod):
            out[(r1, r2)] = prod.astype(np.int64)
    return dict_to_cochain(alg, 2, out)


def gerstenhaber_bracket_matrix(alg, p, q, prime):
    """Structure constants of the induced bracket HH^p (x) HH^q -> HH^{p+q-1} over F_p.

    Returns (C, dp, dq, dpq) where C has shape (dpq, dp, dq): the class of
    [rep_i, rep_j] in HH^{p+q-1} is sum_k C[k,i,j]·basis_k.  Asserts each bracket of
    cocycles is again a cocycle (the descent check, via `cohomology_classes.coords`)."""
    Hp = cohomology_classes(alg, p, prime)
    Hq = cohomology_classes(alg, q, prime)
    Hpq = cohomology_classes(alg, p + q - 1, prime)
    C = np.zeros((Hpq.dim, Hp.dim, Hq.dim), dtype=np.int64)
    for i in range(Hp.dim):
        for j in range(Hq.dim):
            h = gerstenhaber_bracket_cochain(alg, p, q, Hp.reps[:, i], Hq.reps[:, j])
            C[:, i, j] = Hpq.coords(h)
    return C, Hp.dim, Hq.dim, Hpq.dim
