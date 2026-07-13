# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
"""Connes' B operator and cyclic homology HC_*(A)  (PLAN item D).

This completes the Tamarkin-Tsygan / Connes side of the calculus begun in
`tt_calculus`: the cup ring and cap module action live there, here we add the
**Connes boundary** `B : C_n -> C_{n+1}` on the normalized bar complex and the
**cyclic homology** `HC_*(A)` it computes via the (b, B) bicomplex.

Why this matters for Han.  Cyclic homology is the home of the *proven* base-case
mechanism: Bergh-Madsen / Igusa express the Euler characteristic and the
Cartan-determinant obstruction in cyclic terms, and the Koszul / local / cellular
cases of Han's conjecture are settled there.  Computing `HC_*` directly connects
the engine to that proof machinery (PLAN item D).

Conventions (matching `hh_engine`).

  * Chains.  C_n = A (x) \\bar A^{(x)n}; basis (a_0, r_1, ..., r_n) with a_0 in
    0..m-1 and r_i in R (the reduced basis = "drop the unit coordinate t").
  * Hochschild b.  The cyclic Hochschild differential `differential_matrix`
    (b: C_n -> C_{n-1}); b^2 = 0.
  * Connes B.  On the normalized complex,
        B(a_0 (x) a_1 (x) ... (x) a_n)
            = sum_{i=0}^{n} (-1)^{n i} 1 (x) a_i (x) ... (x) a_n (x) a_0 (x) ... (x) a_{i-1},
    i.e. insert the unit 1_A in the A-slot and cyclically rotate all n+1 entries
    (with the original a_0 now reduced mod k.1, so any term carrying the unit in a
    bar slot dies).  This is the standard `B = (1-t) s N` made explicit on the
    normalized complex.  It satisfies B^2 = 0 and the mixed-complex identity
    bB + Bb = 0, so (C_*, b, B) is a mixed complex.
  * Cyclic homology.  HC_n = H_n of the total complex of the first-quadrant
    (b, B) bicomplex:  Tot_n = C_n (+) C_{n-2} (+) C_{n-4} (+) ... with total
    differential D = b + B.  Then HC_n = dim Tot_n - rank D_n - rank D_{n+1}.

The mixed-complex identities (b^2 = 0, B^2 = 0, bB + Bb = 0) are asserted directly
on the matrices in tests/test_cyclic_homology.py -- they are the rigorous oracle
that the (b, B) bicomplex is valid; HC_*(k) = [1,0,1,0,...] and the separable
HC_*(k x k) = [2,0,2,0,...] (periodicity pushing HH_0 into even degrees) pin the
specific Connes convention, and HC_0 = HH_0 is checked across the zoo.
"""

import numpy as np

from quiverlab.engine.hh_engine import (
    cn_basis,
    differential_matrix,
    rank_mod_p,
)


# ---------------------------------------------------------------------------
# Connes' B operator  B : C_n -> C_{n+1}
# ---------------------------------------------------------------------------
def connes_B_matrix(alg, n, basis_n=None, index_np1=None):
    """Matrix of Connes' boundary B : C_n -> C_{n+1}.  Shape (dim C_{n+1}, dim C_n).

    Integer (un-reduced) entries -- the engine reduces mod p downstream.
    """
    t = alg.t
    if basis_n is None:
        basis_n = cn_basis(alg, n)
    basis_np1 = cn_basis(alg, n + 1)
    if index_np1 is None:
        index_np1 = {g: i for i, g in enumerate(basis_np1)}
    rows = len(basis_np1)
    cols = len(basis_n)
    M = np.zeros((rows, cols), dtype=np.int64)
    for cidx, gen in enumerate(basis_n):
        entries = gen  # (a_0, r_1, ..., r_n), length n+1
        for i in range(0, n + 1):
            rotated = entries[i:] + entries[:i]   # cyclic shift, length n+1
            # the only entry that can be the unit is the original a_0; if it lands
            # in a bar slot it is killed in A/k.1 -> the whole term vanishes.
            if t in rotated:
                continue
            key = (t,) + rotated            # new A-slot = unit, bar slots = rotated
            sign = 1 if (n * i) % 2 == 0 else -1
            M[index_np1[key], cidx] += sign
    return M


# ---------------------------------------------------------------------------
# The (b, B) bicomplex total complex and cyclic homology
# ---------------------------------------------------------------------------
def _tot_degrees(n):
    """Chain degrees appearing in Tot_n = C_n (+) C_{n-2} (+) ... (descending)."""
    return list(range(n, -1, -2))


def _total_differential(alg, n, bmats, Bmats, dims):
    """Block matrix of D = b + B : Tot_n -> Tot_{n-1}, over the integers.

    bmats[k]   = b : C_k -> C_{k-1}   (shape (dim C_{k-1}, dim C_k))
    Bmats[k]   = B : C_k -> C_{k+1}   (shape (dim C_{k+1}, dim C_k))
    dims[k]    = dim C_k
    """
    src = _tot_degrees(n)         # column blocks
    tgt = _tot_degrees(n - 1)     # row blocks
    if not tgt:
        return np.zeros((0, sum(dims[d] for d in src)), dtype=np.int64)
    row_off = {}
    off = 0
    for d in tgt:
        row_off[d] = off
        off += dims[d]
    nrows = off
    col_off = {}
    off = 0
    for d in src:
        col_off[d] = off
        off += dims[d]
    ncols = off
    D = np.zeros((nrows, ncols), dtype=np.int64)
    for d in src:
        c0 = col_off[d]
        # b : C_d -> C_{d-1}  (same bicomplex column, lands in Tot_{n-1} block d-1)
        if d - 1 in row_off and d >= 1:
            r0 = row_off[d - 1]
            blk = bmats[d]
            D[r0:r0 + blk.shape[0], c0:c0 + blk.shape[1]] += blk
        # B : C_d -> C_{d+1}  (lands in Tot_{n-1} block d+1)
        if d + 1 in row_off:
            r0 = row_off[d + 1]
            blk = Bmats[d]
            D[r0:r0 + blk.shape[0], c0:c0 + blk.shape[1]] += blk
    return D


def cyclic_homology_dims(alg, N, primes=(32003, 2, 3, 5)):
    """Return dict p -> [dim HC_0, ..., dim HC_N] via the (b, B) bicomplex.

    Builds the bar complex up to C_{N+1} (the usual dim C_n = m*(m-1)^n blow-up
    caps N); HC_n = dim Tot_n - rank D_n - rank D_{n+1}.
    """
    maxdeg = N + 1
    bases = {k: cn_basis(alg, k) for k in range(0, maxdeg + 1)}
    indices = {k: {g: i for i, g in enumerate(bases[k])} for k in range(0, maxdeg + 1)}
    dims = {k: len(bases[k]) for k in range(0, maxdeg + 1)}
    bmats = {}
    Bmats = {}
    for k in range(0, maxdeg + 1):
        if k >= 1:
            bmats[k] = differential_matrix(alg, k, bases[k], indices[k - 1])
        else:
            bmats[k] = np.zeros((0, dims[0]), dtype=np.int64)
        if k <= maxdeg - 1:
            Bmats[k] = connes_B_matrix(alg, k, bases[k], indices[k + 1])
        else:
            Bmats[k] = np.zeros((dims.get(k + 1, 0), dims[k]), dtype=np.int64)
    Dmats = {n: _total_differential(alg, n, bmats, Bmats, dims) for n in range(0, N + 2)}
    out = {}
    for p in primes:
        ranks = {n: rank_mod_p(Dmats[n], p) for n in range(0, N + 2)}
        dimsTot = {n: sum(dims[d] for d in _tot_degrees(n)) for n in range(0, N + 1)}
        out[p] = [dimsTot[n] - ranks[n] - ranks[n + 1] for n in range(0, N + 1)]
    return out


# ---------------------------------------------------------------------------
# Mixed-complex identity checks (the rigorous descent oracle)
# ---------------------------------------------------------------------------
def check_mixed_complex(alg, N, prime=32003):
    """Verify b^2 = 0, B^2 = 0 and bB + Bb = 0 on C_0..C_N over F_p.

    Returns the first failing (identity, degree) or None if all hold.
    """
    bases = {k: cn_basis(alg, k) for k in range(0, N + 3)}
    idx = {k: {g: i for i, g in enumerate(bases[k])} for k in range(0, N + 3)}
    b = {k: differential_matrix(alg, k, bases[k], idx[k - 1]) for k in range(1, N + 3)}
    B = {k: connes_B_matrix(alg, k, bases[k], idx[k + 1]) for k in range(0, N + 2)}
    for n in range(0, N + 1):
        # b^2 : C_{n+1} -> C_{n-1}
        if n >= 1:
            bb = (b[n] @ b[n + 1]) % prime
            if np.any(bb):
                return ("b^2", n)
        # B^2 : C_n -> C_{n+2}
        BB = (B[n + 1] @ B[n]) % prime
        if np.any(BB):
            return ("B^2", n)
        # bB + Bb : C_n -> C_n
        bB = b[n + 1] @ B[n]
        Bb = (B[n - 1] @ b[n]) if n >= 1 else np.zeros_like(bB)
        if np.any((bB + Bb) % prime):
            return ("bB+Bb", n)
    return None
