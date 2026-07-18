# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
"""Numba-compiled hot kernels for the minimal Aᵉ-resolution engine.

These mirror, bit-for-bit, the pure-Python functions in resolutions_minimal.py and
hh_engine.py — which remain the permanent oracle and the fallback when numba is
absent. Dispatch is via USE_KERNELS (read at call time so tests can toggle it).

OVERFLOW: a single product a*b of two reduced residues (each < p) is < p**2. The
nullspace / rank / independent-modulo kernels reduce mod p after every such product, so
they are safe for any prime p < 2**31 (p**2 < ~4.6e18 < 2**63). The radK_kernel and
Dn_kernel matvec, by contrast, accumulate up to m2 products BEFORE reducing once per
output entry, so their bound is m2*(p-1)**2 < 2**63 — e.g. p <= 32003 with m2 <= ~400
gives ~4.1e17, far within int64 for the dim <= ~20 target. Do NOT raise the prime toward
2**31 for those two kernels without reducing the accumulation per step.
"""
import os

import numpy as np

try:
    from numba import njit as _njit, prange
    HAS_NUMBA = True
except Exception:                       # numba/llvmlite missing or broken ABI
    HAS_NUMBA = False
    prange = range                      # so @pnjit kernel bodies import/run under pure Python

# QUIVERLAB_NO_NUMBA=1 forces the pure-Python path (used by the parity test and as an escape
# hatch on hosts where numba misbehaves).
_DISABLED = os.environ.get("QUIVERLAB_NO_NUMBA", "") not in ("", "0", "false", "False")
USE_KERNELS = HAS_NUMBA and not _DISABLED


def njit(func):
    """@njit shim: real numba (cache=True) when available, else identity.

    Kernels decorated with @njit still *run* under pure Python when numba is absent
    (they use only numpy + loops), but the consumers won't call them in that case —
    they fall back to the pure-Python twins via USE_KERNELS.
    """
    if HAS_NUMBA:
        return _njit(cache=True)(func)
    return func


def pnjit(func):
    """@njit(parallel=True, cache=True) shim for prange kernels; identity if numba absent.

    Same contract as njit: when numba is absent the kernel is never called (USE_KERNELS is
    False -> consumers use the pure-Python twins), but the function must still import, so
    `prange` degrades to `range` in the import fallback above.
    """
    if HAS_NUMBA:
        return _njit(cache=True, parallel=True)(func)
    return func


@njit
def pow_mod(a, e, p):
    a = a % p
    r = 1 % p
    while e > 0:
        if e & 1:
            r = (r * a) % p
        a = (a * a) % p
        e >>= 1
    return r


@njit
def inv_mod(a, p):
    # p prime => a^(p-2) = a^{-1} (Fermat). a must be nonzero mod p.
    return pow_mod(a % p, p - 2, p)


@njit
def rank_mod_p_kernel(M, p):
    rows = M.shape[0]
    cols = M.shape[1]
    if rows == 0 or cols == 0:
        return 0
    A = np.empty((rows, cols), dtype=np.int64)
    for i in range(rows):
        for j in range(cols):
            A[i, j] = M[i, j] % p
    r = 0
    for c in range(cols):
        piv = -1
        for i in range(r, rows):
            if A[i, c] != 0:
                piv = i
                break
        if piv == -1:
            continue
        if piv != r:
            for j in range(cols):
                tmp = A[r, j]; A[r, j] = A[piv, j]; A[piv, j] = tmp
        inv = inv_mod(A[r, c], p)
        for j in range(cols):
            A[r, j] = (A[r, j] * inv) % p
        for i in range(rows):
            if i != r and A[i, c] != 0:
                f = A[i, c]
                for j in range(cols):
                    A[i, j] = (A[i, j] - f * A[r, j]) % p
        r += 1
        if r == rows:
            break
    return r


@njit
def nullspace_kernel(M, p):
    rows = M.shape[0]
    cols = M.shape[1]
    A = np.empty((rows, cols), dtype=np.int64)
    for i in range(rows):
        for j in range(cols):
            A[i, j] = M[i, j] % p
    where = -np.ones(cols, dtype=np.int64)     # pivot col -> its row in A
    pivots = np.empty(cols, dtype=np.int64)
    npiv = 0
    r = 0
    for c in range(cols):
        piv = -1
        for i in range(r, rows):
            if A[i, c] != 0:
                piv = i
                break
        if piv < 0:
            continue
        if piv != r:
            for j in range(cols):
                tmp = A[r, j]; A[r, j] = A[piv, j]; A[piv, j] = tmp
        inv = inv_mod(A[r, c], p)
        for j in range(cols):
            A[r, j] = (A[r, j] * inv) % p
        for i in range(rows):
            if i != r and A[i, c] != 0:
                f = A[i, c]
                for j in range(cols):
                    A[i, j] = (A[i, j] - f * A[r, j]) % p
        where[c] = r
        pivots[npiv] = c
        npiv += 1
        r += 1
        if r == rows:
            break
    is_piv = np.zeros(cols, dtype=np.int64)
    for k in range(npiv):
        is_piv[pivots[k]] = 1
    nfree = 0
    for c in range(cols):
        if is_piv[c] == 0:
            nfree += 1
    basis = np.zeros((nfree, cols), dtype=np.int64)
    row = 0
    for fcol in range(cols):
        if is_piv[fcol] == 1:
            continue
        basis[row, fcol] = 1
        for k in range(npiv):
            c = pivots[k]
            basis[row, c] = (-A[where[c], fcol]) % p
        row += 1
    return basis


@njit
def independent_modulo_kernel(span, cand, p):
    nspan = span.shape[0]
    L = span.shape[1]
    ncand = cand.shape[0]
    ech = np.zeros((L, L), dtype=np.int64)     # rank <= L pivot rows
    ech_piv = np.empty(L, dtype=np.int64)
    nech = 0
    work = np.empty(L, dtype=np.int64)
    # reduce span into an echelon basis
    for si in range(nspan):
        for j in range(L):
            work[j] = span[si, j] % p
        for e in range(nech):
            pc = ech_piv[e]
            if work[pc] != 0:
                f = work[pc]
                for j in range(L):
                    work[j] = (work[j] - f * ech[e, j]) % p
        lead = -1
        for j in range(L):
            if work[j] != 0:
                lead = j
                break
        if lead >= 0 and nech < L:
            inv = inv_mod(work[lead], p)
            for j in range(L):
                ech[nech, j] = (work[j] * inv) % p
            ech_piv[nech] = lead
            nech += 1
    chosen = np.empty(ncand, dtype=np.int64)
    nchosen = 0
    for ci in range(ncand):
        for j in range(L):
            work[j] = cand[ci, j] % p
        for e in range(nech):
            pc = ech_piv[e]
            if work[pc] != 0:
                f = work[pc]
                for j in range(L):
                    work[j] = (work[j] - f * ech[e, j]) % p
        lead = -1
        for j in range(L):
            if work[j] != 0:
                lead = j
                break
        if lead >= 0 and nech < L:
            inv = inv_mod(work[lead], p)
            for j in range(L):
                ech[nech, j] = (work[j] * inv) % p
            ech_piv[nech] = lead
            nech += 1
            chosen[nchosen] = ci
            nchosen += 1
    return chosen[:nchosen]


@pnjit
def radK_kernel(Lb, rad_ab, ker, cur_r, m2, p):
    # rad(A^e).ker: output row (ki, ri) = rad_ab[ri] applied to ker[ki].  Rows are
    # independent -- row = ki*nrad + ri writes a disjoint slice of `out` -- so the outer
    # loop runs in parallel (prange).  This is the measured deep-build hotspot (~77%).
    nnull = ker.shape[0]
    nrad = rad_ab.shape[0]
    L = m2 * cur_r
    out = np.zeros((nnull * nrad, L), dtype=np.int64)
    for ki in prange(nnull):
        for ri in range(nrad):
            a = rad_ab[ri, 0]
            b = rad_ab[ri, 1]
            row = ki * nrad + ri            # position-computed (no running counter -> prange-safe)
            for blk in range(cur_r):
                base = blk * m2
                for ii in range(m2):
                    # accumulate m2 products, reduce once per entry: safe while
                    # m2*(p-1)**2 < 2**63 (see module OVERFLOW note)
                    s = 0
                    for jj in range(m2):
                        c = Lb[a, b, ii, jj]
                        if c != 0:
                            s += c * ker[ki, base + jj]
                    out[row, base + ii] = s % p
    return out


@pnjit
def Dn_kernel(Lb, gens, cur_r, m, m2, p):
    # each generator j fills a disjoint column block (col = j*m2 + a*m + b) -> prange over j
    r_n = gens.shape[0]
    Dn = np.zeros((m2 * cur_r, m2 * r_n), dtype=np.int64)
    for j in prange(r_n):
        for a in range(m):
            for b in range(m):
                col = j * m2 + (a * m + b)
                for blk in range(cur_r):
                    base = blk * m2
                    for ii in range(m2):
                        # accumulate m2 products, reduce once per entry: safe while
                        # m2*(p-1)**2 < 2**63 (see module OVERFLOW note)
                        s = 0
                        for jj in range(m2):
                            c = Lb[a, b, ii, jj]
                            if c != 0:
                                s += c * gens[j, base + jj]
                        Dn[base + ii, col] = s % p
    return Dn
