# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
"""
Search III -- the homology/cohomology ASYMMETRY.

Searches I and II computed Hochschild HOMOLOGY only.  Han's conjecture is about homology,
and the reason it is hard (and the reason Happel's original question was FALSE) is that
Hochschild COHOMOLOGY HH^n can vanish in high degrees while HH_n does not.  This module adds a
cohomology engine and exhibits the asymmetry on the same algebras.

Cohomology cochain complex (normalised):
    C^n = Hom_k( \bar A^{(x)n}, A),   coboundary  d: C^n -> C^{n+1},
    (d f)(a_1,...,a_{n+1}) = a_1 f(a_2,...,a_{n+1})
        + sum_{i=1}^n (-1)^i f(a_1,...,a_i a_{i+1},...,a_{n+1})
        + (-1)^{n+1} f(a_1,...,a_n) a_{n+1},
with a_i in the reduced basis \bar A = A/k1 (interior products projected by dropping the unit
coordinate, exactly as in the homology engine).  This coboundary is NOT the transpose of the
homology boundary; the two complexes agree only when A is symmetric (A = A^* as bimodules).

dim HH^n = dim C^n - rank d^n - rank d^{n-1}.
"""

import json
import sys
import numpy as np
import itertools
from quiverlab.engine.hh_engine import (Algebra, hochschild_homology_dims, check_associative, rank_mod_p,
                       truncated_polynomial, two_gen_local)
from quiverlab.engine.linalg_fast import rank_mod_p_auto


def cochain_basis(alg, n):
    """basis of C^n = Hom(\bar A^{(x)n}, A): pairs (word, j), word in R^n, j in 0..m-1."""
    R = alg.R
    words = [()] if n == 0 else list(itertools.product(R, repeat=n))
    return [(w, j) for w in words for j in range(alg.m)]


def coboundary_matrix(alg, n, basis_n, index_np1):
    """Matrix of d^n: C^n -> C^{n+1}. shape (dim C^{n+1}, dim C^n)."""
    m = alg.m
    t = alg.t
    R = alg.R
    rows = len(index_np1)
    cols = len(basis_n)
    M = np.zeros((rows, cols), dtype=np.int64)
    # enumerate input words u in R^{n+1}
    u_words = list(itertools.product(R, repeat=n + 1))
    ej = np.eye(m, dtype=np.int64)
    for col, (w, j) in enumerate(basis_n):
        for u in u_words:
            out = np.zeros(m, dtype=np.int64)
            # T0: a_1 * f(a_2..a_{n+1});  nonzero iff u[1:] == w
            if u[1:] == w:
                out += alg.mult_full(u[0], j)
            # interior terms i = 1..n  (contract u[i-1], u[i])
            for i in range(1, n + 1):
                P = alg.mult_full(u[i - 1], u[i])
                sign = -1 if (i % 2 == 1) else 1
                nz = np.nonzero(P)[0]
                for c in nz:
                    if c == t:
                        continue
                    merged = u[:i - 1] + (int(c),) + u[i + 1:]
                    if merged == w:
                        out += sign * P[c] * ej[j]
            # T_{n+1}: f(a_1..a_n) * a_{n+1};  nonzero iff u[:n] == w
            if u[:n] == w:
                sign = -1 if ((n + 1) % 2 == 1) else 1
                out += sign * alg.mult_full(j, u[n])
            nzp = np.nonzero(out)[0]
            for p in nzp:
                M[index_np1[(u, int(p))], col] += out[p]
    return M


def hochschild_cohomology_dims(alg, N, primes=(32003, 2, 3), resolution=None):
    """Return dict p -> [dim HH^0,...,dim HH^N].

    resolution: a Resolution backend supplying cochain_basis / coboundary_matrix
    (default: the dual normalized bar complex).
    """
    from quiverlab.engine.resolutions import _default
    resolution = _default(resolution)
    bases = {n: resolution.cochain_basis(alg, n) for n in range(0, N + 2)}
    index = {n: {g: i for i, g in enumerate(bases[n])} for n in range(0, N + 2)}
    dmats = {}
    for n in range(0, N + 1):
        dmats[n] = resolution.coboundary_matrix(alg, n, bases[n], index[n + 1])
    ranks = {p: {} for p in primes}
    for p in primes:
        for n in range(0, N + 1):
            ranks[p][n] = rank_mod_p_auto(dmats[n], p)
    out = {}
    for p in primes:
        dims = []
        for n in range(0, N + 1):
            dimCn = len(bases[n])
            rn = ranks[p].get(n, 0)          # rank d^n  (out of C^n)
            rnm1 = ranks[p].get(n - 1, 0)    # rank d^{n-1} (into C^n)
            dims.append(dimCn - rn - rnm1)
        out[p] = dims
    return out


def complexity_of(seq_full):
    """apparent complexity of a (co)homology sequence: 1 + poly growth degree; 0 if eventually 0."""
    seq = seq_full[1:]
    if len(seq) == 0:
        return None
    if all(d == 0 for d in seq):
        return 0
    if len(seq) < 3:
        return None
    # eventual vanishing: trailing zeros (cohomology can die in high degree)
    if seq[-1] == 0 and seq[-2] == 0:
        return 0
    last = seq[-4:] if len(seq) >= 4 else seq[-3:]
    if len(set(last)) == 1:
        return 0 if last[0] == 0 else 1
    cur = np.array(seq, dtype=np.int64)
    for k in range(1, 5):
        cur = np.diff(cur)
        if cur.size >= 2 and np.all(cur == 0):
            return k
        if cur.size >= 2 and np.all(cur == cur[0]) and cur[0] != 0:
            return k + 1
    return ">=2"


def quantum_ci(c):
    """k<x,y>/(x^2, y^2, yx - c*xy). dim 4 quantum complete intersection."""
    return two_gen_local([0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, c],
                         name=f"qCI(yx={c}xy)")
