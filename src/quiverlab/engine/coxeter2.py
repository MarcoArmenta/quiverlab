# Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
# github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
# Mechanical changes only: package-relative imports, __main__ blocks removed,
# float literals eradicated (quiverlab AST gate), env guard renamed.
"""
Search V -- ELEMENTARY (multi-vertex) algebras: testing Thm B, and the AR/Coxeter action.

Armenta's Thm B is stated for ELEMENTARY algebras: the automorphism Theta of the TT calculus
induced by the Auslander-Reiten / Serre bimodule acts on Hochschild homology with characteristic
polynomial equal to the COXETER POLYNOMIAL.  A local algebra gives the degenerate Coxeter polynomial
t+1; the content is multi-vertex.

This module:
  (1) builds multi-vertex elementary algebras:
        - self-injective Nakayama algebras  kZ_n/rad^ell  (Theta = Nakayama automorphism = a rotation),
        - hereditary path algebras of Dynkin quivers  kA_n, kD_4  (finite gl.dim; Theta = Serre functor),
  (2) computes the Cartan matrix and, where invertible, the classical Coxeter matrix
        Phi = -C^{-T} C  and its Coxeter polynomial,
  (3) for the self-injective case, computes the induced Theta-action on HH_n and HH^n via the
        Search IV engine (Theta is an honest algebra automorphism there) and reports
        char poly of Theta on HH_* (test of Thm B) and Theta = id on HH^* (Thm C).
"""

import json
import sys
import numpy as np
import sympy as sp
from quiverlab.engine.hh_engine import Algebra, hochschild_homology_dims
from quiverlab.engine.scan3 import hochschild_cohomology_dims
from quiverlab.engine.coxeter import (induced_on_HH_homology, induced_on_HH_cohomology, is_identity)


# ====================================================================
# builders for multi-vertex elementary algebras
# ====================================================================
def cyclic_nakayama(n, ell):
    """kZ_n / rad^ell : cyclic quiver on n vertices, paths of length < ell. dim = n*ell. self-injective."""
    m = n * ell
    def idx(i, L):
        return (i % n) * ell + L
    T = np.zeros((m, m, m), dtype=np.int64)
    for i in range(n):
        for s in range(ell):
            for j in range(n):
                for tt in range(ell):
                    if (i + s) % n == j and s + tt < ell:
                        T[idx(i, s), idx(j, tt), idx(i, s + tt)] = 1
    unit = np.zeros(m, dtype=np.int64)
    for i in range(n):
        unit[idx(i, 0)] = 1
    return Algebra(m, T, unit, name=f"kZ_{n}/rad^{ell}"), idx


def rotation_matrix_orig(n, ell, k):
    """rotation by k steps on kZ_n/rad^ell, in the ORIGINAL basis p(i,L) -> p(i+k, L)."""
    m = n * ell
    def idx(i, L):
        return (i % n) * ell + L
    S = np.zeros((m, m), dtype=np.int64)
    for i in range(n):
        for L in range(ell):
            S[idx(i + k, L), idx(i, L)] = 1
    return S


def linear_path_algebra(n):
    """kA_n : linear quiver 0->1->...->(n-1), all paths. hereditary, finite gl.dim. dim = n(n+1)/2."""
    # basis: paths p(i,j) = i->...->j for 0<=i<=j<=n-1
    pairs = [(i, j) for i in range(n) for j in range(i, n)]
    index = {pr: k for k, pr in enumerate(pairs)}
    m = len(pairs)
    T = np.zeros((m, m, m), dtype=np.int64)
    for (i, j) in pairs:
        for (a, b) in pairs:
            if j == a:  # compose p(i,j) then p(a,b)=p(j,b)
                T[index[(i, j)], index[(a, b)], index[(i, b)]] = 1
    unit = np.zeros(m, dtype=np.int64)
    for i in range(n):
        unit[index[(i, i)]] = 1
    return Algebra(m, T, unit, name=f"kA_{n}"), index, pairs


def quiver_path_algebra(n, arrows, name="kQ"):
    """Path algebra kQ of a finite acyclic quiver Q (no oriented cycles).

    n       -- number of vertices (0..n-1);
    arrows  -- list of directed edges (s, e), s != e (Q must be acyclic).

    Returns (Algebra, paths) where `paths` is the ordered list of basis paths, each a
    tuple of visited vertices (v0, ..., vk) (trivial path (i,) = e_i).  Composition is
    concatenation left-to-right: p(i..j) * q(j..b) = p(i..b) (matching
    `linear_path_algebra`).  Assumes Q is a quiver of a tree / simply-laced Dynkin
    diagram (at most one arrow between two vertices), so a path is determined by its
    vertex sequence.
    """
    out_arrows = {i: [] for i in range(n)}
    for (s, e) in arrows:
        out_arrows[s].append(e)
    # enumerate all directed paths by DFS (acyclic => finite)
    paths = []
    def dfs(seq):
        paths.append(tuple(seq))
        for nxt in out_arrows[seq[-1]]:
            dfs(seq + [nxt])
    for i in range(n):
        dfs([i])
    index = {p: k for k, p in enumerate(paths)}
    m = len(paths)
    T = np.zeros((m, m, m), dtype=np.int64)
    for p in paths:
        for q in paths:
            if p[-1] == q[0]:                 # p ends where q starts
                comp = p + q[1:]
                if comp in index:             # comp is a genuine path of Q
                    T[index[p], index[q], index[comp]] = 1
    unit = np.zeros(m, dtype=np.int64)
    for i in range(n):
        unit[index[(i,)]] = 1
    return Algebra(m, T, unit, name=name), paths


# ====================================================================
# Dynkin quivers (simply-laced): A_n, D_n, E_6, E_7, E_8
# ====================================================================
# Coxeter number h and the known classical Coxeter polynomial (orientation-
# independent) are the oracle: the Coxeter transformation Phi = -C^{-T}C has order
# exactly h, its characteristic polynomial is palindromic with all roots on the unit
# circle, and equals the tabulated product of cyclotomic polynomials.
def dynkin_quiver(typ, n=None):
    """Return (n_vertices, arrows, coxeter_number) for a simply-laced Dynkin diagram.

    `typ` in {'A','D','E'}; for 'A'/'D' pass the rank `n`; for 'E' pass n in {6,7,8}.
    A linear/standard orientation is chosen (the Coxeter polynomial is orientation-
    independent, so any acyclic orientation gives the same answer).
    """
    if typ == "A":
        arrows = [(i, i + 1) for i in range(n - 1)]
        return n, arrows, n + 1
    if typ == "D":
        # chain 0->1->...->(n-3) with two extra sinks (n-2),(n-1) off vertex (n-3)
        assert n >= 4
        arrows = [(i, i + 1) for i in range(n - 3)]
        arrows += [(n - 3, n - 2), (n - 3, n - 1)]
        return n, arrows, 2 * (n - 1)
    if typ == "E":
        assert n in (6, 7, 8)
        # chain 0-1-2-3-...-(n-2) with a branch vertex (n-1) attached to vertex 2
        arrows = [(i, i + 1) for i in range(n - 2)]
        arrows += [(2, n - 1)]
        h = {6: 12, 7: 18, 8: 30}[n]
        return n, arrows, h
    raise ValueError(f"unknown Dynkin type {typ!r}")


def coxeter_element_order(C):
    """Order of the Coxeter transformation Phi = -C^{-T}C (smallest k with Phi^k = I),
    or None if Phi is not of finite order.  Exact over the rationals via sympy."""
    Cs = sp.Matrix(C.tolist())
    if Cs.det() == 0:
        return None
    Phi = -(Cs.inv().T) * Cs
    I = sp.eye(C.shape[0])
    M = I.copy()
    for k in range(1, 1000):
        M = M * Phi
        if M == I:
            return k
    return None


# ====================================================================
# Cartan matrix and classical Coxeter polynomial
# ====================================================================
# The public `cartan_matrix` (hanlab.cartan_matrix) is an alias for `cartan_from_raw`
# below. A general "Cartan matrix from idempotent indices in the original basis"
# builder was drafted here but never needed; the dead NotImplementedError stub
# (which shadowed nothing but was a footgun for direct `from coxeter2 import
# cartan_matrix` callers) has been removed.


def cartan_from_raw(n, basis_vertex_of, paths):
    """Cartan matrix for a path algebra given (start,end) of each basis path.
       C[i,j] = number of basis paths from vertex j to vertex i  (mult of S_i in P_j = e_j A)."""
    C = np.zeros((n, n), dtype=np.int64)
    for (s, e) in paths:
        # path from s to e contributes to e_e A e_s? path p: s->e satisfies e_e p e_s = p (left mult by e_e, right by e_s)
        # P_j = e_j A has basis paths starting at j; composition factor S_i counts paths j->i
        # so mult of S_i in P_j = #{paths j -> i}
        C[e, s] += 1   # path s->e : a path starting at s ending at e -> in P_s, factor S_e
    return C


def coxeter_polynomial_from_cartan(C):
    """classical Coxeter polynomial det(t*I - Phi), Phi = -C^{-T} C, if det C != 0. Exact via sympy."""
    Cs = sp.Matrix(C.tolist())
    d = Cs.det()
    if d == 0:
        return None, None
    Phi = -(Cs.inv().T) * Cs
    t = sp.symbols('t')
    poly = (t * sp.eye(C.shape[0]) - Phi).det()
    return sp.factor(sp.expand(poly)), Phi


# ====================================================================
# char poly of an F_p matrix (lift to balanced integers; entries are small for permutation actions)
# ====================================================================
def balanced_lift(M, p):
    M = np.asarray(M, dtype=np.int64) % p
    return ((M + p // 2) % p) - p // 2


def charpoly_of_induced(M, p):
    if M.shape[0] == 0:
        return sp.Integer(1)
    L = balanced_lift(M, p)
    t = sp.symbols('t')
    Ms = sp.Matrix(L.tolist())
    return sp.factor(sp.expand((t * sp.eye(M.shape[0]) - Ms).det()))


def auto_to_f_basis(alg, unit, S_orig):
    """transform an automorphism matrix from the original basis to the unit-adapted f-basis."""
    m = alg.m
    t = alg.t
    B = np.eye(m, dtype=np.int64)
    B[:, t] = np.asarray(unit, dtype=np.int64)
    Binv = np.round(np.linalg.inv(B.astype(float))).astype(np.int64)
    assert np.array_equal(B @ Binv, np.eye(m, dtype=np.int64)), (
        "auto_to_f_basis: change-of-basis B is not unimodular (need unit[t] == 1); "
        "the float-rounded integer inverse is unreliable here."
    )
    return (Binv @ S_orig @ B)
