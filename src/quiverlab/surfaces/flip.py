"""Quadrilateral flips + the flip<->mutation certificate (Plan 48). Flipping an interior
arc replaces it with the other diagonal of the quadrilateral formed by its two adjacent
triangles. The Fomin-Shapiro-Thurston theorem: quiver_of(flip(T, a)) equals the Fomin-
Zelevinsky matrix mutation of quiver_of(T) at vertex a. We certify that identity per
instance (exact skew-symmetric integer combinatorics); full DWZ right-equivalence of the
potentials is a named successor (P48.1), not attempted here. Float-free.

The quad flip reads only the anticlockwise SIDE tuples: for the two triangles sharing arc
``a``, T1 = (.., b, a, c, ..) contributes its sides ``b`` (before ``a``) and ``c`` (after);
T2 = (.., a, d, e, ..) contributes ``d`` (after ``a``) and ``e`` (next). The quadrilateral
is anticlockwise (b, d, e, c) with old diagonal ``a``; the new diagonal ``a'`` (reusing
``a``'s id, so the exchange matrices stay index-aligned) splits it into the anticlockwise
triangles (b, d, a') and (a', e, c)."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.surfaces.qp import quiver_of
from quiverlab.surfaces.triangulation import Triangulation, _is_arc


def flip(T, arc) -> Triangulation:
    """The combinatorial quadrilateral flip of interior arc ``arc``. Loud on a boundary
    segment (not flippable) or a non-flippable arc (does not border 2 distinct triangles --
    a self-folded configuration, out of v1 scope)."""
    if not _is_arc(arc):
        raise QuiverlabError(f"flip: {arc!r} is a boundary segment, not flippable",
                             hint="only interior arcs flip")
    idxs = [i for i, tri in enumerate(T.triangles) if arc in tri]
    if len(idxs) != 2:
        raise QuiverlabError(f"flip: arc {arc!r} does not border 2 distinct triangles",
                             hint="non-flippable (self-folded configuration, out of scope)")
    i1, i2 = idxs
    T1, T2 = T.triangles[i1], T.triangles[i2]
    p1, p2 = T1.index(arc), T2.index(arc)
    b = T1[(p1 - 1) % 3]                                  # side before `arc` in T1 (ccw)
    c = T1[(p1 + 1) % 3]                                  # side after `arc` in T1 (ccw)
    d = T2[(p2 + 1) % 3]                                  # side after `arc` in T2 (ccw)
    e = T2[(p2 + 2) % 3]                                  # next side in T2 (= before `arc`)
    new1 = (b, d, arc)                                    # anticlockwise, reuse arc's id
    new2 = (arc, e, c)
    kept = [tri for i, tri in enumerate(T.triangles) if i not in (i1, i2)]
    return Triangulation(T.surface, tuple(kept) + (new1, new2))


def _exchange_on_order(Q, order):
    """The skew-symmetric exchange matrix B (B[i][j] = #(i->j) - #(j->i)) with vertices in
    the supplied ``order``."""
    idx = {v: i for i, v in enumerate(order)}
    n = len(order)
    B = [[0] * n for _ in range(n)]
    for a in Q.arrows:
        i, j = idx[Q.source(a)], idx[Q.target(a)]
        B[i][j] += 1
        B[j][i] -= 1
    return B


def exchange_matrix(Q):
    """The skew-symmetric exchange matrix of Q with vertices in Q's own order."""
    return _exchange_on_order(Q, list(Q.vertices))


def matrix_mutation(B, k):
    """Fomin-Zelevinsky mu_k on a skew-symmetric integer matrix (Cluster algebras I,
    2002). Indices are VERTEX POSITIONS -- the caller aligns k to the arc's position."""
    n = len(B)
    Bp = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == k or j == k:
                Bp[i][j] = -B[i][j]
            else:
                Bp[i][j] = B[i][j] + (abs(B[i][k]) * B[k][j]
                                      + B[i][k] * abs(B[k][j])) // 2
    return Bp


def certify_flip_mutation(T, arc) -> bool:
    """True iff exchange_matrix(quiver_of(flip(T, arc))) equals the Fomin-Zelevinsky
    matrix mutation of exchange_matrix(quiver_of(T)) at vertex ``arc``. The per-instance
    flip<->mutation certificate (FST 2008) -- at the QUIVER level; DWZ potential
    right-equivalence is deferred (P48.1)."""
    Q = quiver_of(T)
    order = list(Q.vertices)
    k = order.index(arc)
    B_mut = matrix_mutation(_exchange_on_order(Q, order), k)
    Qf = quiver_of(flip(T, arc))
    # quiver_of sorts vertices by arc id and flip reuses arc's id, so both share the same
    # vertex ordering; align Qf onto `order` explicitly to be safe.
    B_flip = _exchange_on_order(Qf, order)
    return B_flip == B_mut
