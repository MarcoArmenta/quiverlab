"""Quiver with potential of an ideal triangulation (Plan 48). For an UNPUNCTURED surface
with boundary: quiver_of(T) has one vertex per arc and one arrow per anticlockwise angle
between consecutive arcs in a triangle (then 2-cycle reduced); potential_of(T) is the sum
of 3-cycles over internal triangles; jacobian_of(T) is P44's JacobianAlgebra. By ABCP 2010
/ Labardini 2009 the Jacobian is a GENTLE algebra (self-certified via P46's is_gentle).
Punctures / closed surfaces / self-folded triangles are OUT of v1 scope (loud). Float-free.

ORIENTATION (arbitrated by the disc oracle, not assumed). For each anticlockwise triangle
(s0, s1, s2) and each cyclically-consecutive pair (a, b) = (s_i, s_{i+1}) of sides that are
BOTH arcs, we add the arrow ``b -> a`` (the later ccw side to the earlier). The disc arbiter
``test_disc_fan_is_linear_A_n`` FIXES this direction: the fan of the (n+3)-gon then gives
the linear A_n quiver 1 -> 2 -> ... -> n. The naive ``a -> b`` reading gives the REVERSED
chain n -> ... -> 1 (documented deviation, flipped once per the house oracle-decides rule)."""
from __future__ import annotations

from collections import Counter

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.families.jacobian import JacobianAlgebra, Potential


def _require_v1_scope(T):
    S = T.surface
    if S.punctures > 0:
        raise QuiverlabError(
            "surfaces: punctured surfaces are out of v1 scope",
            hint="v1 handles unpunctured surfaces with boundary; puncture potentials "
                 "and self-folded triangles are the successor (P48.1)")
    if S.num_boundary_components == 0:
        raise QuiverlabError(
            "surfaces: closed surfaces are out of v1 scope",
            hint="v1 requires non-empty boundary (successor P48.1)")


def _is_arc(side):
    return isinstance(side, int) and not isinstance(side, bool)


def _raw_arrows(T):
    """The angle arrows before 2-cycle reduction: for each anticlockwise triangle, each
    cyclically-consecutive ordered pair of ARC sides (a, b) yields the arrow ``b -> a``
    (the disc-arbitrated orientation -- see the module docstring)."""
    edges = []                                            # list of (src_arc, tgt_arc)
    for tri in T.triangles:
        for i in range(3):
            a, b = tri[i], tri[(i + 1) % 3]
            if _is_arc(a) and _is_arc(b):
                edges.append((b, a))                      # arbitrated: later ccw -> earlier
    return edges


def _reduce_two_cycles(edges):
    """Cancel oriented 2-cycles pairwise: while both (a,b) and (b,a) occur, drop one of
    each. Returns the surviving multiset of directed edges (FST reduced quiver). For an
    in v1's unpunctured-with-boundary scope 2-cycles provably cannot arise (they need self-folded triangles) -- this is a defensive no-op kept for the P48.1 punctured successor; C(1,1) yields PARALLEL arrows, not a 2-cycle."""
    cnt = Counter(edges)
    for (a, b) in list(cnt):
        if a < b:                                         # handle each unordered pair once
            k = min(cnt[(a, b)], cnt[(b, a)])
            cnt[(a, b)] -= k
            cnt[(b, a)] -= k
    return [e for e, c in cnt.items() for _ in range(c)]


def quiver_of(T) -> Quiver:
    """One vertex per interior arc, one arrow per anticlockwise arc-to-arc angle (2-cycle
    reduced). Loud on out-of-v1-scope T (punctures / closed / self-folded)."""
    _require_v1_scope(T)
    arcs = list(T.arcs())
    edges = _reduce_two_cycles(_raw_arrows(T))
    arrows, k = {}, 0
    for (a, b) in sorted(edges):
        k += 1
        arrows[f"m{k}"] = (a, b)
    return Quiver(arcs, arrows)


def _internal_triangles(T):
    return [tri for tri in T.triangles if all(_is_arc(s) for s in tri)]


def _cycle_word(Q, arcs):
    """Order the arrows of Q whose endpoints all lie in the arc set ``arcs`` (an internal
    triangle) into a composable cyclic word. The internal triangle contributes exactly a
    3-cycle among its arcs; we walk source -> target to build the word."""
    out = {}                                              # src -> (arrow name, tgt)
    for a in Q.arrows:
        s, t = Q.arrows[a]
        if s in arcs and t in arcs:
            if s in out:
                raise QuiverlabError(
                    "potential_of: internal triangle is not a simple cycle",
                    hint="a degenerate (self-folded / double-arc) internal triangle is "
                         "out of v1 scope")
            out[s] = (a, t)
    if len(out) != len(arcs):
        raise QuiverlabError(
            "potential_of: internal-triangle arrows do not form a full cycle",
            hint="expected one arrow per arc of the internal triangle")
    start = next(iter(out))
    word, cur = [], start
    for _ in range(len(out)):
        name, nxt = out[cur]
        word.append(name)
        cur = nxt
    if cur != start:
        raise QuiverlabError("potential_of: internal-triangle walk did not close",
                             hint="internal triangle must be an oriented cycle")
    return tuple(word)


def potential_of(T) -> Potential:
    """A P44 Potential over quiver_of(T): one +1 3-cycle per internal triangle (a triangle
    all three of whose sides are arcs). For the fan (no internal triangle) W = 0."""
    Q = quiver_of(T)
    terms = []
    for tri in _internal_triangles(T):
        terms.append((1, _cycle_word(Q, set(tri))))
    return Potential(Q, terms)


def jacobian_of(T, field=None, degree_bound=None):
    """JacobianAlgebra(quiver_of(T), potential_of(T), field). NotFiniteDimensionalError
    propagates from P44 (it should not fire in v1 scope -- gentle => finite -- but stays
    honest). The potential is built on its own quiver object, which P44 requires to be the
    SAME object handed to JacobianAlgebra, so we read Q back off the potential."""
    W = potential_of(T)
    return JacobianAlgebra(W.quiver, W, field=field, degree_bound=degree_bound)
