"""Commutative ladders CL(n) = A_n [] A_2 (Plan 69 / R33; Escolar-Hiraoka 2016).

A **commutative ladder** is the box product of a type-A base quiver ``A_n`` (any
orientation ``tau``) with the rung quiver ``A_2`` (``1 -> 2``): two parallel rows
of ``n`` vertices joined by ``n`` rungs, with **every square commutative**. By
Escolar-Hiraoka (*Persistence Modules on Commutative Ladders of Finite Type*,
Discrete & Computational Geometry 55(1) 2016, 100-157) ``CL(n)`` is
**representation-finite iff n <= 4**; ``n >= 5`` is representation-infinite and
refused loudly (never a partial/silent diagram).

Representation theory FIRST: a *persistence module* on a commutative ladder is a
representation of this bound quiver; its *generalized persistence diagram* is the
Krull-Schmidt decomposition indexed by the AR quiver (see ``modules/barcode.py``).

VERTICES ARE SCALAR (H1, a hard build gate): the natural ``[n]x[2]`` grid labels
``(i, j)`` are TUPLES, which crash ``knit_ar_quiver`` at ``complex_reps.py:179``
(``"e_%s" % (i, j)`` -> ``TypeError``); ``CommutativeLadder`` therefore names its
vertices ``"i_j"``. The forward ``CL(n)`` is exactly the incidence algebra of the
scalar-relabelled ``[n]x[2]`` grid poset (``dim = 3*C(n+1, 2)``: 9 / 18 / 30 for
n = 2 / 3 / 4, live-verified).
"""
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError

_MAX_FINITE = 4          # Escolar-Hiraoka: CL(n) rep-finite iff n <= 4


def _orientation_list(n, base_orientation):
    """The base-edge orientation as a length-``n-1`` list of ``'f'``/``'b'``.

    ``base_orientation`` is ``"forward"`` (all ``'f'``), ``"backward"`` (all
    ``'b'``), or a ``dict`` keyed by the edge's lower vertex ``i`` in ``1..n-1``
    (or by the pair ``(i, i+1)``) mapping to ``'f'``/``'b'`` (default ``'f'``)."""
    m = n - 1
    if base_orientation == "forward":
        return ["f"] * m
    if base_orientation == "backward":
        return ["b"] * m
    if isinstance(base_orientation, dict):
        out = []
        for e in range(1, n):
            d = base_orientation.get(e, base_orientation.get((e, e + 1), "f"))
            if d not in ("f", "b"):
                raise QuiverlabError(
                    f"CommutativeLadder base_orientation[{e}] = {d!r} must be 'f' or 'b'",
                    hint="a base_orientation dict maps edge i (1..n-1) -> 'f' or 'b'")
            out.append(d)
        return out
    raise QuiverlabError(
        f"CommutativeLadder base_orientation must be 'forward', 'backward', or a dict, "
        f"got {base_orientation!r}",
        hint="e.g. base_orientation='forward' or {1: 'f', 2: 'b'}")


def CommutativeLadder(n, base_orientation="forward", field=None, expected_dim=None):
    """The commutative ladder ``CL(n) = A_n [] A_2`` as a presented ``kQ/I``.

    ``n`` is the base length; ``base_orientation in {"forward","backward", <dict>}``
    orients the base ``A_n`` edges (the rung direction ``(i,1) -> (i,2)`` is fixed).
    Every square commutes (the ``IncidenceAlgebra`` ``"a*b - c*d"`` idiom, left-to-
    right composition). **Refuses ``n >= 5`` LOUDLY** (Escolar-Hiraoka:
    representation-infinite) and ``n < 1``.

    Vertices are SCALAR (``"i_j"``) -- tuple labels crash the AR knit (H1). For the
    fully-forward orientation the dimension certificate ``3*C(n+1, 2)`` is asserted
    (9 / 18 / 30 for n = 2 / 3 / 4); a non-forward ``tau`` accepts the engine's
    ``dim`` (or the optional ``expected_dim=``, loud on mismatch)."""
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise QuiverlabError(
            f"CommutativeLadder needs an integer n >= 1, got {n!r}",
            hint="n is the base length; CL(n) = A_n [] A_2, rep-finite for n <= 4")
    if n > _MAX_FINITE:
        raise QuiverlabError(
            f"CommutativeLadder: CL({n}) is representation-infinite for n >= 5 "
            "(Escolar-Hiraoka 2016, DCG 55(1)); the AR-indexed generalized persistence "
            "diagram is undefined -- there is no complete AR quiver",
            hint="use n <= 4 (the representation-finite commutative ladders)")
    orient = _orientation_list(n, base_orientation)

    verts = [f"{i}_{j}" for i in range(1, n + 1) for j in (1, 2)]
    arrows = {}
    for i in range(1, n + 1):                        # rungs (i,1) -> (i,2)
        arrows[f"r{i}"] = (f"{i}_1", f"{i}_2")
    for e in range(1, n):                            # base arrows in both rows
        d = orient[e - 1]
        for j in (1, 2):
            arrows[f"h{e}_{j}"] = ((f"{e}_{j}", f"{e + 1}_{j}") if d == "f"
                                   else (f"{e + 1}_{j}", f"{e}_{j}"))
    rels = []                                        # one commutativity relation per square
    for e in range(1, n):
        d = orient[e - 1]
        h1, h2, rl, rr = f"h{e}_1", f"h{e}_2", f"r{e}", f"r{e + 1}"
        # forward square:  (e,1)->(e+1,1)->(e+1,2)  ==  (e,1)->(e,2)->(e+1,2)
        # backward square: (e+1,1)->(e,1)->(e,2)    ==  (e+1,1)->(e+1,2)->(e,2)
        rels.append(f"{h1}*{rr} - {rl}*{h2}" if d == "f"
                    else f"{h1}*{rl} - {rr}*{h2}")

    A = Quiver(verts, arrows).algebra(relations=rels, field=field)

    if all(d == "f" for d in orient):
        want = 3 * (n * (n + 1) // 2)                # 3 * C(n+1, 2)
    else:
        want = expected_dim
    if want is not None and A.dim != want:
        raise QuiverlabError(
            f"CommutativeLadder(n={n}, base_orientation={base_orientation!r}) built "
            f"dim {A.dim}, expected {want}",
            hint="pass expected_dim= for a non-forward orientation, or report a bug if "
                 "the forward dimension certificate 3*C(n+1,2) drifted")
    A._family_citations = ("escolar_hiraoka", "botnan_crawley_boevey", "assem_book")
    return A


def persistence_line(n, orientation="forward"):
    """The persistence-line base ``A_n`` as a :class:`Quiver` (no relations).

    ``orientation`` is ``"forward"`` (``1 -> 2 -> ... -> n``, ordinary
    persistence), ``"zigzag"`` (``1 -> 2 <- 3 -> 4 <- ...``, a zigzag module), or a
    ``dict`` keyed by the edge's lower vertex ``i`` mapping to ``'f'``/``'b'``."""
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise QuiverlabError(f"persistence_line needs an integer n >= 1, got {n!r}",
                             hint="n is the number of filtration steps")
    if orientation == "zigzag":
        orient = ["f" if e % 2 == 1 else "b" for e in range(1, n)]
    else:
        orient = _orientation_list(n, orientation)
    verts = list(range(1, n + 1))
    arrows = {}
    for e in range(1, n):
        arrows[f"a{e}"] = ((e, e + 1) if orient[e - 1] == "f" else (e + 1, e))
    return Quiver(verts, arrows)


def _undirected(Q):
    """``(vertices, adjacency-dict, simple-edge-set)`` of the underlying graph."""
    V = list(Q.vertices)
    adj = {v: set() for v in V}
    edges = set()
    for (s, t) in Q.arrows.values():
        if s != t:
            adj[s].add(t)
            adj[t].add(s)
            edges.add(frozenset((s, t)))
    return V, adj, edges


def _connected(V, adj):
    if not V:
        return False
    seen = {V[0]}
    stack = [V[0]]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == len(V)


def _is_ladder(V, adj, n):
    """Exact test that the undirected graph is the ladder ``P_n [] K_2`` -- a
    backtracking bijection ``V -> {1..n} x {1,2}`` preserving adjacency. With the
    already-checked equal edge count (``3n-2``) this is an isomorphism certificate
    (every G-edge maps to a ladder-edge; equal counts force onto)."""
    targets = [(i, j) for i in range(1, n + 1) for j in (1, 2)]

    def t_adj(a, b):
        (i, j), (k, l) = a, b
        return (abs(i - k) == 1 and j == l) or (i == k and abs(j - l) == 1)

    # place highest-degree vertices first to prune fast
    order = sorted(V, key=lambda v: -len(adj[v]))
    assign, used = {}, set()

    def bt(idx):
        if idx == len(order):
            return True
        v = order[idx]
        for t in targets:
            if t in used:
                continue
            if all(t_adj(t, assign[u]) for u in adj[v] if u in assign):
                assign[v] = t
                used.add(t)
                if bt(idx + 1):
                    return True
                del assign[v]
                used.discard(t)
        return False

    return bt(0)


def is_commutative_ladder(A):
    """Recognize ``A_n [] A_2`` by the underlying **graph shape** -> ``(True, n)``,
    else ``(False, 0)``. Refuses loudly on a presentation-less (structure-constant)
    algebra. Accepts a hand-built ``CL(>=5)`` (returns ``(True, n>=5)``) so
    ``barcode`` can then refuse it (rep-infinite)."""
    Q = getattr(A, "quiver", None)
    if Q is None:
        raise QuiverlabError(
            "is_commutative_ladder needs a quiver-presented algebra (kQ/I); the given "
            "algebra is presentation-less (structure constants)",
            hint="build via CommutativeLadder(n) or Quiver(...).algebra(...)")
    V, adj, edges = _undirected(Q)
    m = len(V)
    if m < 4 or m % 2 != 0:                          # ladder has 2n >= 4 vertices
        return (False, 0)
    n = m // 2
    if len(edges) != 3 * n - 2:                      # n rungs + 2(n-1) rails
        return (False, 0)
    degs = sorted(len(adj[v]) for v in V)
    if degs != sorted([2] * 4 + [3] * (2 * n - 4)):  # 4 corners deg 2, rest deg 3
        return (False, 0)
    if not _connected(V, adj):
        return (False, 0)
    if not _is_ladder(V, adj, n):
        return (False, 0)
    return (True, n)
