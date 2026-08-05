"""Orientation-blind Dynkin / Euclidean type detection (Plan 38 / C2).

A quiver is classified by its UNDERLYING undirected multigraph (arrow
directions are irrelevant to the diagram type). Returns one of

    ("A", n) | ("D", n) | ("E", 6|7|8)          -- finite ADE
    ("~A", n) | ("~D", n) | ("~E", 6|7|8)       -- Euclidean (affine)
    None                                        -- anything else

Conventions match families/dynkin.py: A_n / D_n / E_n have n vertices;
~A_n is a cycle on n+1 vertices (~A_1 = the double edge / Kronecker),
~D_n / ~E_n have n+1 vertices. Loops, multi-edges (other than ~A_1),
disconnected graphs, and non-diagram shapes all give None.

Pure combinatorics; no field, no floats."""
from collections import deque


def _undirected(quiver):
    """Underlying multigraph: (vertices, adjacency v -> {w: multiplicity},
    has_loop). Loops are flagged; the caller returns None on any loop."""
    verts = list(quiver.vertices)
    adj = {v: {} for v in verts}
    has_loop = False
    for (s, t) in quiver.arrows.values():
        if s == t:
            has_loop = True
            continue
        adj[s][t] = adj[s].get(t, 0) + 1
        adj[t][s] = adj[t].get(s, 0) + 1
    return verts, adj, has_loop


def is_connected(quiver) -> bool:
    """True iff the underlying undirected graph is connected (BFS). An empty
    quiver (no vertices) is vacuously not connected."""
    verts, adj, _ = _undirected(quiver)
    if not verts:
        return False
    seen = {verts[0]}
    q = deque([verts[0]])
    while q:
        v = q.popleft()
        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                q.append(w)
    return len(seen) == len(verts)


def _degrees(verts, adj):
    return {v: sum(adj[v].values()) for v in verts}


def _arm_lengths(center, adj):
    """Lengths (edge counts) of the arms hanging off `center` in a tree whose
    only branch vertex is `center` (all other vertices have degree <= 2). Sorted
    ascending."""
    arms = []
    for start in adj[center]:
        length, prev, cur = 1, center, start
        while True:
            nbrs = [w for w in adj[cur] if w != prev]
            if len(nbrs) != 1:              # leaf (0) or branch (>=2): arm ends
                break
            prev, cur = cur, nbrs[0]
            length += 1
        arms.append(length)
    return tuple(sorted(arms))


def _classify_tree(verts, adj, deg):
    v = len(verts)
    degvals = sorted(deg.values())
    if degvals[-1] <= 2:                    # a path (or single vertex)
        return ("A", v)
    if degvals[-1] >= 5:
        return None
    deg3 = [u for u in verts if deg[u] == 3]
    deg4 = [u for u in verts if deg[u] == 4]

    # ~D4: a single degree-4 vertex with four leaves.
    if len(deg4) == 1 and not deg3:
        c = deg4[0]
        if all(deg[w] == 1 for w in adj[c]):
            return ("~D", 4)
        return None
    if deg4:                                # any other degree-4 config: not a diagram
        return None

    if len(deg3) == 1:                      # single fork: finite D / E or affine ~E
        arms = _arm_lengths(deg3[0], adj)
        a, b, c = arms
        if a == 1 and b == 1:
            return ("D", v)
        return {(1, 2, 2): ("E", 6), (1, 2, 3): ("E", 7), (1, 2, 4): ("E", 8),
                (2, 2, 2): ("~E", 6), (1, 3, 3): ("~E", 7),
                (1, 2, 5): ("~E", 8)}.get(arms)

    if len(deg3) == 2:                      # two forks: ~D_n (each with two leaves)
        for f in deg3:
            if sum(1 for w in adj[f] if deg[w] == 1) != 2:
                return None
        return ("~D", v - 1)

    return None


def dynkin_type(quiver):
    """The diagram type of `quiver` (see module docstring) or None."""
    verts, adj, has_loop = _undirected(quiver)
    if has_loop or not verts:
        return None
    if not is_connected(quiver):
        return None
    v = len(verts)
    e = sum(1 for (s, t) in quiver.arrows.values() if s != t)  # edges (no loops)
    deg = _degrees(verts, adj)
    if e == v - 1:                          # connected + E = V-1 => a tree
        return _classify_tree(verts, adj, deg)
    if e == v:                              # connected + E = V => a single cycle
        if all(d == 2 for d in deg.values()) and v >= 2:
            return ("~A", v - 1)            # ~A_1 = double edge, ~A_{>=2} = simple cycle
        return None
    return None
