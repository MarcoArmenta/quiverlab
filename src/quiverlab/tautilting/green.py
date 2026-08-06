"""Maximal green sequences (Plan 45 / C4). A maximal green sequence is a maximal chain of
GREEN (downward = left-mutation) edges from the top ``(A, 0)`` to the bottom ``(0, A)`` in
the oriented exchange (Hasse) quiver -- a source-to-sink monotone path. The count is a
subtle invariant (kA2 has exactly 2, the two sides of the pentagon). Honest cap: the DFS
over directed maximal chains stops LOUDLY at ``cap`` (``status="budget"``); an infinite
exchange graph short-circuits to ``status="budget"`` (unbounded / infinitely many MGSs)."""
from __future__ import annotations


def maximal_green_sequences(A, cap=4096):
    """All maximal green sequences of ``A`` as ``{"sequences": [[pair_id, ...], ...],
    "count": int, "complete": bool, "status": "complete"|"budget"}``. Each sequence is a
    directed maximal chain ``(A,0) -> ... -> (0,A)`` of downward (left-mutation) edges."""
    from quiverlab.tautilting.mutation import exchange_graph
    from quiverlab.tautilting.torsion import hasse_orientation
    eg = exchange_graph(A, budget_pairs=cap)
    if not eg.is_complete:                             # infinite exchange graph => unbounded
        return {"sequences": [], "count": 0, "complete": False, "status": "budget"}
    orient = hasse_orientation(eg)
    nverts = len(eg.vertices)
    down = {i: [] for i in range(nverts)}
    indeg = {i: 0 for i in range(nverts)}
    for (i, j), d in orient.items():
        a, b = (i, j) if d == "down" else (j, i)      # a -> b downward
        down[a].append(b)
        indeg[b] += 1
    source = next(i for i in range(nverts) if indeg[i] == 0)
    sink = next(i for i in range(nverts) if not down[i])
    sequences = []
    overflow = [False]

    def dfs(node, path):
        if overflow[0]:
            return
        if node == sink:
            if len(sequences) >= cap:
                overflow[0] = True
                return
            sequences.append(list(path))
            return
        for nxt in down[node]:
            dfs(nxt, path + [nxt])

    dfs(source, [source])
    if overflow[0]:
        return {"sequences": [], "count": 0, "complete": False, "status": "budget"}
    return {"sequences": sequences, "count": len(sequences),
            "complete": True, "status": "complete"}
