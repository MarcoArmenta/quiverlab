"""Layered layout for a quiver with relations, in EXACT int/Fraction coordinates
(spec §3.7). Vertices are placed in columns by their longest-path depth on the
strongly-connected-component condensation (so loops and cycles are handled --
all members of a cycle share a column); within a column they are centered on
integer/half-integer rows. Parallel arrows fan out with symmetric Fraction
bends; loops become LoopRoutes at integer base angles. NO floats: matplotlib
receives these exact numbers and coerces them itself (see viz/__init__.py)."""
from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class EdgeRoute:
    name: str
    src: object
    tgt: object
    kind: str          # "straight" | "parallel"
    bend: Fraction     # 0 for straight; symmetric offsets for parallel bundles


@dataclass(frozen=True)
class LoopRoute:
    name: str
    at: object
    angle_deg: int     # integer base angle of the self-arc


@dataclass(frozen=True)
class LayoutData:
    positions: dict     # vertex -> (x:int, y:Fraction)
    edges: tuple        # tuple[EdgeRoute, ...]  (non-loop arrows)
    loops: tuple        # tuple[LoopRoute, ...]
    relations: tuple    # tuple[str, ...]
    columns: tuple      # tuple[tuple[vertex, ...], ...] by ascending depth


def _sccs(quiver):
    """Strongly-connected components via iterative Tarjan; returns
    comp: dict[vertex,int] mapping each vertex to its component id."""
    index = {}
    low = {}
    onstack = {}
    stack = []
    comp = {}
    counter = [0]
    ncomp = [0]
    adj = {v: [] for v in quiver.vertices}
    for _n, (s, t) in quiver.arrows.items():
        adj[s].append(t)
    for root in quiver.vertices:
        if root in index:
            continue
        work = [(root, 0)]                 # (vertex, next-neighbor-index)
        while work:
            v, i = work[-1]
            if i == 0:
                index[v] = low[v] = counter[0]
                counter[0] += 1
                stack.append(v)
                onstack[v] = True
            recursed = False
            neigh = adj[v]
            while i < len(neigh):
                w = neigh[i]
                i += 1
                if w not in index:
                    work[-1] = (v, i)
                    work.append((w, 0))
                    recursed = True
                    break
                if onstack.get(w):
                    low[v] = min(low[v], index[w])
            if recursed:
                continue
            work[-1] = (v, i)
            if low[v] == index[v]:
                while True:
                    w = stack.pop()
                    onstack[w] = False
                    comp[w] = ncomp[0]
                    if w == v:
                        break
                ncomp[0] += 1
            work.pop()
            if work:
                p, _pi = work[-1]
                low[p] = min(low[p], low[v])
    return comp


def layer(quiver):
    """Longest-path depth per vertex on the SCC condensation (handles cycles/loops)."""
    comp = _sccs(quiver)
    ncomp = max(comp.values()) + 1 if comp else 0
    cadj = {c: set() for c in range(ncomp)}
    cindeg = {c: 0 for c in range(ncomp)}
    for _name, (s, t) in quiver.arrows.items():
        cs, ct = comp[s], comp[t]
        if cs != ct and ct not in cadj[cs]:
            cadj[cs].add(ct)
            cindeg[ct] += 1
    depth = {c: 0 for c in range(ncomp)}
    queue = [c for c in range(ncomp) if cindeg[c] == 0]
    order = []
    while queue:
        c = queue.pop()
        order.append(c)
        for d in cadj[c]:
            if depth[c] + 1 > depth[d]:
                depth[d] = depth[c] + 1
            cindeg[d] -= 1
            if cindeg[d] == 0:
                queue.append(d)
    return {v: depth[comp[v]] for v in quiver.vertices}


def layout(quiver, relations=()):
    depth = layer(quiver)
    maxd = max(depth.values(), default=0)
    columns = []
    positions = {}
    for d in range(maxd + 1):
        col = [v for v in quiver.vertices if depth[v] == d]  # Quiver.vertices order
        columns.append(tuple(col))
        k = len(col)
        for i, v in enumerate(col):
            positions[v] = (d, Fraction(k - 1, 2) - i)
    # bundle arrows by (src, tgt); loops split off; parallels get symmetric bends
    bundles = {}
    loops = []
    for name, (s, t) in quiver.arrows.items():
        if s == t:
            loops.append(name)
        else:
            bundles.setdefault((s, t), []).append(name)
    loop_routes = []
    loop_counter = {}
    for name in loops:
        s = quiver.source(name)
        j = loop_counter.get(s, 0)
        loop_counter[s] = j + 1
        loop_routes.append(LoopRoute(name=name, at=s, angle_deg=90 - 60 * j))
    edges = []
    for (s, t), names in bundles.items():
        k = len(names)
        for i, name in enumerate(names):
            if k == 1:
                edges.append(EdgeRoute(name=name, src=s, tgt=t, kind="straight",
                                       bend=Fraction(0)))
            else:
                edges.append(EdgeRoute(name=name, src=s, tgt=t, kind="parallel",
                                       bend=Fraction(k - 1 - 2 * i, 4)))
    return LayoutData(
        positions=positions,
        edges=tuple(edges),
        loops=tuple(loop_routes),
        relations=tuple(repr(r) if not isinstance(r, str) else r for r in relations),
        columns=tuple(columns),
    )
