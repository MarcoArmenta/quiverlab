"""Quiver = finite directed multigraph with named arrows. Paths are tuples of
arrow names read LEFT TO RIGHT: ('a', 'b') means first a, then b, and requires
target(a) == source(b) (Assem-Simson-Skowronski convention)."""
import re

from quiverlab.errors import RelationError

_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


class Quiver:
    def __init__(self, vertices, arrows):
        vertices = list(vertices)
        if len(set(vertices)) != len(vertices):
            raise RelationError("duplicate vertices", hint="each vertex must appear once")
        self.vertices = vertices
        vset = set(vertices)
        self.arrows = {}
        for name, ends in dict(arrows).items():
            if not (isinstance(name, str) and _NAME.match(name)):
                raise RelationError(
                    f"bad arrow name {name!r}",
                    hint="arrow names must be identifiers like a, b2, alpha (they appear in relation strings)",
                )
            try:
                s, t = ends
            except (TypeError, ValueError):
                raise RelationError(f"arrow {name!r}: endpoints must be a (source, target) pair",
                                    hint="e.g. arrows={'a': (1, 2)}") from None
            if s not in vset or t not in vset:
                raise RelationError(f"arrow {name!r}: endpoint not a vertex",
                                    hint=f"vertices are {vertices}")
            self.arrows[name] = (s, t)

    # -- accessors -----------------------------------------------------------
    def source(self, name):
        return self.arrows[name][0]

    def target(self, name):
        return self.arrows[name][1]

    def word_source(self, word):
        return self.arrows[word[0]][0] if word else None

    def word_target(self, word):
        return self.arrows[word[-1]][1] if word else None

    def compose_ok(self, word) -> bool:
        return all(self.target(a) == self.source(b) for a, b in zip(word, word[1:]))

    def is_acyclic(self) -> bool:
        adj = {v: [] for v in self.vertices}
        for s, t in self.arrows.values():
            adj[s].append(t)
        WHITE, GRAY, BLACK = 0, 1, 2
        color = {v: WHITE for v in self.vertices}
        for start in self.vertices:
            if color[start] != WHITE:
                continue
            stack = [(start, iter(adj[start]))]
            color[start] = GRAY
            while stack:
                v, it = stack[-1]
                nxt = next(it, None)
                if nxt is None:
                    color[v] = BLACK
                    stack.pop()
                elif color[nxt] == GRAY:
                    return False
                elif color[nxt] == WHITE:
                    color[nxt] = GRAY
                    stack.append((nxt, iter(adj[nxt])))
        return True

    def is_connected(self) -> bool:
        """True iff the underlying undirected graph is connected (Plan 38)."""
        from quiverlab.invariants.dynkin_type import is_connected
        return is_connected(self)

    # -- graph primitives for pi1 / coverings (Plan 56) ----------------------
    # Pure combinatorics; NO invariants import (they are used from combinat).
    def predecessors(self, v) -> set:
        """Proper immediate predecessors {s : (s -> v) arrow, s != v} (loops excluded)."""
        return {s for (s, t) in self.arrows.values() if t == v and s != v}

    def successors(self, v) -> set:
        """Proper immediate successors {t : (v -> t) arrow, t != v} (loops excluded)."""
        return {t for (s, t) in self.arrows.values() if s == v and t != v}

    def undirected_components(self, vertices=None) -> list:
        """Connected components (list of frozensets) of the underlying UNDIRECTED
        graph induced on `vertices` (default: all vertices). Loops are ignored for
        adjacency. Component order follows the quiver's vertex order (deterministic)."""
        order = list(self.vertices) if vertices is None else [
            v for v in self.vertices if v in set(vertices)]
        vset = set(order)
        # a caller may pass vertices not in this quiver: keep them as singletons
        extra = [v for v in (vertices or ()) if v not in vset]
        order = order + extra
        vset = set(order)
        adj = {v: set() for v in vset}
        for (s, t) in self.arrows.values():
            if s == t or s not in vset or t not in vset:
                continue
            adj[s].add(t)
            adj[t].add(s)
        seen = set()
        comps = []
        for v in order:
            if v in seen:
                continue
            seen.add(v)
            stack = [v]
            comp = {v}
            while stack:
                x = stack.pop()
                for y in adj[x]:
                    if y not in seen:
                        seen.add(y)
                        comp.add(y)
                        stack.append(y)
            comps.append(frozenset(comp))
        return comps

    def spanning_forest(self, root=None):
        """One BFS spanning tree per connected component of the underlying
        undirected graph. Returns (tree_arrow_names:set, parent_map:dict). A tree
        arrow is chosen once per newly reached vertex; parallel/back arrows and
        loops stay non-tree, so pi1 generators = set(arrows) - tree_arrow_names.
        Deterministic: components are entered in vertex order (``root`` first if
        given), incident arrows are scanned in sorted-name order.

        parent_map[v] is (parent_vertex, arrow_name, forward) for a non-root v
        (forward = True iff the arrow points parent -> v) and None for a root."""
        from collections import deque
        inc = {v: [] for v in self.vertices}
        for name in sorted(self.arrows):
            s, t = self.arrows[name]
            inc[s].append((name, t, True))
            if t != s:
                inc[t].append((name, s, False))
        order = list(self.vertices)
        if root is not None:
            order = [root] + [v for v in self.vertices if v != root]
        tree = set()
        parent = {}
        visited = set()
        for start in order:
            if start in visited:
                continue
            visited.add(start)
            parent[start] = None
            dq = deque([start])
            while dq:
                v = dq.popleft()
                for (name, w, forward) in inc[v]:
                    if w not in visited:
                        visited.add(w)
                        tree.add(name)
                        parent[w] = (v, name, forward)
                        dq.append(w)
        return tree, parent

    def induced_subquiver(self, vertices) -> "Quiver":
        """The FULL subquiver on `vertices`: keep the vertex order, keep exactly the
        arrows whose BOTH endpoints lie in the set (boundary arrows are dropped)."""
        vset = set(vertices)
        verts = [v for v in self.vertices if v in vset]
        arrows = {name: (s, t) for name, (s, t) in self.arrows.items()
                  if s in vset and t in vset}
        return Quiver(verts, arrows)

    def algebra(self, relations=(), field=None, degree_bound=None, trace=None):
        """Build kQ/I over the field (default CC). Monomial presentations route
        through the Plan-01 monomial path; general (non-monomial) relations route
        through the Groebner engine (Plan 03). Paths compose LEFT TO RIGHT.

        degree_bound: cap on ambiguity length during completion (default: adaptive).
        trace: optional list; step events (Dispatch, ReductionStep) are appended
        (inert hooks -- formal trace rendering is Plan 07)."""
        from quiverlab.combinat.relations import parse_relations

        if field is None:
            from quiverlab.fields import CC
            field = CC
        rels = parse_relations(list(relations), self)
        if all(r.is_monomial for r in rels):
            if trace is not None:
                from quiverlab.groebner.events import Dispatch
                trace.append(Dispatch(route="monomial",
                                      reason="every relation is a single monomial",
                                      n_relations=len(rels)))
            from quiverlab.core.monomial import build_monomial_algebra
            return build_monomial_algebra(self, rels, field)
        if trace is not None:
            from quiverlab.groebner.events import Dispatch
            trace.append(Dispatch(route="groebner",
                                  reason="at least one relation is non-monomial",
                                  n_relations=len(rels)))
        from quiverlab.groebner.lower import groebner_algebra
        return groebner_algebra(self, rels, field, degree_bound=degree_bound, trace=trace)

    def __repr__(self):
        lines = [f"Quiver with vertices {self.vertices} and arrows:"]
        lines += [f"  {s} --{a}--> {t}" for a, (s, t) in self.arrows.items()]
        return "\n".join(lines)
