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

    def algebra(self, relations=(), field=None):
        """Build kQ/I over the field (default CC). Monomial presentations are
        certified and lowered; general relations arrive with the Groebner engine."""
        from quiverlab.combinat.relations import parse_relations
        from quiverlab.core.monomial import build_monomial_algebra

        if field is None:
            from quiverlab.fields import CC
            field = CC
        rels = parse_relations(list(relations), self)
        if all(r.is_monomial for r in rels):
            return build_monomial_algebra(self, rels, field)
        raise NotImplementedError(
            "general (non-monomial) relations arrive with the Groebner engine "
            "(Plan 03); this build certifies monomial presentations only"
        )

    def __repr__(self):
        lines = [f"Quiver with vertices {self.vertices} and arrows:"]
        lines += [f"  {s} --{a}--> {t}" for a, (s, t) in self.arrows.items()]
        return "\n".join(lines)
