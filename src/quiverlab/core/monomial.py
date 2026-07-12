"""Monomial kQ/I -> Algebra with a certified finite basis (spec §3.3, §5 component 4).
Basis = trivial paths e_v + irreducible paths (words avoiding every forbidden
word as a contiguous subword). Finiteness is decided by a suffix-window
automaton; an infinite family is reported with an explicit arrow cycle."""
from collections import deque

from quiverlab.core.algebra import Algebra
from quiverlab.errors import AdmissibilityError, NotFiniteDimensionalError, RelationError


def _contains_forbidden(word, forbidden):
    return any(
        word[i: i + len(f)] == f
        for f in forbidden
        for i in range(len(word) - len(f) + 1)
    )


def _automaton(quiver, forbidden):
    """States (vertex, window) reachable by irreducible words; window = last r-1 arrows."""
    r = max((len(f) for f in forbidden), default=1)
    starts = [(v, ()) for v in quiver.vertices]
    graph = {}
    seen = set(starts)
    dq = deque(starts)
    while dq:
        state = dq.popleft()
        v, w = state
        outs = []
        for a, (s, t) in quiver.arrows.items():
            if s != v:
                continue
            word = w + (a,)
            if any(len(f) <= len(word) and word[-len(f):] == f for f in forbidden):
                continue
            ns = (t, word[-(r - 1):] if r > 1 else ())
            outs.append((a, ns))
            if ns not in seen:
                seen.add(ns)
                dq.append(ns)
        graph[state] = outs
    return starts, graph


def _find_cycle(graph):
    """Return the arrow labels of a cycle in the state graph, or None."""
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {s: WHITE for s in graph}
    for root in graph:
        if color[root] != WHITE:
            continue
        stack = [(root, iter(graph[root]), None)]
        path = []  # arrows along the current DFS path
        color[root] = GRAY
        while stack:
            state, it, _ = stack[-1]
            step = next(it, None)
            if step is None:
                color[state] = BLACK
                stack.pop()
                if path:
                    path.pop()
                continue
            arrow, ns = step
            if color[ns] == GRAY:
                # unwind to where ns sits on the stack
                idx = next(i for i, fr in enumerate(stack) if fr[0] == ns)
                return path[idx:] + [arrow]
            if color[ns] == WHITE:
                color[ns] = GRAY
                stack.append((ns, iter(graph[ns]), arrow))
                path.append(arrow)
    return None


def irreducible_paths(quiver, forbidden):
    """All irreducible words, sorted by (length, word). Raises
    NotFiniteDimensionalError (with a cycle) when infinitely many exist."""
    forbidden = [tuple(f) for f in forbidden]
    starts, graph = _automaton(quiver, forbidden)
    cycle = _find_cycle(graph)
    if cycle is not None:
        raise NotFiniteDimensionalError(
            "kQ/I is infinite-dimensional: irreducible paths grow forever along the "
            "cycle " + " -> ".join(cycle),
            hint="add relations killing a power of this cycle (monomial), or check your quiver",
        )
    words = []
    for st in starts:
        stack = [(st, ())]
        while stack:
            state, word = stack.pop()
            for arrow, ns in graph[state]:
                nw = word + (arrow,)
                words.append(nw)
                stack.append((ns, nw))
    return sorted(set(words), key=lambda w: (len(w), w))


def build_monomial_algebra(quiver, relations, field):
    """relations: parsed Relation objects, all monomial, lengths >= 2."""
    for rel in relations:
        if rel.min_length < 2:
            raise AdmissibilityError(
                f"relation {rel!r} has a path of length {rel.min_length}: the ideal is "
                "not inside the square of the arrow ideal",
                hint="admissible relations use paths of length >= 2",
            )
    parsed_pool = [field.parse_entry(0), field.parse_entry(1)]
    dom = field.make_domain(parsed_pool)
    for rel in relations:
        (coeff, word), = rel.terms  # monomial: exactly one term
        if dom.is_zero(dom.coerce(coeff)):
            raise RelationError(
                f"coefficient {coeff} of {rel!r} vanishes in {dom.name}, so the relation is 0 = 0",
                hint="a monomial relation's coefficient must be nonzero in the chosen field; "
                     "drop the relation or change the field",
            )
    forbidden = [w for rel in relations for _, w in rel.terms]  # monomial: one word each
    words = irreducible_paths(quiver, forbidden)
    trivial = [("e", v) for v in quiver.vertices]
    basis = trivial + [("p", w) for w in words]
    index = {b: i for i, b in enumerate(basis)}
    m = len(basis)

    def src(b):
        return b[1] if b[0] == "e" else quiver.word_source(b[1])

    def tgt(b):
        return b[1] if b[0] == "e" else quiver.word_target(b[1])

    def prod(x, y):
        if tgt(x) != src(y):
            return None
        if x[0] == "e":
            return y
        if y[0] == "e":
            return x
        w = x[1] + y[1]
        if _contains_forbidden(w, [tuple(f) for f in forbidden]):
            return None
        return ("p", w)

    zero, one = dom.zero(), dom.one()
    T = [[[zero] * m for _ in range(m)] for _ in range(m)]
    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            p = prod(bi, bj)
            vec = [zero] * m
            if p is not None:
                vec[index[p]] = one
            T[i][j] = vec
    unit = [zero] * m
    for v in quiver.vertices:
        unit[index[("e", v)]] = one
    labels = [f"e_{b[1]}" if b[0] == "e" else "*".join(b[1]) for b in basis]
    return Algebra(dom, T, unit, basis_labels=labels, _quiver=quiver,
                   _relations=list(relations))
