"""Positive roots of the Tits form for hereditary Dynkin quivers (Plan 38 / C2).

For a Dynkin quiver kQ (no relations, underlying graph an ADE diagram) the
positive roots of the Tits form are exactly the dimension vectors of the
indecomposable representations (Gabriel's theorem). They form a FINITE set,
enumerated here as the reflection closure of the simple roots under the Weyl
group: s_i(d)_i = -d_i + sum_{j~i} m_ij d_j (m_ij = number of edges between i
and j in the underlying graph), s_i(d)_j = d_j otherwise. Only non-negative
vectors are roots; the closure is finite for Dynkin type (no search cap needed).

Affine / wild quivers have infinitely many positive roots and are refused
loudly."""
from quiverlab.errors import QuiverlabError


def positive_roots(A):
    """All positive roots (dimension vectors, vertex order) of the Tits form of a
    hereditary Dynkin algebra A. Raises QuiverlabError for anything else."""
    if A.quiver is None:
        raise QuiverlabError(
            "positive_roots needs the quiver presentation",
            hint="build the algebra via Quiver.algebra(...)")
    from quiverlab.invariants.dynkin_type import dynkin_type
    from quiverlab.invariants.recognizers import is_hereditary
    dt = dynkin_type(A.quiver)
    if not (is_hereditary(A) and dt is not None and dt[0] in ("A", "D", "E")):
        raise QuiverlabError(
            "positive roots enumerated only for Dynkin hereditary type; "
            "affine/wild have infinitely many",
            hint="A must be a path algebra kQ with Q an ADE diagram "
                 f"(here: relations={'present' if A.relations else 'none'}, "
                 f"type={dt})")

    verts = list(A.quiver.vertices)
    n = len(verts)
    idx = {v: k for k, v in enumerate(verts)}
    m = [[0] * n for _ in range(n)]
    for (s, t) in A.quiver.arrows.values():
        if s == t:
            continue
        m[idx[s]][idx[t]] += 1
        m[idx[t]][idx[s]] += 1

    def reflect(d, i):
        e = list(d)
        e[i] = -d[i] + sum(m[i][j] * d[j] for j in range(n))
        return tuple(e)

    seen, frontier = set(), []
    for i in range(n):
        r = tuple(1 if j == i else 0 for j in range(n))
        if r not in seen:
            seen.add(r)
            frontier.append(r)
    while frontier:
        d = frontier.pop()
        for i in range(n):
            e = reflect(d, i)
            if all(x >= 0 for x in e) and any(x > 0 for x in e) and e not in seen:
                seen.add(e)
                frontier.append(e)
    return sorted(seen)
