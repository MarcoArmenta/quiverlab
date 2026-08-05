"""Brauer graph algebras (Plan 46 / C5).

A ribbon (multi)graph with vertex multiplicities presents a symmetric special-biserial
``kQ/I``: the quiver has one VERTEX per EDGE of the graph; the arrows are the
consecutive edge-pairs around each graph-vertex (the "special cycle" ``C_v`` of length
``val(v)``); the relations are the standard Brauer presentation (special-cycle gluing,
zero-after-``m_v`` cycles, and cycle-switch quadratics). Presented via
``Quiver.algebra``; per-instance dimension-certified ``dim = sum_v m_v * val(v)^2``
(loud otherwise). Schroll 2018 (survey); Wald-Waschbuesch 1985. Float-free / exact.

Arrow-direction convention (ARBITRATED): the cyclic edge order around a vertex is read
so that consecutive edges ``(e_i, e_{i+1})`` give an arrow ``e_i -> e_{i+1}``; the
Brauer-STAR special case is then byte-Cartan-equal to
``NakayamaAlgebra(n, m*n+1, cyclic=True)`` (the cross-oracle). If that Cartan
disagreed, the reading would be reversed once -- it does not."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError


@dataclass(frozen=True)
class BrauerGraph:
    """A connected ribbon multigraph. ``edges`` is a tuple of unordered endpoint pairs
    ``(u, v)`` (loops ``u == v`` allowed), indexed ``0..E-1``. ``cyclic_order`` maps each
    graph-vertex to a tuple of ``(edge_index, end_tag)`` in ribbon (cyclic) order; a
    loop contributes two ends at its vertex."""
    edges: tuple
    cyclic_order: dict

    def vertices(self):
        vs = set(self.cyclic_order)
        for u, v in self.edges:
            vs.add(u)
            vs.add(v)
        return sorted(vs, key=repr)

    def valency(self, v):
        return len(self.cyclic_order.get(v, ()))

    def validate(self):
        # every edge-end must appear exactly once in exactly one cyclic order (a loop
        # appears twice at its vertex): total occurrences of edge k == 2.
        counts = {k: 0 for k in range(len(self.edges))}
        vertex_ends = {}
        for v, order in self.cyclic_order.items():
            for (k, _tag) in order:
                if k not in counts:
                    raise QuiverlabError(
                        f"BrauerGraph: cyclic order at {v!r} names edge index {k} "
                        f"which does not exist (edges 0..{len(self.edges)-1})",
                        hint="edge indices in cyclic_order must index `edges`")
                counts[k] += 1
                vertex_ends.setdefault(k, []).append(v)
        for k, (u, v) in enumerate(self.edges):
            if counts[k] != 2:
                raise QuiverlabError(
                    f"BrauerGraph: edge {k} = {(u, v)!r} has {counts[k]} ribbon ends "
                    "(need exactly 2 -- one per endpoint, both for a loop)",
                    hint="every edge must appear in the cyclic order of each endpoint")
            ends = sorted(vertex_ends[k], key=repr)
            expect = sorted([u, v], key=repr)
            if ends != expect:
                raise QuiverlabError(
                    f"BrauerGraph: edge {k} = {(u, v)!r} is ordered at vertices {ends} "
                    f"but connects {expect}",
                    hint="the cyclic-order placements must match the edge endpoints")
        if not self.edges:
            raise QuiverlabError("BrauerGraph: need at least one edge")


def _present(graph, mult):
    """Build (Quiver, relations) for the Brauer graph algebra. Q-vertex ``k+1`` is the
    ``k``-th edge (1-indexed, in `edges` order)."""
    edges = graph.edges
    E = len(edges)
    vertices = list(range(1, E + 1))

    def qv(edge_idx):
        return edge_idx + 1

    arrows = {}
    special = {}                       # graph-vertex -> [arrow names] (the cycle C_v)
    for v in graph.vertices():
        order = graph.cyclic_order.get(v, ())
        d = len(order)
        mv = mult[v]
        if d == 0:
            special[v] = []
            continue
        if d == 1 and mv == 1:
            special[v] = []            # truncated leaf: no arrow, no loop
            continue
        cyc = []
        for i in range(d):
            k_i = order[i][0]
            k_next = order[(i + 1) % d][0]
            name = f"c{v}_{i}"
            arrows[name] = (qv(k_i), qv(k_next))
            cyc.append(name)
        special[v] = cyc

    Q = Quiver(vertices, arrows)
    rels = []

    # (1) zero relations: going m_v times around C_v then one more arrow vanishes.
    for v, cyc in special.items():
        if not cyc:
            continue
        d = len(cyc)
        mv = mult[v]
        for start in range(d):
            path = [cyc[(start + k) % d] for k in range(mv * d + 1)]
            rels.append("*".join(path))

    def _start_index(cyc, edge_idx):
        """Position in ``cyc`` of the arrow leaving Q-vertex ``qv(edge_idx)``."""
        target = qv(edge_idx)
        for i, nm in enumerate(cyc):
            if Q.source(nm) == target:
                return i
        return None

    # (2) gluing: for an edge e = {u, w} non-truncated at both ends, the two special
    # cycles agree: C_u^{m_u} (from e) = C_w^{m_w} (from e).
    for k, (u, w) in enumerate(edges):
        cu, cw = special.get(u, []), special.get(w, [])
        if not cu or not cw or u == w:
            continue
        su, sw = _start_index(cu, k), _start_index(cw, k)
        if su is None or sw is None:
            continue
        pu = [cu[(su + i) % len(cu)] for i in range(mult[u] * len(cu))]
        pw = [cw[(sw + i) % len(cw)] for i in range(mult[w] * len(cw))]
        rels.append("*".join(pu) + " - " + "*".join(pw))

    # (3) switch quadratics: at a shared edge, composing across the two cycles is 0.
    def _in_arrow(cyc, edge_idx):
        target = qv(edge_idx)
        for nm in cyc:
            if Q.target(nm) == target:
                return nm
        return None

    for k, (u, w) in enumerate(edges):
        cu, cw = special.get(u, []), special.get(w, [])
        if not cu or not cw or u == w:
            continue
        su, sw = _start_index(cu, k), _start_index(cw, k)
        if su is None or sw is None:
            continue
        out_u, out_w = cu[su], cw[sw]
        in_u, in_w = _in_arrow(cu, k), _in_arrow(cw, k)
        if in_u is not None and out_w is not None:
            rels.append(f"{in_u}*{out_w}")     # ... into e via C_u, out via C_w = 0
        if in_w is not None and out_u is not None:
            rels.append(f"{in_w}*{out_u}")
    return Q, rels


def BrauerGraphAlgebra(graph, multiplicities, field=None):
    """The Brauer graph algebra of a ribbon graph with vertex ``multiplicities``.
    Presented ``kQ/I``, per-instance dimension-certified ``dim = sum_v m_v*val(v)^2``
    (loud ``QuiverlabError`` otherwise). BGAs are symmetric."""
    graph.validate()
    for v in graph.vertices():
        if v not in multiplicities or int(multiplicities[v]) < 1:
            raise QuiverlabError(
                f"BrauerGraphAlgebra: missing/invalid multiplicity at vertex {v!r}",
                hint="every graph vertex needs an integer multiplicity >= 1")
    mult = {v: int(multiplicities[v]) for v in graph.vertices()}
    Q, rels = _present(graph, mult)
    max_proj = 1
    for k, (u, w) in enumerate(graph.edges):
        du = mult[u] * graph.valency(u) if graph.valency(u) else 0
        dw = mult[w] * graph.valency(w) if graph.valency(w) else 0
        max_proj = max(max_proj, du + dw)
    A = Q.algebra(relations=rels, field=field, degree_bound=max_proj + 2)
    expected = sum(mult[v] * graph.valency(v) ** 2 for v in graph.vertices())
    if A.dim != expected:
        raise QuiverlabError(
            f"BrauerGraphAlgebra: dim {A.dim} != certificate {expected} "
            "(= sum_v m_v*val(v)^2)",
            hint="the ribbon ordering or the multiplicities do not present a Brauer "
                 "graph algebra; check the cyclic edge order around each vertex")
    A._family_citations = ("schroll_brauer", "wald_waschbusch", "assem_book")
    return A
