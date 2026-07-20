"""Dynkin/Euclidean diagrams -> quivers with a chosen orientation (spec §3.4).
Edges are undirected pairs; orientation turns each into an arrow. Default 'linear':
i->j for an edge {i,j} with i<j in the standard labeling."""
import re

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError

_TYPE = re.compile(r"^(~|t)?([ADE])(\d+)$")


def _edges(letter, n):
    if letter == "A":
        return [(i, i + 1) for i in range(1, n)]
    if letter == "D":
        if n < 4:
            raise QuiverlabError(f"D{n} needs n >= 4", hint="D4, D5, ...")
        return [(i, i + 1) for i in range(1, n - 1)] + [(n - 2, n)]
    if letter == "E":
        if n not in (6, 7, 8):
            raise QuiverlabError(f"E{n} is not a diagram", hint="E6, E7, E8")
        chain = [(i, i + 1) for i in range(1, n - 1)]     # 1-2-...-(n-1)
        return chain + [(3, n)]                           # branch node 3 -> extra vertex n
    raise QuiverlabError(f"unknown Dynkin letter {letter!r}", hint="A, D, E")


def _euclidean_edges(letter, n):
    if letter == "A":                                     # ~A_n: cycle on 0..n
        return [(i, i + 1) for i in range(n)] + [(0, n)]
    # ~D, ~E labelings per Kac; implement the ones used in the family tour as needed.
    raise QuiverlabError(f"Euclidean ~{letter}{n} not yet tabulated",
                         hint="use ~A_n, or pass an explicit Quiver")


def dynkin_quiver(type_str, orientation="linear"):
    m = _TYPE.match(type_str)
    if not m:
        raise QuiverlabError(f"cannot parse diagram type {type_str!r}",
                             hint="examples: 'A5', 'D4', 'E6', '~A3'")
    euclid, letter, n = bool(m.group(1)), m.group(2), int(m.group(3))
    edges = _euclidean_edges(letter, n) if euclid else _edges(letter, n)
    verts = sorted({v for e in edges for v in e})
    arrows = {}
    for k, (u, v) in enumerate(edges):
        name = f"e{u}{v}"
        rev = f"e{v}{u}"
        # An orientation dict keys each edge by an arrow name; accept the edge under
        # EITHER endpoint ordering (e{u}{v} or e{v}{u}) so a user may pass, e.g.,
        # "e20": (2, 0) for the canonical edge (0, 2). Without this the reversed key
        # silently misses and the requested (e.g. cyclic) orientation is lost.
        if isinstance(orientation, dict) and (name in orientation or rev in orientation):
            s, t = orientation.get(name, orientation.get(rev))
        elif orientation == "reverse":
            s, t = (v, u) if u < v else (u, v)
        else:                                             # "linear"
            s, t = (u, v) if u < v else (v, u)
        arrows[name] = (s, t)
    return Quiver(verts, arrows)
