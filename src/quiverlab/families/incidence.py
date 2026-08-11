"""Incidence algebra kP of a finite poset = Hasse quiver modulo commutativity
(all parallel paths equal). General (non-monomial) route."""
from quiverlab.families.poset import Poset


def _all_paths(Q):
    """Every directed path (arrow-name tuple) in the Hasse quiver (it is acyclic)."""
    paths = []
    for a in Q.arrows:
        stack = [(a,)]
        while stack:
            w = stack.pop()
            paths.append(w)
            for b in Q.arrows:
                if Q.target(w[-1]) == Q.source(b):
                    stack.append(w + (b,))
    return paths


def IncidenceAlgebra(poset_or_covers, elements=None, field=None):
    P = poset_or_covers if isinstance(poset_or_covers, Poset) else \
        Poset(poset_or_covers, elements)
    Q, _ = P.hasse_quiver()
    by_ends = {}
    for w in _all_paths(Q):
        if len(w) >= 2:
            by_ends.setdefault((Q.word_source(w), Q.word_target(w)), []).append(w)
    rels = []
    for (_st), group in by_ends.items():
        base = group[0]
        for other in group[1:]:                       # base == other (commutativity)
            rels.append("*".join(base) + " - " + "*".join(other))
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("incidence", "assem_book")
    # Provenance for the Plan-75 fast path: HH^*(kP) = H^*(Delta(P); k) (Cibils 1989)
    # only applies when we KNOW the algebra came from a poset. Stashing P here is what
    # lets `incidence_cohomology` route to the order complex WITHOUT ever having to
    # RECOGNIZE an arbitrary kQ/I as an incidence algebra (which it must never guess).
    A._poset = P
    return A
