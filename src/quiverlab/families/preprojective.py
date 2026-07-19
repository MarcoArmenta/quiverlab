"""Preprojective algebra Pi(Q) of a Dynkin quiver (spec §3.4). Double every edge,
impose one mesh relation per vertex. Paths read left-to-right. General route."""
from quiverlab.combinat.quiver import Quiver
from quiverlab.families.dynkin import dynkin_quiver


def PreprojectiveAlgebra(type_or_quiver, field=None):
    base = type_or_quiver if isinstance(type_or_quiver, Quiver) else \
        dynkin_quiver(type_or_quiver, "linear")
    arrows = {}
    star = {}
    for name, (s, t) in base.arrows.items():
        arrows[name] = (s, t)
        arrows[name + "s"] = (t, s)          # a* : t -> s
        star[name] = name + "s"
    Q = Quiver(list(base.vertices), arrows)
    # mesh relation at vertex v: sum_{a: s(a)=v} a a*  -  sum_{b: t(b)=v} b* b  = 0
    rels = []
    for v in base.vertices:
        pos = [f"{a}*{star[a]}" for a, (s, t) in base.arrows.items() if s == v]   # a a*
        neg = [f"{star[b]}*{b}" for b, (s, t) in base.arrows.items() if t == v]   # b* b
        terms = pos + [f"-{p}" for p in neg]
        if not terms:
            continue
        rel = " + ".join(terms).replace("+ -", "- ")
        rels.append(rel)
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("preprojective", "chouhy_solotar", "assem_book")
    return A
