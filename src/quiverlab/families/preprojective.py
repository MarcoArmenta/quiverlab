"""Preprojective algebra Pi(Q) of a Dynkin quiver (spec §3.4). Double every edge,
impose one mesh relation per vertex. Paths read left-to-right. General route.

Pi(Q) is finite-dimensional IFF Q is (an orientation of) a FINITE Dynkin diagram
of type A/D/E. For an affine (Euclidean, ~/t-prefixed) diagram -- or any other
non-Dynkin quiver -- Pi(Q) is INFINITE-dimensional and has no finite Groebner
presentation, so this constructor refuses it rather than launching a doomed
completion:

  * Given a TYPE STRING, an affine prefix (``~`` or ``t``, e.g. ``"~A2"``) is
    rejected up front with a loud, honest infinite-dimensional QuiverlabError
    naming the type; non-A/D/E letters are refused by ``dynkin_quiver``'s parser.
  * Given an explicit QUIVER, the caller warrants the PRECONDITION that its
    underlying graph is a finite Dynkin A/D/E diagram. There is no cheap Dynkin
    classifier to check this (``dynkin.py`` only *generates* diagrams), so the
    backstop is the finiteness certificate downstream: a non-Dynkin quiver is
    refused loudly (NotFiniteDimensionalError, naming the offending arrow cycle,
    or AdmissibilityError if completion outgrows ``degree_bound``) -- never
    silently or unboundedly computed. ``degree_bound`` caps that completion.
"""
import re

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.families.dynkin import dynkin_quiver

# Affine/Euclidean marker on a type string: a leading '~' or 't' (as in dynkin.py's
# _TYPE = ^(~|t)?([ADE])(\d+)$). These diagrams give INFINITE-dimensional Pi(Q).
_AFFINE_TYPE = re.compile(r"^\s*(~|t)([ADE]\d+)\s*$")


def _reject_infinite_type(type_str):
    """Refuse an affine (Euclidean) type string loudly and honestly: its
    preprojective algebra is infinite-dimensional (no finite Groebner basis)."""
    m = _AFFINE_TYPE.match(type_str)
    if m:
        raise QuiverlabError(
            f"the preprojective algebra of the affine (Euclidean) type {type_str!r} is "
            f"infinite-dimensional, so it has no finite Groebner presentation to compute",
            hint="Pi(Q) is finite-dimensional only for finite Dynkin types: pass A_n "
                 "(n>=1), D_n (n>=4), or E6/E7/E8 (drop the '~'/'t' affine prefix)",
        )


def PreprojectiveAlgebra(type_or_quiver, field=None, degree_bound=None):
    """Preprojective algebra Pi(Q).

    type_or_quiver: a finite Dynkin type string ('A5', 'D4', 'E6', ...) or an
        explicit Quiver that is (an orientation of) a finite Dynkin A/D/E diagram.
        Affine/non-Dynkin input is infinite-dimensional and is refused (see module
        docstring for exactly how each form is refused).
    degree_bound: forwarded to the Groebner engine to cap completion (default:
        adaptive). Raise it if certification refuses a large but finite Dynkin case.
    """
    if isinstance(type_or_quiver, Quiver):
        base = type_or_quiver                          # PRECONDITION: finite Dynkin A/D/E
    else:
        _reject_infinite_type(type_or_quiver)          # loud refusal of affine ~/t types
        base = dynkin_quiver(type_or_quiver, "linear")
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
    A = Q.algebra(relations=rels, field=field, degree_bound=degree_bound)
    A._family_citations = ("preprojective", "chouhy_solotar", "assem_book")
    return A
