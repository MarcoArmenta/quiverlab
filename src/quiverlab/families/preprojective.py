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

# A finite Dynkin type string, split into (letter, n) exactly as dynkin.py parses it
# (sans the already-rejected affine prefix) -- used only for the auto-bound lookup.
_FINITE_TYPE = re.compile(r"^([ADE])(\d+)$")

# Per-type Groebner degree_bound needed to CERTIFY completion when the caller passes
# no explicit degree_bound. Pi(Q)'s mesh relations are quadratic, so the shipped
# adaptive default (max(8, 2*maxlen+4) = 8) certifies only the smallest cases
# (A2-A4 build at 8); larger Dynkin diagrams complete to longer leading words and
# need a larger bound. Values are the feasibility-probe working bounds (each is >= the
# minimum certifiable bound and, since the reduction system is bound-STABLE above that
# minimum, yields the same certified system any larger bound would). Types NOT listed
# fall through to the adaptive default: A2-A4 build there; D6+ and E6/E7/E8 hit the
# loud AdmissibilityError (with its "raise degree_bound to at least N" hint) -- those
# exceptional/large cases are cluster-scale and deliberately NOT auto-lifted here.
_AUTO_DEGREE_BOUND = {
    ("A", 5): 12, ("A", 6): 14, ("A", 7): 16,
    ("D", 4): 12, ("D", 5): 16,
}


def _auto_degree_bound(type_str):
    """The auto degree_bound for a finite Dynkin type string, or None to leave the
    adaptive default in force (both for unlisted small types and for the large
    exceptional types that must keep raising AdmissibilityError)."""
    m = _FINITE_TYPE.match(type_str)
    if not m:
        return None
    return _AUTO_DEGREE_BOUND.get((m.group(1), int(m.group(2))))


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
    degree_bound: forwarded to the Groebner engine to cap completion. If None (the
        default), a per-type auto-bound is applied for the type strings that the
        shipped adaptive default cannot certify (A5:12, A6:14, A7:16, D4:12, D5:16),
        so PreprojectiveAlgebra("A5"), ("D5"), etc. build with no kwarg. An explicit
        degree_bound always overrides the auto-bound. Large exceptional types (D6+,
        E6/E7/E8) are NOT auto-lifted: they keep the loud AdmissibilityError with its
        bound hint -- certifying them is cluster-scale. An explicit Quiver never gets
        an auto-bound (there is no type string to key on); pass degree_bound yourself.
    """
    if isinstance(type_or_quiver, Quiver):
        base = type_or_quiver                          # PRECONDITION: finite Dynkin A/D/E
    else:
        _reject_infinite_type(type_or_quiver)          # loud refusal of affine ~/t types
        if degree_bound is None:                       # fill in the certifiable bound for A5+/D4+
            degree_bound = _auto_degree_bound(type_or_quiver)
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
