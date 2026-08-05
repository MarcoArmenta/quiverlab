"""Jacobian algebra of a quiver with potential (Q, W) (Plan 44 / C7).

A potential ``W`` is a k-linear combination of CYCLIC words (oriented cycles in Q). The
cyclic derivative ``d_a W`` sums over the cyclic rotations that place ``a`` first: for each
occurrence of ``a`` in a cyclic word ``w`` at position ``r``, the deletion
``w[r+1:] + w[:r]`` (delete ``a``, keep the rest of the cycle). The Jacobian (DWZ) algebra
is ``Jac(Q, W) = kQ / (d_a W : a in Q_1)``. It may be infinite-dimensional --
``NotFiniteDimensionalError`` is the honest refusal (never a fabricated finiteness claim).
Derksen-Weyman-Zelevinsky, Selecta Math. 2008; surface case Labardini-Fragoso, Proc. LMS
2009. Float-free; coefficients are exact (int / Fraction / sanctioned scalars)."""
from quiverlab.errors import QuiverlabError
from quiverlab.families.trivial_extension import _relation_string   # reuse the emitter

_CITATIONS = ("derksen_weyman_zelevinsky", "labardini", "assem_book")


class Potential:
    """A potential ``W`` on a quiver: a list of ``(coeff, word)`` terms, each ``word`` an
    oriented CYCLE (composable, ``source(first) == target(last)``). Loud on a non-cyclic
    word."""

    def __init__(self, quiver, terms):
        self.quiver = quiver
        self.terms = []
        for coeff, word in terms:
            word = tuple(word)
            if (not word or not quiver.compose_ok(word)
                    or quiver.word_source(word) != quiver.word_target(word)):
                raise QuiverlabError(
                    f"Potential: {word!r} is not a cyclic word in Q",
                    hint="each potential term must be an oriented cycle (composable, "
                         "with source of the first arrow == target of the last)")
            self.terms.append((coeff, word))


def cyclic_derivative(W, arrow):
    """``d_arrow W`` as a relation-grammar string (the parallel sum of the cyclic
    deletions), or ``""`` when ``d_arrow W = 0``. All deletions of a cyclic word share
    source ``target(arrow)`` and target ``source(arrow)`` -- a genuine parallel relation."""
    Q = W.quiver
    out_terms = []
    for coeff, word in W.terms:
        for r in range(len(word)):
            if word[r] == arrow:
                out_terms.append((coeff, word[r + 1:] + word[:r]))   # delete arrow, rotate
    if not out_terms:
        return ""                                                    # d_a W = 0: no relation
    srcs = {Q.word_source(w) for _c, w in out_terms}
    tgts = {Q.word_target(w) for _c, w in out_terms}
    if len(srcs) > 1 or len(tgts) > 1:                               # cannot happen for cyclic W
        raise QuiverlabError(
            f"cyclic_derivative: d_{arrow} W is not a parallel relation "
            f"(sources {sorted(srcs, key=str)}, targets {sorted(tgts, key=str)})",
            hint="the potential terms must be genuine oriented cycles")
    return _relation_string(out_terms)


def _auto_degree_bound(Q):
    """A per-family degree-bound heuristic (``None`` leaves the Groebner route's own
    adaptive bound in force). Mirrors ``preprojective.py``: the shipped adaptive default
    certifies the small cases; non-monomial potentials that need a larger bound raise the
    Groebner ``AdmissibilityError`` with its explicit 'raise degree_bound' hint."""
    return None


def JacobianAlgebra(Q, W, field=None, degree_bound=None):
    """The Jacobian algebra ``kQ / (d_a W : a in Q_1)``. ``NotFiniteDimensionalError`` is
    the honest Jacobian-INFINITE refusal (propagated from ``Quiver.algebra``)."""
    if not isinstance(W, Potential):
        raise QuiverlabError(
            f"JacobianAlgebra needs a Potential, got {type(W).__name__}",
            hint="wrap W in Potential(Q, terms)")
    if W.quiver is not Q:
        raise QuiverlabError(
            "JacobianAlgebra: the potential is over a different quiver than Q",
            hint="build the Potential on the same Quiver object passed here")
    rels = [s for s in (cyclic_derivative(W, a) for a in Q.arrows) if s]
    if degree_bound is None:
        degree_bound = _auto_degree_bound(Q)
    A = Q.algebra(relations=rels, field=field, degree_bound=degree_bound)
    A._family_citations = _CITATIONS
    return A
