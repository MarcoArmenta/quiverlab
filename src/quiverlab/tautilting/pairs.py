"""Support tau-tilting pairs, certified per instance (Plan 45 / C4, AIR 2014 Def 0.1).

A pair (M, P): M a basic tau-rigid module, P = (+)_{v in support} P_v the killed
projectives, Hom(P, M) = 0, and |M| + |support| = |Q_0| (the "support tau-tilting"
full-rank condition). (M, emptyset) with |M| = n is a genuine tau-tilting module; the
initial pair is (A_A, emptyset); the terminal pair (0, Q_0). Every constructor validates
all four axioms and refuses loudly (never a silent wrong pair). Over QQ / GF(32003): the
axioms lean on decompose + is_isomorphic (char caveat rigorous over char 0 or char > dim).

Identity note: a :class:`Module` is not hashable by value, so although the dataclass is
``frozen=True`` we NEVER rely on dataclass ``__hash__`` over the modules. The canonical
identity of a pair is :meth:`SupportTauTiltingPair.g_key` -- a ``frozenset`` of int
tuples, fully hashable -- which (AIR: g-vectors determine the pair) is a complete
invariant. The exchange-graph BFS memoizes on ``g_key()``, never on the dataclass."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.modules.decompose import is_indecomposable
from quiverlab.modules.hom import hom_dim, is_isomorphic
from quiverlab.modules.morphism import direct_sum
from quiverlab.tautilting.rigid import g_columns
from quiverlab.tautilting.rigid import g_matrix as _g_matrix


@dataclass(frozen=True)
class SupportTauTiltingPair:
    algebra: object
    summands: tuple
    support: frozenset

    def g_matrix(self):
        return _g_matrix(self)

    def g_key(self):
        """The canonical dedup key: the ``frozenset`` of column tuples of the g-matrix
        (order-independent; AIR: the SET of g-vectors determines the pair)."""
        return frozenset(tuple(c) for c in g_columns(self))


def _direct_sum_or_zero(mods):
    if not mods:
        return None
    D, _, _ = direct_sum(*mods)
    return D


def make_pair(A, summands, support, *, budget=512, check=True):
    """Build a :class:`SupportTauTiltingPair`, validating the four AIR axioms
    (``check=True``, the public/user path): (1) M is tau-rigid; (2) each summand is
    indecomposable and pairwise non-iso (basic); (3) Hom(P_v, M) = 0 for every support
    vertex; (4) ``|M| + |support| = |Q_0|`` (full rank). Refuses loudly otherwise.
    ``check=False`` is the fast path for the BFS AFTER the mutation exchange has already
    certified the neighbour (the mutation sequence is its own certificate)."""
    summands = tuple(summands)
    support = frozenset(support)
    n = len(list(A.quiver.vertices))
    if not check:
        return SupportTauTiltingPair(A, summands, support)
    # (4) rank
    if len(summands) + len(support) != n:
        raise QuiverlabError(
            f"support tau-tilting: |M|={len(summands)} + |support|={len(support)} "
            f"!= |Q_0|={n}", hint="a support tau-tilting pair has full rank n")
    # (2) each summand indecomposable, pairwise non-iso
    for Mi in summands:
        if not is_indecomposable(Mi, budget=budget):     # loud off char scope / decomposable
            raise QuiverlabError("support tau-tilting: a summand is decomposable",
                                 hint="pass indecomposable M_i (decompose first)")
    for i in range(len(summands)):
        for j in range(i + 1, len(summands)):
            if is_isomorphic(summands[i], summands[j]):
                raise QuiverlabError("support tau-tilting: repeated summand (not basic)")
    M = _direct_sum_or_zero(summands)
    # (1) tau-rigid: Hom(M, tau M) = 0
    if M is not None and hom_dim(M, M.tau()) != 0:
        raise QuiverlabError("support tau-tilting: M is not tau-rigid "
                             "(Hom(M, tau M) != 0)")
    # (3) support axiom: Hom(P_v, M) = 0 for v in support
    if M is not None:
        for v in support:
            if hom_dim(A.projective(v), M) != 0:
                raise QuiverlabError(
                    f"support tau-tilting: Hom(P_{v}, M) != 0 -- vertex {v} cannot be "
                    "in the support")
    return SupportTauTiltingPair(A, summands, support)


def initial_pair(A):
    """The initial pair ``(A_A, emptyset)`` = the free module (top of the lattice)."""
    verts = list(A.quiver.vertices)
    return make_pair(A, [A.projective(v) for v in verts], support=[])


def terminal_pair(A):
    """The terminal pair ``(0, Q_0)`` = all projectives killed (bottom of the lattice)."""
    verts = list(A.quiver.vertices)
    return make_pair(A, [], support=verts)
