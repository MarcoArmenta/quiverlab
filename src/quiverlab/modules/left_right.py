"""Left and right parts of the module category, and the support algebras (Plan 55 / R15,
Assem-Coelho-Trepode J. Algebra 281 (2004); Assem-Castonguay-Lanzilotta-Vargas arXiv:1102.1188).

L_A = { M in ind A : pd L <= 1 for every predecessor L of M }, closed under predecessors;
R_A dually (successors, id <= 1). A predecessor of M is the source of a path of non-zero
morphisms between indecomposables ending at M -- the reflexive-transitive closure of
"hom_dim(X, Y) > 0" over the finite (representation-finite) indecomposable universe knitted
by P41. Exact over every Domain; representation-finite scope only (the knit refuses
self-injective and rep-infinite input -- surfaced, never a partial atlas)."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.modules.hom import hom_dim, is_isomorphic
from quiverlab.modules.injective import injective_resolution


@dataclass
class LeftRightAtlas:
    algebra: object
    left: list           # [ {"index","name","dimvec"} , ... ]  X in L_A
    right: list
    intersection: list
    complement: list     # ind A \ (L_A u R_A)          -- P61 laura datum (need NOT be empty)
    ext_injectives_left: list = ()   # populated in Task B (default empty)
    ext_projectives_right: list = ()  # populated in Task B
    left_support: object = None       # SupportAlgebra, populated in Task C
    right_support: object = None      # populated in Task C
    projective_placement: dict = None  # {v: "L"|"R"|"both"|"neither"}   -- P61 ada check
    injective_placement: dict = None
    universe_size: int = 0
    is_complete: bool = False
    status: str = "error"
    note: str = ""
    pd_le_1: tuple = ()  # addendum (P61): pd_le_1[i] iff pd(U[i]) <= 1, index-aligned
    id_le_1: tuple = ()  # addendum (P61): id_le_1[i] iff id(U[i]) <= 1, index-aligned
    _modules: tuple = ()  # actual Module universe in index order (P60/P61 in-process)
    _leq: tuple = ()      # reachability matrix (self-cert hook)

    def in_union(self, index) -> bool:
        return any(r["index"] == index for r in self.left) or \
               any(r["index"] == index for r in self.right)


def _universe(A, budget):
    ar = A.ar_quiver(budget_modules=budget)
    U = [v["module"] for v in ar.vertices]
    names = [v["name"] for v in ar.vertices]
    recs = [{"index": i, "name": v["name"], "dimvec": v["dimvec"]}
            for i, v in enumerate(ar.vertices)]
    return ar, U, names, recs


def _leq_matrix(U):
    """DEFAULT predecessor relation: reflexive-transitive closure of hom_dim(X,Y) > 0.
    Definition-faithful and exact over every Domain (no is_isomorphic search)."""
    n = len(U)
    leq = [[i == j for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and hom_dim(U[i], U[j]) > 0:
                leq[i][j] = True
    _transitive_closure(leq)
    return leq


def _ar_reachability(ar, n):
    """The CROSS-CHECK predecessor relation: reflexive-transitive closure of the knit's
    irreducible-map arrows. Coincides with the Hom closure for rep-finite (rad^infty = 0)."""
    leq = [[i == j for j in range(n)] for i in range(n)]
    for (i, j) in ar.arrows:
        leq[i][j] = True
    _transitive_closure(leq)
    return leq


def _transitive_closure(leq):
    n = len(leq)
    for k in range(n):
        rk = leq[k]
        for i in range(n):
            if leq[i][k]:
                ri = leq[i]
                for j in range(n):
                    if rk[j]:
                        ri[j] = True


def _pd_le_1(M):
    return M.projective_resolution(2).betti(2) == 0            # P_2 = 0  <=>  pd <= 1


def _id_le_1(M):
    return injective_resolution(M, 2).betti(2) == 0            # E^2 = 0  <=>  id <= 1


def _left_indices(leq, pd_ok):
    n = len(leq)
    return {j for j in range(n)
            if all((not leq[i][j]) or pd_ok[i] for i in range(n))}


def _right_indices(leq, id_ok):
    n = len(leq)
    return {j for j in range(n)
            if all((not leq[j][i]) or id_ok[i] for i in range(n))}


def _index_in_U(U, M):
    """First j with dimvec(U[j]) == dimvec(M) and is_isomorphic(U[j], M); None if absent.

    QQ / char-0 scope (M1): is_isomorphic is decisive over char 0, but over large
    GF(p)/GF(p^n) it is positive-only and RAISES when it cannot exhibit an isomorphism --
    that loud QuiverlabError propagates (the honest refusal, never a silent wrong index)."""
    if M is None or M.dim == 0:
        return None
    dv = M.dimension_vector()
    for j, X in enumerate(U):
        if X.dim == M.dim and X.dimension_vector() == dv and is_isomorphic(X, M):
            return j
    return None


def _placement(A, U, L, R, *, kind):
    """{v: "both"|"L"|"R"|"neither"} for the standard projectives (kind="projective") or
    injectives (kind="injective") of A, by locating P_v / I_v in the universe U and reading
    its L/R membership. QQ-scope (M1); a generator the rep-finite universe must contain but
    which cannot be located is an internal inconsistency -> loud refusal."""
    build = A.projective if kind == "projective" else A.injective
    place = {}
    for v in A.quiver.vertices:
        idx = _index_in_U(U, build(v))
        if idx is None:
            raise QuiverlabError(
                f"left_right_parts: the indecomposable {kind[0].upper()}_{v} is not in the "
                "knitted universe (rep-finite invariant violated)",
                hint="report this presentation; the AR knit should contain every "
                     "indecomposable projective and injective")
        in_l, in_r = idx in L, idx in R
        place[v] = ("both" if in_l and in_r else "L" if in_l
                    else "R" if in_r else "neither")
    return place


def left_right_parts(A, *, budget=256):
    # Fast representation-infinite refusal (hereditary non-Dynkin): the projective-seeded
    # knit would run away on a rep-infinite algebra (a multi-minute hang before the module
    # budget trips). Refuse before it does -- mirrors degeneration.py's guard (2026-08-05).
    try:
        hereditary = A.is_hereditary()
    except QuiverlabError:
        hereditary = False
    if hereditary and A.form_type() != "finite":
        return LeftRightAtlas(
            A, [], [], [], [], is_complete=False, status="unsupported",
            note=f"hereditary of {A.dynkin_type()} type: representation-infinite "
                 "(infinitely many indecomposables) -- no finite left/right parts")
    ar, U, names, recs = _universe(A, budget)
    if not ar.is_complete:
        return LeftRightAtlas(A, [], [], [], [], is_complete=False,
                              status=ar.status, note=ar.note or "")
    leq = _leq_matrix(U)                                        # DEFAULT: Hom closure
    pd_ok = [_pd_le_1(M) for M in U]
    id_ok = [_id_le_1(M) for M in U]
    L = _left_indices(leq, pd_ok)
    R = _right_indices(leq, id_ok)
    sel = lambda S: [recs[i] for i in sorted(S)]
    proj_place = _placement(A, U, L, R, kind="projective")
    inj_place = _placement(A, U, L, R, kind="injective")
    return LeftRightAtlas(
        A, sel(L), sel(R), sel(L & R), sel(set(range(len(U))) - (L | R)),
        projective_placement=proj_place, injective_placement=inj_place,
        universe_size=len(U), is_complete=True, status="complete", note=ar.note or "",
        pd_le_1=tuple(pd_ok), id_le_1=tuple(id_ok),
        _modules=tuple(U), _leq=tuple(tuple(r) for r in leq))
    # ext_injectives_left / ext_projectives_right and left_support / right_support keep their
    # ()/None defaults here; Task B and Task C fill them (M3: no forward reference).
