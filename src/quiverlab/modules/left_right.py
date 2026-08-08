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

from quiverlab.errors import DepthLimitError, QuiverlabError
from quiverlab.modules.hom import hom_dim, is_isomorphic
from quiverlab.modules.injective import injective_resolution

_LR_REFS = ("act_left_right", "aclv_supports", "organising_module_category")


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


@dataclass(frozen=True)
class SupportAlgebra:
    vertices: tuple            # e_lambda (or e_rho) -- the chosen quiver vertices
    dim: int                   # dim_k A_lambda = sum_{x,y in e} hom_dim(P_x, P_y) = dim End(+ P_x)
    algebra: object            # the presented induced-subquiver Algebra (None iff no vertices)
    components: tuple          # tuple of {"vertices": (...), "algebra": Algebra} factors


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


def _ext_injectives_left(U, left_idx):
    """{ i in L_A : tau^{-1} U[i] = 0  OR  tau^{-1} U[i] notin L_A } -- the Ext-injectives of
    add L_A (ACT/[13](3.4): X Ext-injective in add L_A  <=>  tau^{-1}X notin L_A). tau^{-1}X = 0
    IS the injective edge (an injective in L_A is always Ext-injective). QQ-scope (M1)."""
    out = set()
    for i in left_idx:
        Y = U[i].tau_minus()
        if Y.dim == 0 or _index_in_U(U, Y) not in left_idx:
            out.add(i)
    return out


def _ext_projectives_right(U, right_idx):
    """The dual: { i in R_A : tau U[i] = 0  OR  tau U[i] notin R_A } -- the Ext-projectives of
    add R_A (X Ext-projective in add R_A  <=>  tau X notin R_A). tau X = 0 is the projective
    edge (a projective in R_A is always Ext-projective). QQ-scope (M1)."""
    out = set()
    for i in right_idx:
        Y = U[i].tau()
        if Y.dim == 0 or _index_in_U(U, Y) not in right_idx:
            out.add(i)
    return out


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


def _coeff_sign_mag(coeff):
    """(sign in {+1,-1}, magnitude:str) for a Fraction or exact sympy coefficient, as a
    combinat.relations grammar token (coefficient before the arrows, sign folded out)."""
    from fractions import Fraction
    if isinstance(coeff, Fraction):
        neg = coeff < 0
        return (-1 if neg else 1, str(-coeff if neg else coeff))
    s = str(coeff)
    if s.startswith("-"):
        return (-1, s[1:])
    return (1, s)


def _relation_to_string(r):
    """A parsed ``Relation`` (``.terms`` = ``((coeff, word), ...)``) back to a
    combinat.relations grammar string: signs folded into ' + ' / ' - ' separators (a leading
    '-' on a negative first term), coefficient before the arrows. Monomial relations (the
    rad^2 truncation and every kA_n test) round-trip to just the path word."""
    parts = []
    for idx, (coeff, word) in enumerate(r.terms):
        path = "*".join(word)
        sign, mag = _coeff_sign_mag(coeff)
        if idx == 0:
            sep = "" if sign > 0 else "-"
        else:
            sep = " + " if sign > 0 else " - "
        parts.append(f"{sep}{path}" if mag == "1" else f"{sep}{mag}*{path}")
    return "".join(parts)


def _relation_support(r, quiver):
    """The set of vertices a relation touches (source of its first arrow, target of every
    arrow, plus the parallel source/target endpoints)."""
    vs = {r.source, r.target}
    for _coeff, word in r.terms:
        if word:
            vs.add(quiver.source(word[0]))
            for a in word:
                vs.add(quiver.target(a))
    return vs


def _reachability(quiver):
    """reach[a][b] iff there is a directed path a ->* b (reflexive-transitive) in ``quiver``."""
    V = list(quiver.vertices)
    idx = {v: i for i, v in enumerate(V)}
    n = len(V)
    reach = [[i == j for j in range(n)] for i in range(n)]
    for (s, t) in quiver.arrows.values():
        reach[idx[s]][idx[t]] = True
    _transitive_closure(reach)
    return V, idx, reach


def _is_convex(A, verts):
    """``verts`` is convex in A's quiver: every vertex on a directed path between two chosen
    vertices is itself chosen (so the induced full subquiver loses/creates no relation)."""
    V, idx, reach = _reachability(A.quiver)
    vset = set(verts)
    for x in verts:
        for y in verts:
            for z in V:
                if z not in vset and reach[idx[x]][idx[z]] and reach[idx[z]][idx[y]]:
                    return False
    return True


def _induced_subquiver_algebra(A, verts):
    """The induced FULL subquiver on ``verts`` (arrows with both endpoints chosen) with the
    restricted relations (those whose support lies inside ``verts``), presented over A's field
    -- convexity guarantees no relation is lost or created."""
    from quiverlab.combinat.quiver import Quiver
    vset = set(verts)
    new_arrows = {name: (s, t) for name, (s, t) in A.quiver.arrows.items()
                  if s in vset and t in vset}
    new_rels = [_relation_to_string(r) for r in (A.relations or [])
                if _relation_support(r, A.quiver) <= vset]
    B = Quiver(list(verts), new_arrows).algebra(relations=new_rels, field=A.domain)
    return B, new_arrows


def _quiver_components(verts, arrows):
    """Connected components of the induced subquiver (UNDIRECTED adjacency), each a sorted
    tuple of vertices, in first-appearance order over ``verts``."""
    adj = {v: set() for v in verts}
    for (s, t) in arrows.values():
        adj[s].add(t)
        adj[t].add(s)
    seen = set()
    comps = []
    for start in verts:
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        block = []
        while stack:
            v = stack.pop()
            block.append(v)
            for w in adj[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(tuple(sorted(block, key=verts.index)))
    return comps


def _support_algebra(A, U, verts, *, projectives=True):
    """A_lambda = End_A(+ P_x : x in e_lambda) ~= e_lambda A e_lambda, the induced full convex
    subcategory (Task C). REFUSES loudly if A has no quiver presentation (before any vertex is
    read). Certifies per instance: e_lambda convex, and the presented induced-subquiver dim ==
    dim End (= sum_{x,y} hom_dim(gen_x, gen_y), cross-checked against end_algebra)."""
    if A.quiver is None:
        raise QuiverlabError(
            "support algebra: A has no quiver presentation, so the induced subquiver "
            "e A e cannot be built and the vertex projectives P_x are unavailable",
            hint="feed a quiver-presented kQ/I (the GUI/webapp always do); a "
                 "structure-constants-only algebra has no support-algebra presentation")
    verts = sorted(verts)
    if not verts:                                          # no projective/injective in the part
        return SupportAlgebra(vertices=(), dim=0, algebra=None, components=())
    build = A.projective if projectives else A.injective
    gens = [build(x) for x in verts]
    if not _is_convex(A, verts):
        raise QuiverlabError(
            f"support algebra: the support vertices {verts} are not convex in A's quiver "
            "(a directed path between two of them leaves the set) -- the induced-subquiver "
            "presentation would drop or invent relations", hint="report this presentation")
    B, new_arrows = _induced_subquiver_algebra(A, verts)
    from quiverlab.modules.endomorphism import end_algebra
    from quiverlab.modules.morphism import direct_sum
    end = sum(hom_dim(gens[i], gens[j])
              for i in range(len(gens)) for j in range(len(gens)))   # = dim End(+ P_x)
    D = direct_sum(*gens)[0] if len(gens) > 1 else gens[0]
    ealg = end_algebra(D).dim                              # a 2nd, structure-constant dim End
    if not (B.dim == end == ealg):
        raise QuiverlabError(
            f"support algebra certificate failed: presented induced-subquiver dim {B.dim}, "
            f"sum_{{x,y}} hom_dim {end}, end_algebra dim {ealg} disagree -- a convexity or "
            "relation-restriction bug", hint="report this presentation")
    components = []
    for comp_verts in _quiver_components(verts, new_arrows):
        comp_alg, _ = _induced_subquiver_algebra(A, list(comp_verts))
        components.append({"vertices": comp_verts, "algebra": comp_alg})
    return SupportAlgebra(vertices=tuple(verts), dim=B.dim, algebra=B,
                          components=tuple(components))


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
    ext_inj = _ext_injectives_left(U, L)                       # Task B
    ext_proj = _ext_projectives_right(U, R)                    # Task B
    e_lambda = [v for v in A.quiver.vertices if proj_place[v] in ("L", "both")]
    e_rho = [v for v in A.quiver.vertices if inj_place[v] in ("R", "both")]
    left_support = _support_algebra(A, U, e_lambda, projectives=True)    # Task C
    right_support = _support_algebra(A, U, e_rho, projectives=False)     # Task C
    return LeftRightAtlas(
        A, sel(L), sel(R), sel(L & R), sel(set(range(len(U))) - (L | R)),
        ext_injectives_left=sel(ext_inj), ext_projectives_right=sel(ext_proj),
        left_support=left_support, right_support=right_support,
        projective_placement=proj_place, injective_placement=inj_place,
        universe_size=len(U), is_complete=True, status="complete", note=ar.note or "",
        pd_le_1=tuple(pd_ok), id_le_1=tuple(id_ok),
        _modules=tuple(U), _leq=tuple(tuple(r) for r in leq))


def _sorted_dimvec(dv):
    return {str(w): int(n) for w, n in sorted(dv.items(), key=lambda kv: str(kv[0]))}


def _rec_json(r):
    return {"name": r["name"], "dimvec": _sorted_dimvec(r["dimvec"])}


def _support_json(sa):
    if sa is None:
        return None
    return {"vertices": list(sa.vertices), "dim": int(sa.dim),
            "components": [{"vertices": list(c["vertices"])} for c in sa.components]}


def _empty_lr_block(n, status, note):
    return {"kind": "left_right_parts", "n": int(n), "complete": False, "status": status,
            "universe_size": None, "left": [], "right": [], "intersection": [],
            "complement": [], "ext_injectives_left": [], "ext_projectives_right": [],
            "left_support": None, "right_support": None,
            "projective_placement": {}, "injective_placement": {}, "note": note,
            "references": list(_LR_REFS)}


def left_right_parts_block(A, *, budget=256):
    """The JSON block shared byte-for-byte by both runners (Plan 55, Task E). An ALGEBRA-level
    kind (schema v1, no module block): the module-category atlas -- both parts named with
    S_v/P_v/I_v, the complement (P61 laura datum), the Ext-injectives (P60 hooks) and the two
    support algebras. `Module`s and the reachability matrix are stripped. A char-scope
    identification refusal (large GF(p) is_isomorphic) is surfaced as a loud `error` field
    (the Plan-30 per-entry precedent -- never a 500)."""
    n = len(A.quiver.vertices) if A.quiver is not None else 0
    try:
        atlas = left_right_parts(A, budget=budget)
    except (QuiverlabError, DepthLimitError) as e:
        block = _empty_lr_block(n, "error", None)
        block["error"] = str(e)
        return block
    if not atlas.is_complete:
        return _empty_lr_block(n, atlas.status, atlas.note or None)
    verts = list(A.quiver.vertices)
    return {
        "kind": "left_right_parts", "n": int(n), "complete": True, "status": atlas.status,
        "universe_size": int(atlas.universe_size),
        "left": [_rec_json(r) for r in atlas.left],
        "right": [_rec_json(r) for r in atlas.right],
        "intersection": [_rec_json(r) for r in atlas.intersection],
        "complement": [_rec_json(r) for r in atlas.complement],
        "ext_injectives_left": [_rec_json(r) for r in atlas.ext_injectives_left],
        "ext_projectives_right": [_rec_json(r) for r in atlas.ext_projectives_right],
        "left_support": _support_json(atlas.left_support),
        "right_support": _support_json(atlas.right_support),
        "projective_placement": {str(v): atlas.projective_placement[v] for v in verts},
        "injective_placement": {str(v): atlas.injective_placement[v] for v in verts},
        "note": atlas.note or None,
        "references": list(_LR_REFS),
    }
