"""Coverings (Plan 56): the combinatorial fundamental group pi1(Q, I) of a
presentation, its abelianization by exact integer Smith normal form, the
Hom(pi1, k+) Hurewicz count, the three-valued ``is_simply_connected`` verdict,
and the R16 separation-condition / ``is_strongly_simply_connected`` recognizer
(the clean certificate P62 consumes).

pi1(Q, I) is the group of ~_I-classes of closed walks (Martinez-Villa--de la
Pena; Le Meur arXiv:math/0503302). Its generators are the non-tree arrows of a
spanning forest; its relators come from the *minimal relations* of I, which are a
basis of the finite-dimensional space ``I / (rad.I + I.rad)`` -- computed here by
EXACT LINEAR ALGEBRA with lifting, NEVER by a greedy reduction of the stored
generators (which is unsound: a stored generator can be irredundant yet not a
minimal relation). The abelianization pi1^ab is a function of ~_I alone, hence
basis-independent, and is obtained by exact ZZ Smith normal form.

The INTRINSIC pi1 (inverse limit over connected gradings, Cibils--Redondo--Solotar
arXiv:0906.3069) is NOT bounded-computable in general and is refused loudly.

Exact only: the group theory is over ZZ; Hom(pi1, k+) counts free rank plus
char-divisible torsion (an integer/characteristic computation).
"""
from dataclasses import dataclass, field

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import invariant_factors

from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import nullspace, rref

_PATH_BUDGET = 500_000


# ---------------------------------------------------------------------------
# guards
# ---------------------------------------------------------------------------
def _require_quiver(A, what):
    if A.quiver is None or A.relations is None:
        raise QuiverlabError(
            f"{what} needs the quiver presentation",
            hint="build the algebra via Quiver.algebra(...); structure-constant "
                 "algebras carry no quiver presentation, so pi1(Q, I) is undefined")


# ---------------------------------------------------------------------------
# minimal relations = a basis of I / (rad.I + I.rad), per (source, target) block
# ---------------------------------------------------------------------------
def _paths_between(quiver, N):
    """dict (x, y) -> sorted list of composable words x -> y of length in [2, N]
    (single-arrow paths are excluded: an admissible I sits inside rad^2)."""
    from collections import defaultdict
    by_src = defaultdict(list)
    for name in sorted(quiver.arrows):
        s, t = quiver.arrows[name]
        by_src[s].append((name, t))
    result = defaultdict(list)
    cur = [((name,), quiver.arrows[name][0], quiver.arrows[name][1])
           for name in sorted(quiver.arrows)]          # length 1
    total = 0
    length = 1
    while length < N:
        nxt = []
        for (w, s, t) in cur:
            for (name, t2) in by_src.get(t, []):
                nxt.append((w + (name,), s, t2))
        length += 1
        total += len(nxt)
        if total > _PATH_BUDGET:
            raise QuiverlabError(
                "fundamental_group: path enumeration exceeded the internal budget",
                hint="the quiver has an enormous number of bounded-length paths; "
                     "this is a resource guard, not a mathematical refusal")
        for (w, s, t) in nxt:
            result[(s, t)].append(w)
        cur = nxt
    for k in result:
        result[k].sort()
    return dict(result)


def _kernel(M, ncols, dom):
    """Basis of {c : M c = 0} in path-coordinates. M has rows = A-basis words,
    columns = paths; an EMPTY M (every path maps to 0 in A) means the whole space
    is the kernel."""
    if not M:
        out = []
        for j in range(ncols):
            v = [dom.zero()] * ncols
            v[j] = dom.one()
            out.append(v)
        return out
    return nullspace(M, dom)


def _block_W(paths, rs, dom):
    """W_{x,y} = I ∩ span(paths) = kernel of the path -> A map, as basis vectors
    over `paths` (in order)."""
    n = len(paths)
    nfs = [rs.normal_form(p) for p in paths]
    basiswords = sorted({w for nf in nfs for w in nf})
    idx = {w: i for i, w in enumerate(basiswords)}
    M = [[dom.zero()] * n for _ in basiswords]
    for j, nf in enumerate(nfs):
        for w, c in nf.items():
            M[idx[w]][j] = c
    return _kernel(M, n, dom)


def _extend(w, src_paths, dst_idx, name, prepend, N, dom):
    """Prepend (or append) arrow `name` to each path of vector `w` (over
    `src_paths`), landing in the coordinate system `dst_idx` (path -> column).
    Paths that would exceed length N are dropped (they lie in rad^N ⊆ I)."""
    v = [dom.zero()] * len(dst_idx)
    for k, coeff in enumerate(w):
        if dom.is_zero(coeff):
            continue
        p = src_paths[k]
        np_ = (name,) + p if prepend else p + (name,)
        if len(np_) > N:
            continue
        j = dst_idx.get(np_)
        if j is not None:
            v[j] = dom.add(v[j], coeff)
    return v


def _R_rows(x, y, paths_xy, paths, W, quiver, N, dom):
    """(rad.I + I.rad) ∩ span(paths_xy): span of one-arrow extensions of shorter
    minimal relations. Returns generating rows over paths_xy (not row-reduced)."""
    dst_idx = {p: i for i, p in enumerate(paths_xy)}
    rows = []
    for name in sorted(quiver.arrows):
        s, t = quiver.arrows[name]
        if s == x:                                  # prepend name: x -> t, w over (t, y)
            src = paths.get((t, y))
            if src:
                for w in W.get((t, y), []):
                    rows.append(_extend(w, src, dst_idx, name, True, N, dom))
        if t == y:                                  # append name: s -> y, w over (x, s)
            src = paths.get((x, s))
            if src:
                for w in W.get((x, s), []):
                    rows.append(_extend(w, src, dst_idx, name, False, N, dom))
    return rows


def _reduce_against(w, rows, dom):
    """w reduced modulo the row space of `rows` (zeros at every pivot column)."""
    if not rows:
        return list(w)
    R, pivots = rref(rows, dom)
    y = list(w)
    for r, pc in enumerate(pivots):
        c = y[pc]
        if not dom.is_zero(c):
            y = [dom.sub(a, dom.mul(c, b)) for a, b in zip(y, R[r])]
    return y


def _minimal_reps(W, R_rows, dom):
    """Reduced-support basis of span(W) / span(R): reduce W mod R, then RREF. Each
    returned row is a minimal relation whose nonzero coordinates are the parallel
    paths it glues."""
    Wred = [_reduce_against(w, R_rows, dom) for w in W]
    Wred = [w for w in Wred if any(not dom.is_zero(x) for x in w)]
    if not Wred:
        return []
    R, pivots = rref(Wred, dom)
    return [R[i] for i in range(len(pivots))]


@dataclass(frozen=True)
class _RelationData:
    quiver: object
    tree: frozenset
    parent: dict
    generators: tuple
    glue: tuple            # tuple of tuples-of-paths (support >= 2), across all blocks
    counts: dict           # {(src, tgt): #minimal relations} (>= 1, incl. monomial)


def _relation_data(A, base=None):
    """The single shared computation behind fundamental_group and
    minimal_relation_counts: minimal relations of I by exact I/(rad.I + I.rad)
    linear algebra, plus the spanning forest and generators."""
    _require_quiver(A, "fundamental_group")
    from quiverlab.resolutions_cs.build import reduction_system_of
    quiver = A.quiver
    dom = A.domain
    rs = reduction_system_of(A)
    N = A.loewy_length()
    tree, parent = quiver.spanning_forest(root=base)
    generators = tuple(sorted(set(quiver.arrows) - tree))

    paths = _paths_between(quiver, N) if N >= 2 else {}
    W = {blk: _block_W(ps, rs, dom) for blk, ps in paths.items()}

    glue = []
    counts = {}
    for (x, y), ps in paths.items():
        Wxy = W[(x, y)]
        if not Wxy:
            continue
        R_rows = _R_rows(x, y, ps, paths, W, quiver, N, dom)
        reps = _minimal_reps(Wxy, R_rows, dom)
        if not reps:
            continue
        counts[(x, y)] = len(reps)
        for row in reps:
            support = tuple(ps[k] for k, c in enumerate(row) if not dom.is_zero(c))
            if len(support) >= 2:                    # a monomial glues nothing
                glue.append(support)
    return _RelationData(quiver, frozenset(tree), parent, generators,
                         tuple(glue), counts)


# ---------------------------------------------------------------------------
# the fundamental group + abelianization
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class FundamentalGroup:
    generators: tuple          # non-tree arrow labels (the base-point loops)
    relators: tuple            # words in the generators (strings), from minimal relations
    free_rank: int             # rank of the free part of pi1^ab
    invariant_factors: tuple   # the SNF factors d_i > 1 (torsion of pi1^ab)
    components: int            # connected components of Q (pi1^ab is a direct sum over them)
    base_vertex: object
    presentation_note: str
    relator_words: tuple = field(default=(), repr=False, compare=False)
    # relator_words: tuple of tuples of (generator, +1|-1) -- the structured form
    # Tietze-lite reads; the public payload uses `relators` (strings).

    def abelianization_repr(self) -> str:
        parts = []
        if self.free_rank == 1:
            parts.append("Z")
        elif self.free_rank > 1:
            parts.append(f"Z^{self.free_rank}")
        parts += [f"Z/{d}" for d in self.invariant_factors]
        return " (+) ".join(parts) if parts else "0"

    def hom_to_additive_dim(self, char: int) -> int:
        """dim_k Hom(pi1^ab, (k, +)) = free_rank + #{d_i : char > 0 and char | d_i}."""
        extra = 0
        if char and char > 0:
            extra = sum(1 for d in self.invariant_factors if d % char == 0)
        return self.free_rank + extra


def _gen_word(path, genset):
    """The forward generator subsequence of a directed path (tree arrows drop)."""
    return tuple(a for a in path if a in genset)


def _vec(path, generators, genset):
    """Integer vector in Z^{generators}: occurrences of each generator in the path."""
    v = [0] * len(generators)
    pos = {g: i for i, g in enumerate(generators)}
    for a in path:
        if a in genset:
            v[pos[a]] += 1
    return v


def fundamental_group(A, base=None) -> FundamentalGroup:
    """The presentation fundamental group pi1(Q, I): generators = non-tree arrows,
    relators from the minimal relations of I, and pi1^ab by exact ZZ Smith normal
    form. Always computable (no refusal beyond presentation-less input)."""
    data = _relation_data(A, base=base)
    quiver = data.quiver
    generators = data.generators
    genset = set(generators)
    components = len(quiver.undirected_components())
    base_vertex = base if base is not None else (
        quiver.vertices[0] if quiver.vertices else None)

    rows = []
    relators = []
    relator_words = []
    for support in data.glue:
        supp = sorted(support)
        u0 = supp[0]
        for ui in supp[1:]:
            rows.append([a - b for a, b in zip(_vec(u0, generators, genset),
                                               _vec(ui, generators, genset))])
            w0 = _gen_word(u0, genset)
            wi = _gen_word(ui, genset)
            structured = tuple((g, 1) for g in w0) + tuple((g, -1) for g in reversed(wi))
            relator_words.append(structured)
            relators.append(_relator_str(structured))

    if rows:
        M = sp.Matrix(rows)
        rank = M.rank()
        factors = tuple(int(d) for d in invariant_factors(M, domain=ZZ) if int(d) > 1)
        free_rank = len(generators) - rank
    else:
        factors = ()
        free_rank = len(generators)

    return FundamentalGroup(
        generators=generators,
        relators=tuple(relators),
        free_rank=free_rank,
        invariant_factors=factors,
        components=components,
        base_vertex=base_vertex,
        presentation_note="presentation invariant of (Q, I); NOT of the algebra A",
        relator_words=tuple(relator_words),
    )


def _relator_str(structured):
    if not structured:
        return "1"
    return " ".join(g if s > 0 else f"{g}^-1" for g, s in structured)


def intrinsic_fundamental_group(A):
    """The intrinsic fundamental group is NOT bounded-computable -- always refused."""
    _require_quiver(A, "intrinsic_fundamental_group")
    raise QuiverlabError(
        "the intrinsic fundamental group (inverse limit over connected gradings) is "
        "not bounded-computable in general; quiverlab emits the PRESENTATION group "
        "pi1(Q, I) and its abelianization instead",
        hint="e.g. the presentation pi1(k[x]/(x^p)) = Z (monomial loop) differs from "
             "the intrinsic pi1 = Z x C_p in characteristic p (Cibils-Redondo-Solotar, "
             "arXiv:0906.3069) -- the C_p torsion is invisible to any presentation")


# ---------------------------------------------------------------------------
# P62 contract (addendum): per-(src,tgt) minimal-relation counts = Tits r_ij
# ---------------------------------------------------------------------------
def minimal_relation_counts(A) -> dict:
    """{(src, tgt): count} for every ordered vertex pair with >= 1 minimal relation
    (= dim_k e_tgt (I / (rad.I + I.rad)) e_src, the Tits-form r_ij). Presentation-less
    A: loud QuiverlabError. Empty dict for a hereditary (relation-free) algebra."""
    return dict(_relation_data(A).counts)


# ---------------------------------------------------------------------------
# bypasses (Le Meur) + Tietze-lite triviality test
# ---------------------------------------------------------------------------
_BYPASS_BUDGET = 100_000


def _all_paths(quiver, s, t):
    """All oriented paths s -> t of length >= 1 (arrow-name tuples). Acyclic => finite;
    a cyclic quiver is length-capped at |V| and refused loudly if it blows the budget."""
    from collections import defaultdict
    by_src = defaultdict(list)
    for name in sorted(quiver.arrows):
        a, b = quiver.arrows[name]
        by_src[a].append((name, b))
    N = len(quiver.vertices)
    out = []
    frontier = [((name,), b) for name, b in by_src.get(s, [])]
    length = 1
    while frontier:
        nxt = []
        for (w, cur) in frontier:
            if cur == t:
                out.append(w)
            if length < N:
                for (name, b) in by_src.get(cur, []):
                    nxt.append((w + (name,), b))
        if len(out) > _BYPASS_BUDGET or len(nxt) > _BYPASS_BUDGET:
            raise QuiverlabError(
                "bypasses: oriented-path enumeration exceeded the internal budget",
                hint="the quiver has too many bounded-length paths (resource guard)")
        frontier = nxt
        length += 1
    return out


def bypasses(A):
    """[(arrow, path_word), ...]: an arrow alpha and an oriented path parallel to it
    (same endpoints, distinct). Bounded enumeration (loud budget on blow-up)."""
    _require_quiver(A, "bypasses")
    quiver = A.quiver
    out = []
    for name in sorted(quiver.arrows):
        s, t = quiver.arrows[name]
        for u in _all_paths(quiver, s, t):
            if u != (name,):
                out.append((name, u))
    return out


def has_double_bypass(A) -> bool:
    """True iff some (alpha, u, beta, v) with (alpha,u),(beta,v) bypasses and the
    arrow beta appearing inside the path u (Le Meur)."""
    bps = bypasses(A)
    for (_alpha, u) in bps:
        for (beta, _v) in bps:
            if beta in u:
                return True
    return False


def _free_reduce(word):
    out = []
    for tok in word:
        if out and out[-1][0] == tok[0] and out[-1][1] == -tok[1]:
            out.pop()
        else:
            out.append(tok)
    return out


def _inv(word):
    return [(g, -s) for g, s in reversed(word)]


def _tietze_trivial(generators, relator_words):
    """A SUFFICIENT triviality test: repeatedly eliminate a generator occurring EXACTLY
    ONCE in some relator (a sound Tietze/Nielsen transformation), substituting it into
    the others. All generators eliminated => the group is trivial. It NEVER returns True
    for a group it did not genuinely trivialise (it only ever certifies triviality)."""
    from collections import Counter
    gens = set(generators)
    rels = [_free_reduce(list(r)) for r in relator_words]
    rels = [r for r in rels if r]
    progress = True
    while gens and progress:
        progress = False
        for ri in range(len(rels)):
            r = rels[ri]
            cnt = Counter(g for g, _s in r)
            target = next((g for g in gens if cnt.get(g, 0) == 1), None)
            if target is None:
                continue
            pos = next(i for i, (g, _s) in enumerate(r) if g == target)
            _g, ss = r[pos]
            pre, post = r[:pos], r[pos + 1:]
            # r = pre * g^ss * post = 1  =>  g = inv(pre)*inv(post) (ss=+1) or post*pre (ss=-1)
            value = _inv(pre) + _inv(post) if ss == 1 else post + pre
            newrels = []
            for rj, other in enumerate(rels):
                if rj == ri:
                    continue
                w = []
                for (g, s) in other:
                    if g == target:
                        w += value if s == 1 else _inv(value)
                    else:
                        w.append((g, s))
                w = _free_reduce(w)
                if w:
                    newrels.append(w)
            rels = newrels
            gens.discard(target)
            progress = True
            break
    return not gens


# ---------------------------------------------------------------------------
# separation condition (R16) + strongly-simply-connected recognizer
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Separation:
    holds: object          # True / False / None (undecided_char)
    witness: object        # on False: {vertex, summand_supports:[set, set]}
    reason: str


def _transitive_predecessors(quiver, a):
    """All x != a with a directed path x ~> a."""
    preds = set()
    stack = [a]
    while stack:
        v = stack.pop()
        for p in quiver.predecessors(v):
            if p != a and p not in preds:
                preds.add(p)
                stack.append(p)
    return preds


def separation_condition(A) -> Separation:
    """For triangular A: at each vertex a the DISTINCT indecomposable summands of
    rad P_a must have supports in DISTINCT connected components of
    Q_a = the full subquiver on {v : a is UNREACHABLE from v} = Q minus a and its
    transitive predecessor closure. Raises loudly ONLY on non-triangular input; a
    decompose char refusal is CAUGHT as holds=None (undecided_char)."""
    _require_quiver(A, "separation_condition")
    quiver = A.quiver
    if not quiver.is_acyclic():
        raise QuiverlabError(
            "separation_condition requires a triangular algebra (acyclic quiver)",
            hint="the separation condition / strong simple connectivity are defined "
                 "for triangular algebras (Skowronski 1993); this quiver has an "
                 "oriented cycle")
    from quiverlab.modules.decompose import decompose
    verts = list(quiver.vertices)
    for a in verts:
        radP = A.projective(a).radical()
        if radP.dim == 0:
            continue
        try:
            summands = decompose(radP)
        except QuiverlabError:
            return Separation(None, None,
                              "undecided_char (decompose refused over char<=dim)")
        Qa_verts = set(verts) - _transitive_predecessors(quiver, a) - {a}
        comps = quiver.induced_subquiver(Qa_verts).undirected_components()
        comp_of = {}
        for i, comp in enumerate(comps):
            for v in comp:
                comp_of[v] = i
        claimed = {}
        for (M, _mult) in summands:
            dv = M.dimension_vector()
            supp = frozenset(v for v, d in dv.items() if d > 0)
            comp_ids = {comp_of[v] for v in supp if v in comp_of}
            for cid in comp_ids:
                if cid in claimed and claimed[cid] != supp:
                    return Separation(
                        False,
                        {"vertex": a,
                         "summand_supports": [set(claimed[cid]), set(supp)]},
                        f"vertex {a}: summands share a component")
                if cid in claimed:            # same support, distinct summand: still shares
                    return Separation(
                        False,
                        {"vertex": a, "summand_supports": [set(supp), set(supp)]},
                        f"vertex {a}: summands share a component")
                claimed[cid] = supp
    return Separation(True, None, "separated")


@dataclass(frozen=True)
class StrongSimpleConnectivity:
    verdict: object        # True / False / None (budget_exceeded OR undecided_char)
    witness: object        # on False: {convex_subset, vertex, summand_supports}
    reason: str
    checked_convex: int


def _reachable(quiver):
    """x -> frozenset of vertices reachable from x (including x)."""
    from collections import deque
    succ = {v: set() for v in quiver.vertices}
    for (s, t) in quiver.arrows.values():
        if s != t:
            succ[s].add(t)
    reach = {}
    for v in quiver.vertices:
        seen = {v}
        dq = deque([v])
        while dq:
            x = dq.popleft()
            for y in succ[x]:
                if y not in seen:
                    seen.add(y)
                    dq.append(y)
        reach[v] = frozenset(seen)
    return reach


def _convex_closure(quiver, T, reach):
    S = set(T)
    changed = True
    while changed:
        changed = False
        add = set()
        for x in S:
            rx = reach[x]
            for y in S:
                for u in rx:                        # x ~> u and u ~> y  =>  u in [x, y]
                    if u not in S and y in reach[u]:
                        add.add(u)
        if add:
            S |= add
            changed = True
    return frozenset(S)


def _convex_subsets(quiver, budget):
    """Every non-empty convex vertex subset (complete: each convex C is reached by
    adding its vertices under convex closure, staying inside C). Returns None if the
    count exceeds `budget`."""
    from collections import deque
    verts = list(quiver.vertices)
    reach = _reachable(quiver)
    seen = set()
    dq = deque()
    for v in verts:
        cc = _convex_closure(quiver, {v}, reach)
        if cc not in seen:
            seen.add(cc)
            dq.append(cc)
    while dq:
        if len(seen) > budget:
            return None
        S = dq.popleft()
        for v in verts:
            if v in S:
                continue
            cc = _convex_closure(quiver, S | {v}, reach)
            if cc not in seen:
                if len(seen) >= budget:
                    return None
                seen.add(cc)
                dq.append(cc)
    return seen


def _restrict_relations(A, S):
    """repr-strings of the stored relations whose every path stays inside S."""
    quiver = A.quiver
    out = []
    for rel in (A.relations or []):
        ok = True
        for _c, w in rel.terms:
            for name in w:
                s, t = quiver.arrows[name]
                if s not in S or t not in S:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            out.append(repr(rel))
    return out


def _subalgebra(A, S):
    from quiverlab.resolutions_cs._fieldshim import field_for_domain
    sub_q = A.quiver.induced_subquiver(S)
    field = field_for_domain(A.domain)
    return sub_q.algebra(relations=_restrict_relations(A, S), field=field)


def is_strongly_simply_connected(A, convex_budget=20000) -> StrongSimpleConnectivity:
    """R16 (Skowronski 1993): triangular A is strongly simply connected iff every
    full convex subcategory satisfies the separation condition. Three-valued: the
    first HARD separation failure => False (+ witness); a per-subcategory char refusal
    (undecided_char) or exceeding convex_budget => None; else True. Raises loudly ONLY
    on non-triangular TOP input -- NEVER from inside the sweep."""
    _require_quiver(A, "is_strongly_simply_connected")
    if not A.quiver.is_acyclic():
        raise QuiverlabError(
            "is_strongly_simply_connected requires a triangular algebra (acyclic quiver)",
            hint="strong simple connectivity is defined for triangular algebras "
                 "(Skowronski 1993); this quiver has an oriented cycle")
    subsets = _convex_subsets(A.quiver, convex_budget)
    if subsets is None:
        return StrongSimpleConnectivity(
            None, None, "budget_exceeded", convex_budget)
    ordered = sorted(subsets, key=lambda s: (len(s), sorted(map(str, s))))
    checked = 0
    undecided = None
    for S in ordered:
        sub = _subalgebra(A, S)
        sep = separation_condition(sub)
        checked += 1
        if sep.holds is False:
            w = dict(sep.witness or {})
            w["convex_subset"] = sorted(S, key=str)
            return StrongSimpleConnectivity(
                False, w, f"{sorted(S, key=str)}@vertex {w.get('vertex')} fails", checked)
        if sep.holds is None and undecided is None:
            undecided = sorted(S, key=str)
    if undecided is not None:
        return StrongSimpleConnectivity(
            None, None, f"undecided_char@{undecided}", checked)
    return StrongSimpleConnectivity(True, None, "all convex separated", checked)


# ---------------------------------------------------------------------------
# three-valued is_simply_connected
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SimpleConnectivity:
    verdict: object                 # True / False / None (inconclusive -- Adian-Rabin)
    reason: str
    witness: object                 # on False: {"kind": ...}
    abelianization: FundamentalGroup
    strongly: object                # StrongSimpleConnectivity or None


def is_simply_connected(A, strong="auto", convex_budget=20000) -> SimpleConnectivity:
    """Three-valued simple-connectivity verdict (sound: True is NEVER derived from a
    failed search). False guards first (disconnected / oriented cycle / nontrivial
    pi1^ab), then the True routes cheapest-first R1 (tree) -> R3 (no-bypass Le Meur,
    Tietze-lite) -> R2 (separation / R16). `strong` controls the R16 certificate:
    "auto" computes it only if R1/R3 do not decide; True always computes it (the
    P62 gate); False never computes it (cheap routes only, so this NEVER raises on an
    R1/R3-decidable input). Loud only on presentation-less input."""
    _require_quiver(A, "is_simply_connected")
    quiver = A.quiver
    g = fundamental_group(A)

    def _cert():
        return is_strongly_simply_connected(A, convex_budget)

    # -- decidable False witnesses ------------------------------------------
    if g.components > 1:
        return SimpleConnectivity(
            False, "the quiver is not connected (simple connectivity requires a "
            "connected quiver)", {"kind": "disconnected"}, g, None)
    if not quiver.is_acyclic():
        return SimpleConnectivity(
            False, "the quiver has an oriented cycle (simple connectivity requires "
            "no oriented cycles)", {"kind": "oriented_cycle"}, g, None)
    if g.free_rank > 0 or g.invariant_factors:
        return SimpleConnectivity(
            False, "the stored presentation has pi1(Q,I)^ab != 0, so pi1 != 1 for "
            "this presentation", {"kind": "nontrivial_pi1ab",
                                  "free_rank": g.free_rank,
                                  "invariant_factors": list(g.invariant_factors)},
            g, None)

    # -- True routes, cheapest first ----------------------------------------
    # R1: tree (Betti 0) => pi1(Q, J) = 1 for EVERY presentation.
    if len(g.generators) == 0:
        cert = _cert() if strong is True else None
        return SimpleConnectivity(
            True, "the underlying graph is a tree (first Betti number 0): pi1 is "
            "trivial for every presentation (R1)", None, g, cert)

    # R3: triangular + no bypasses => pi1 presentation-independent (Le Meur); trivialise.
    bps = bypasses(A)
    tietze_ok = None
    if len(bps) == 0:
        tietze_ok = _tietze_trivial(g.generators, g.relator_words)
        if tietze_ok:
            cert = _cert() if strong is True else None
            return SimpleConnectivity(
                True, "triangular with no bypasses (pi1 is presentation-independent, "
                "Le Meur) and the presentation trivialises (R3)", None, g, cert)

    # R2: separation / strong simple connectivity (the expensive convex sweep).
    if strong in (True, "auto"):
        cert = _cert()
        if cert.verdict is True:
            return SimpleConnectivity(
                True, "strongly simply connected via the separation condition (R2)",
                None, g, cert)
        return SimpleConnectivity(
            None, _none_reason(bps, tietze_ok, cert), None, g, cert)

    # strong is False: cheap routes only -> honest None.
    return SimpleConnectivity(None, _none_reason(bps, tietze_ok, None), None, g, None)


def _none_reason(bps, tietze_ok, cert):
    if bps:
        base = ("bypasses are present, so pi1(Q,I) is presentation-dependent and no "
                "privileged-presentation route applies")
    elif tietze_ok is False:
        base = ("pi1(Q,I)^ab = 0 but the finite presentation did not trivialise "
                "(a perfect-pi1 corner)")
    else:
        base = "no decidable sufficient criterion fired"
    if cert is not None and cert.verdict is None:
        base += f"; separation inconclusive ({cert.reason})"
    return (base + "; verdict undecided (triviality of a finitely presented group is "
            "undecidable -- Adian-Rabin)")
