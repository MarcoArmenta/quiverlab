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
