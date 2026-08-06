"""Torsion-class lattice + Hasse orientation + bricks / semibricks (Plan 45 / C4,
AIR 2014; DIRRT / Asai brick labelling). Each support tau-tilting pair (M, P) corresponds
to the functorially finite torsion class ``Gen(M)`` (AIR bijection); mutation is a covering
relation of the torsion-class lattice, and each cover is labelled by a BRICK (the wall
normal, King). The Hasse quiver is oriented so ``(A, 0)`` is the unique source (top =
mod A) and ``(0, A)`` the unique sink (bottom = 0). Semibricks are the Asai counterparts:
the set of down-labels of each pair, a bijection with the pairs (the four-way count
identity ``#s-tau-tilt = #f.f. torsion = #semibricks``).

Rigorous over the base field (bricks decide ``End(B) = k`` over an algebraically closed /
char-0 base; the GF(p^n) proper-division-ring caveat is honest-scope). Batteries run over
QQ. Float-free; the wall normals are exact integer vectors (:func:`mutation.wall_normal`)."""
from __future__ import annotations

import weakref

from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.hom import end_dim, is_isomorphic
from quiverlab.modules.morphism import direct_sum


# --------------------------------------------------------------------------- #
# the module universe (indecomposable tau-rigid summands across the exchange graph)
# --------------------------------------------------------------------------- #
# Keyed by the algebra OBJECT (WeakKeyDictionary), NOT id(A): an id can be recycled after
# a GC'd algebra, which would hand a stale universe to a DIFFERENT algebra. The weak key
# is identity-based and drops with the algebra, so there is no cross-instance bleed.
_UNIVERSE_CACHE: "weakref.WeakKeyDictionary" = weakref.WeakKeyDictionary()


def _torsion_universe(A, budget=512):
    """The indecomposable modules appearing as summands across every support tau-tilting
    pair of ``A`` (a sufficient universe to fingerprint the torsion classes -- AIR: the
    Ext-projectives of a class generate it). Deduped by ``is_isomorphic``; cached per
    algebra (weak-keyed). Raises loudly (via the BFS) if ``A`` is tau-tilting-infinite."""
    cached = _UNIVERSE_CACHE.get(A)
    if cached is not None:
        return cached
    from quiverlab.tautilting.mutation import exchange_graph
    eg = exchange_graph(A, budget_pairs=budget)
    if not eg.is_complete:
        from quiverlab.errors import QuiverlabError
        raise QuiverlabError(
            f"torsion universe: exchange graph did not close (status={eg.status}); "
            "tau-tilting-infinite or budget-capped", hint="raise budget / finite algebra")
    uni = []
    for rec in eg.vertices:
        for M in rec["pair"].summands:
            if not any(U.dim == M.dim and U.dimension_vector() == M.dimension_vector()
                       and is_isomorphic(U, M) for U in uni):
                uni.append(M)
    _UNIVERSE_CACHE[A] = uni
    return uni


def _in_gen(M, X):
    """``X in Gen(M)``: the trace of ``M`` in ``X`` (the span of the images of all
    ``M -> X``) fills ``X``. ``M is None`` (the zero module) generates only 0."""
    if X.dim == 0:
        return True
    if M is None:
        return False
    from quiverlab.modules.hom import hom_space
    dom = X.domain
    cols = []
    for h in hom_space(M, X):                          # X.dim x M.dim matrices
        for j in range(len(h[0]) if h else 0):
            cols.append([h[i][j] for i in range(len(h))])
    if not cols:
        return False
    return lm.mat_rank(lm.cols_to_matrix(cols), dom) == X.dim


def _module_of(pair):
    if not pair.summands:
        return None
    D, _, _ = direct_sum(*pair.summands)
    return D


def _dvtuple(dv):
    return tuple((str(v), int(dv[v])) for v in sorted(dv, key=str))


def torsion_class_data(pair):
    """A distinguishing fingerprint of ``Gen(M)``: ``{"gen_dimvecs": <sorted dim-vectors of
    the universe indecomposables in Gen(M)>, "is_full": bool, "is_zero": bool}``. Injective
    on the pairs (the AIR bijection) -- the pair -> torsion-class map is one-to-one -- which
    the tests pin. NOT claimed as the full torsion-class enumeration (honest scope)."""
    A = pair.algebra
    n = len(list(A.quiver.vertices))
    M = _module_of(pair)
    universe = _torsion_universe(A)
    members = [X for X in universe if _in_gen(M, X)]
    gen_dimvecs = sorted(_dvtuple(X.dimension_vector()) for X in members)
    return {
        "gen_dimvecs": gen_dimvecs,
        "is_full": (not pair.support and len(pair.summands) == n),
        "is_zero": (not pair.summands),
    }


def _gen_size(pair):
    return len(torsion_class_data(pair)["gen_dimvecs"])


# --------------------------------------------------------------------------- #
# Hasse orientation
# --------------------------------------------------------------------------- #
def hasse_orientation(eg):
    """Orient every exchange edge ``(i, j)`` (i<j) as ``"down"`` (i -> j downward) or
    ``"up"`` (j -> i downward) by torsion-class INCLUSION -- the endpoint with the larger
    ``Gen`` is higher. ``(A, 0)`` becomes the unique source, ``(0, A)`` the unique sink
    (covers are strict inclusions, so the two endpoints always differ in size)."""
    sizes = {i: _gen_size(rec["pair"]) for i, rec in enumerate(eg.vertices)}
    out = {}
    for (i, j) in eg.arrows:
        out[(i, j)] = "down" if sizes[i] > sizes[j] else "up"
    return out


# --------------------------------------------------------------------------- #
# bricks + semibricks
# --------------------------------------------------------------------------- #
# The wall/edge brick is identified by ISO-CLASS, not by dim-vector. The dim-vector of a
# wall (its King normal) does NOT determine the brick on a non-thin algebra: kZ2/rad^2 has
# two non-isomorphic bricks P1, P2 both of dim-vector (1,1), on OPPOSITE King-stability rays
# of the same hyperplane. A first-dim-vector-match label would collapse them (bricks() would
# undercount, and two genuinely distinct semibricks {P1}, {P2} would merge -- breaking the
# AIR four-way identity #pairs = #semibricks). We disambiguate by torsion-class membership.
def _bricks_with_dimvec(A, dimvec, universe):
    """Every brick (``end_dim == 1``) in ``universe`` whose dim-vector equals ``dimvec`` (a
    vertex-keyed dict), in universe order. Usually a singleton -- but a NON-THIN algebra
    carries several non-isomorphic bricks of one dim-vector, so the wall label is NOT
    determined by the dim-vector alone."""
    if dimvec is None:
        return []
    target = {v: int(dimvec.get(v, 0)) for v in A.quiver.vertices}
    return [X for X in universe if X.dimension_vector() == target and end_dim(X) == 1]


def _edge_brick(A, pair_i, pair_j, dimvec, universe):
    """The brick MODULE labelling the exchange cover between adjacent pairs ``pair_i`` and
    ``pair_j`` (wall normal ``dimvec``), identified by ISO-CLASS -- not merely by dim-vector.

    A single universe brick of that dim-vector IS the label (the thin case -- kA_n, Kronecker,
    every algebra whose bricks have distinct dim-vectors; byte-identical to the old
    first-match). When several NON-isomorphic bricks share the wall normal, disambiguate by
    the DIRRT brick labelling: the cover ``T > T'`` (T the larger torsion class) is labelled
    by the unique brick ``B`` with ``B in Gen(M_bigger)`` and ``B not in Gen(M_smaller)``.
    The wall normal alone cannot tell P1 from P2 -- their King-stability rays are the two
    opposite halves of the same hyperplane -- but exactly one of them FLIPS its torsion-class
    membership across THIS wall. LOUD ``QuiverlabError`` if the membership does not single out
    exactly one candidate (never a silent wrong wall label / count)."""
    cands = _bricks_with_dimvec(A, dimvec, universe)
    if len(cands) <= 1:
        return cands[0] if cands else None
    si, sj = _gen_size(pair_i), _gen_size(pair_j)
    if si == sj:
        from quiverlab.errors import QuiverlabError
        raise QuiverlabError(
            "brick labelling: adjacent pairs have equal torsion-class size, but a mutation "
            "cover must be a strict inclusion (report the algebra + edge)")
    big, small = (pair_i, pair_j) if si > sj else (pair_j, pair_i)
    Mbig, Msmall = _module_of(big), _module_of(small)
    hits = [B for B in cands if _in_gen(Mbig, B) and not _in_gen(Msmall, B)]
    if len(hits) == 1:
        return hits[0]
    from quiverlab.errors import QuiverlabError
    raise QuiverlabError(
        "brick labelling: %d non-isomorphic bricks share the wall normal %s and "
        "torsion-class membership matched %d of them -- refusing to guess the wall label"
        % (len(cands), _dvtuple(cands[0].dimension_vector()), len(hits)),
        hint="report the algebra + edge; the DIRRT brick label is ambiguous here")


def _dedup_iso(mods):
    """``mods`` deduped by iso-class (dim-vector prefilter, then ``is_isomorphic``)."""
    out = []
    for B in mods:
        if not any(C.dim == B.dim and C.dimension_vector() == B.dimension_vector()
                   and is_isomorphic(C, B) for C in out):
            out.append(B)
    return out


def _same_iso_multiset(s1, s2):
    """True iff the two brick lists are equal as MULTISETS of iso-classes (dim-vector
    prefilter, then ``is_isomorphic`` within a bucket) -- the semibrick equality that keeps
    ``{P1}`` and ``{P2}`` distinct though their dim-vector fingerprints coincide."""
    if len(s1) != len(s2):
        return False
    rem = list(s2)
    for B in s1:
        for k, C in enumerate(rem):
            if (C.dim == B.dim and C.dimension_vector() == B.dimension_vector()
                    and is_isomorphic(C, B)):
                rem.pop(k)
                break
        else:
            return False
    return not rem


def _brick_of_edge(pair_i, pair_j):
    """The brick labelling the cover between two adjacent pairs: ``{"module": B|None,
    "dimvec": {...}, "name": str|None}``. ``B`` is identified by iso-class via
    :func:`_edge_brick` (disambiguated by torsion-class membership on non-thin algebras)."""
    from quiverlab.modules.hom import identify_standard
    from quiverlab.tautilting.mutation import wall_normal
    A = pair_i.algebra
    dv = wall_normal(pair_i, pair_j)
    B = _edge_brick(A, pair_i, pair_j, dv, _torsion_universe(A))
    name = None
    if B is not None:
        std = identify_standard(B)
        if std is not None:
            kind, v = std
            name = {"simple": "S", "projective": "P", "injective": "I"}[kind] + str(v)
    return {"module": B, "dimvec": dv, "name": name}


def bricks(A, budget=512):
    """The bricks labelling the exchange edges of ``A`` (each ``end_dim == 1``), deduped by
    ``is_isomorphic``. Each edge's brick is identified by ISO-CLASS (:func:`_edge_brick`), so
    two non-isomorphic bricks of the same dim-vector (e.g. kZ2/rad^2's P1, P2) are counted
    SEPARATELY -- never collapsed to a first dim-vector match."""
    from quiverlab.tautilting.mutation import exchange_graph
    eg = exchange_graph(A, budget_pairs=budget)
    universe = _torsion_universe(A, budget=budget)
    out = []
    for (i, j) in eg.arrows:
        dv = eg.arrows[(i, j)]["brick"]
        B = _edge_brick(A, eg.vertices[i]["pair"], eg.vertices[j]["pair"], dv, universe)
        if B is None:
            continue
        if not any(X.dim == B.dim and X.dimension_vector() == B.dimension_vector()
                   and is_isomorphic(X, B) for X in out):
            out.append(B)
    return out


def semibricks(A, budget=512):
    """The semibricks of ``A``: for each pair, the set of bricks labelling the edges going
    DOWN out of it (the Asai labelling), deduped by ISO-CLASS. ``#semibricks == #pairs`` in
    the finite case (the four-way identity) -- which REQUIRES the iso-class dedup: two
    semibricks differing only by a same-dim-vector brick (kZ2/rad^2's ``{P1}`` vs ``{P2}``)
    are DISTINCT, though their dim-vector fingerprints coincide. Each returned semibrick is a
    list of pairwise Hom-orthogonal bricks."""
    from quiverlab.tautilting.mutation import exchange_graph
    eg = exchange_graph(A, budget_pairs=budget)
    universe = _torsion_universe(A, budget=budget)
    orient = hasse_orientation(eg)
    down = {i: [] for i in range(len(eg.vertices))}
    for (i, j), d in orient.items():
        a, b = (i, j) if d == "down" else (j, i)      # a -> b downward
        dv = eg.arrows[(i, j)]["brick"]
        B = _edge_brick(A, eg.vertices[i]["pair"], eg.vertices[j]["pair"], dv, universe)
        if B is not None:
            down[a].append(B)
    result = []
    for i in range(len(eg.vertices)):
        sb = _dedup_iso(down[i])
        if not any(_same_iso_multiset(sb, prev) for prev in result):
            result.append(sb)
    return result
