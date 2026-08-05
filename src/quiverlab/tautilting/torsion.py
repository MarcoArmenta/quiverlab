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

from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.hom import end_dim, is_isomorphic
from quiverlab.modules.morphism import direct_sum


# --------------------------------------------------------------------------- #
# the module universe (indecomposable tau-rigid summands across the exchange graph)
# --------------------------------------------------------------------------- #
_UNIVERSE_CACHE = {}


def _torsion_universe(A, budget=512):
    """The indecomposable modules appearing as summands across every support tau-tilting
    pair of ``A`` (a sufficient universe to fingerprint the torsion classes -- AIR: the
    Ext-projectives of a class generate it). Deduped by ``is_isomorphic``; cached per
    algebra. Raises loudly (via the BFS) if ``A`` is tau-tilting-infinite."""
    key = id(A)
    cached = _UNIVERSE_CACHE.get(key)
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
    _UNIVERSE_CACHE[key] = uni
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
def _brick_module(A, dimvec, universe):
    """A brick (``end_dim == 1``) with the given dim-vector, found in ``universe``, or
    None. (For the tested tau-tilting-finite algebras every brick is a pair summand.)"""
    if dimvec is None:
        return None
    target = {v: int(dimvec.get(v, 0)) for v in A.quiver.vertices}
    for X in universe:
        if X.dimension_vector() == target and end_dim(X) == 1:
            return X
    return None


def _brick_of_edge(pair_i, pair_j):
    """The brick labelling the cover between two adjacent pairs: ``{"module": B|None,
    "dimvec": {...}, "name": str|None}``. ``B`` is the brick (``end_dim == 1``) with the
    wall-normal dim-vector, found in the module universe."""
    from quiverlab.modules.hom import identify_standard
    from quiverlab.tautilting.mutation import wall_normal
    A = pair_i.algebra
    dv = wall_normal(pair_i, pair_j)
    B = _brick_module(A, dv, _torsion_universe(A))
    name = None
    if B is not None:
        std = identify_standard(B)
        if std is not None:
            kind, v = std
            name = {"simple": "S", "projective": "P", "injective": "I"}[kind] + str(v)
    return {"module": B, "dimvec": dv, "name": name}


def bricks(A, budget=512):
    """The bricks labelling the exchange edges of ``A`` (each ``end_dim == 1``), deduped by
    ``is_isomorphic``."""
    from quiverlab.tautilting.mutation import exchange_graph
    eg = exchange_graph(A, budget_pairs=budget)
    universe = _torsion_universe(A, budget=budget)
    out = []
    for (i, j) in eg.arrows:
        dv = eg.arrows[(i, j)]["brick"]
        B = _brick_module(A, dv, universe)
        if B is None:
            continue
        if not any(X.dim == B.dim and X.dimension_vector() == B.dimension_vector()
                   and is_isomorphic(X, B) for X in out):
            out.append(B)
    return out


def semibricks(A, budget=512):
    """The semibricks of ``A``: for each pair, the set of bricks labelling the edges going
    DOWN out of it (the Asai labelling), deduped. ``#semibricks == #pairs`` in the finite
    case (the four-way identity). Each returned semibrick is a list of pairwise
    Hom-orthogonal bricks."""
    from quiverlab.tautilting.mutation import exchange_graph
    eg = exchange_graph(A, budget_pairs=budget)
    universe = _torsion_universe(A, budget=budget)
    orient = hasse_orientation(eg)
    down = {i: [] for i in range(len(eg.vertices))}
    for (i, j), d in orient.items():
        a, b = (i, j) if d == "down" else (j, i)      # a -> b downward
        dv = eg.arrows[(i, j)]["brick"]
        B = _brick_module(A, dv, universe)
        if B is not None:
            down[a].append(B)
    seen = {}
    for i in range(len(eg.vertices)):
        key = frozenset(_dvtuple(B.dimension_vector()) for B in down[i])
        if key not in seen:
            seen[key] = down[i]
    return list(seen.values())
