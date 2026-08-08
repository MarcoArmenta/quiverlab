"""Wall-and-chamber structure via bricks (Plan 63 / R25, Brustle-Smith-Treffinger 2019 +
Demonet-Iyama-Jasso 2019 + King 1994; Kaipel-Treffinger 2023 worked examples). A thin
exact-rational layer over the merged Plan-45 tau-tilting engine -- NO new math engines.

For every brick ``B`` (``end_dim(B) == 1``) the WALL ``D(B)`` is the King theta-semistable
locus computed as an EXACT rational inequality system over the submodule dim-vectors::

    D(B) = {theta : theta . dim B = 0  and  theta . dim N <= 0 for every submodule N <= B}

-- a rational polyhedral cone of codim >= 1 in the stability space
``K0(proj A)_R ~= R^{Q0}`` (BST Def 3.1-3.3, King with (submodule, <= 0)). A SIMPLE brick
gives the full hyperplane ``{theta . dim B = 0}``; a non-simple brick a PROPER face -- e.g.
over kA2, ``D(P1) = {theta1 + theta2 = 0, theta2 <= 0}`` is a RAY (direction ``(1,-1)``), NOT
the full line (the record's headline). The inequality reuses the shipped
:func:`quiverlab.tautilting.stability._submodule_dimvecs` + King convention verbatim, so the
wall is byte-consistent with :func:`~quiverlab.tautilting.stability.is_theta_semistable`.

The CHAMBERS are the g-vector cones of the support tau-tilting pairs (the Plan-45 exchange
graph); the whole structure is a fan (BST 2019): chambers <-> maximal g-cones, walls =
codim-1 boundaries, ``#chambers = #support tau-tilting``. Certified COMPLETE iff ``A`` is
brick-finite <=> tau-tilting-finite (DIJ; decided by the Plan-45 exchange-graph BFS);
otherwise a BOUNDED region with honest truncation (the P62 discipline -- exact as far as it
goes, never claimed complete, no count asserted).

Rigorous over char 0 / char > dim (the Plan-45 brick / is_isomorphic caveat -- QQ default);
off scope the engine inherits the loud ``QuiverlabError`` refusal unchanged. All geometry is
exact ``fractions.Fraction`` -- no floats in ``src/`` (the JS renderer does the only
fraction->pixel conversion). Plan 45 is consumed READ-ONLY and left byte-unchanged."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd

from quiverlab.errors import QuiverlabError

_CITATIONS = ["brustle_smith_treffinger", "demonet_iyama_jasso", "asai_semibricks",
              "king_stability", "kaipel_treffinger"]


@dataclass(frozen=True)
class Wall:
    """The wall ``D(B)`` of a brick ``B`` as an exact inequality system (Plan 63)."""
    brick_dimvec: tuple          # dim B in vertex order (the equality NORMAL): theta.dim B=0
    brick_name: str | None       # "S1"/"P2"/... via identify_standard, else None
    equality: tuple              # == brick_dimvec (the hyperplane {theta . dim B = 0})
    inequalities: tuple          # sorted tuple of dim-vector tuples d with theta.d <= 0,
                                 #   = submodule dim-vectors of B MINUS 0 and the full dim B
    is_full_hyperplane: bool     # True iff `inequalities` empty => D(B) is the whole
                                 #   hyperplane {theta . dim B = 0}
    codim: int                   # 1 (D(B) lies in a codim-1 hyperplane)
    rays: tuple | None           # n <= 2 ONLY: exact extreme rays as vertex-order fraction
                                 #   strings; a full line = two opposite rays, a ray = one;
                                 #   () for n == 1 (D(B) = {0}). None for n >= 3.


def _std_name(B):
    """The standard name ``S_v``/``P_v``/``I_v`` of a brick via ``identify_standard``, else
    None (an unnamed brick is still a valid wall -- shown by dim-vector)."""
    from quiverlab.modules.hom import identify_standard
    std = identify_standard(B)
    if std is None:
        return None
    kind, v = std
    return {"simple": "S", "projective": "P", "injective": "I"}[kind] + str(v)


def _prim2(vec):
    """The primitive integer 2-direction of ``vec`` (divide out the gcd; SIGN kept)."""
    g = gcd(abs(vec[0]), abs(vec[1]))
    if g == 0:
        return (vec[0], vec[1])
    return (vec[0] // g, vec[1] // g)


def _rays_n2(dv, ineqs):
    """The exact extreme rays of ``D(B)`` in ``R^2`` (no floats). The hyperplane
    ``theta . dv = 0`` is the line ``{s * t}`` with direction ``t = (dv[1], -dv[0])`` (the
    normal rotated 90 deg). Each inequality ``theta . d <= 0`` becomes ``s * (t . d) <= 0``,
    a sign constraint on ``s``. No sign constraint => the full line (both +-t); one feasible
    half => a single ray; both halves forced to zero => degenerate (raise -- never for a
    genuine brick)."""
    t = (dv[1], -dv[0])
    need_nonpos = False                               # some d forces s <= 0 (t . d > 0)
    need_nonneg = False                               # some d forces s >= 0 (t . d < 0)
    for d in ineqs:
        c = t[0] * d[0] + t[1] * d[1]
        if c > 0:
            need_nonpos = True
        elif c < 0:
            need_nonneg = True
    if need_nonpos and need_nonneg:
        raise QuiverlabError(
            f"wall_of_brick: the King inequalities force theta = 0 for dim-vector {dv} "
            "(a degenerate wall -- this should not happen for a genuine brick)",
            hint="report the algebra + brick")
    plus = _prim2(t)
    minus = _prim2((-t[0], -t[1]))
    if not need_nonpos and not need_nonneg:
        rays = (plus, minus)                          # full line: two opposite rays
    elif need_nonpos:
        rays = (minus,)                               # s <= 0 feasible => ray -t
    else:
        rays = (plus,)                                # s >= 0 feasible => ray +t
    return tuple((str(r[0]), str(r[1])) for r in rays)


def wall_of_brick(B, *, budget=4096):
    """The wall ``D(B) = {theta : theta . dim B = 0 and theta . dim N <= 0 for every
    submodule N <= B}`` of a brick ``B`` as an exact :class:`Wall` (Plan 63 / R25). Reads the
    shipped exact submodule enumerator (:func:`stability._submodule_dimvecs`), drops the 0
    vector and the full ``dim B`` (both implied by the equality), computes the exact extreme
    rays for ``n <= 2``. Loud ``QuiverlabError`` past ``budget`` (a large brick)."""
    from quiverlab.tautilting.stability import _submodule_dimvecs
    verts = list(B.algebra.quiver.vertices)
    n = len(verts)
    dvdict = B.dimension_vector()
    dv = tuple(int(dvdict[v]) for v in verts)
    subs = _submodule_dimvecs(B, budget=budget)
    ineqs = tuple(sorted(d for d in subs if any(d) and tuple(d) != dv))
    is_full = not ineqs
    name = _std_name(B)
    if n == 2:
        rays = _rays_n2(dv, ineqs)
    elif n == 1:
        # the only stability space is R; theta . (d,) = 0 with d >= 1 forces theta = 0, so
        # D(B) = {0} -- a DEGENERATE point, not the "full hyperplane" (kept honest per plan).
        rays = ()
        is_full = False
    else:
        rays = None                                   # n >= 3: geometry via grouped facets
    return Wall(brick_dimvec=dv, brick_name=name, equality=dv, inequalities=ineqs,
                is_full_hyperplane=is_full, codim=1, rays=rays)


# --------------------------------------------------------------------------- #
# the full wall-and-chamber structure payload
# --------------------------------------------------------------------------- #
def _build_chambers(eg, verts):
    """One chamber per support tau-tilting pair (a maximal g-cone): the g-matrix, its
    g-vector columns as exact rays, the pair label/support, and (n == 3) the L1/octahedron
    projection (mirroring the Plan-45 fan)."""
    from quiverlab.tautilting.stability import _l1_project
    n = len(verts)
    chambers = []
    for i, rec in enumerate(eg.vertices):
        G = rec["g_matrix"]
        rays = [[G[r][c] for r in range(n)] for c in range(n)]     # columns = g-vectors
        ch = {"id": i, "g_matrix": G,
              "rays": [[str(Fraction(x)) for x in ray] for ray in rays],
              "label": rec["label"], "support": list(rec["support"]),
              "is_initial": rec["is_initial"]}
        if n == 3:
            proj = [_l1_project(ray) for ray in rays]
            ch["rays_l1"] = [p[0] for p in proj]
            ch["faces"] = [p[1] for p in proj]
            ch["net2d"] = [p[2] for p in proj]
        chambers.append(ch)
    return chambers


def _shared_gvectors(eg, i, j):
    """The g-vectors IN the wall between adjacent pairs ``i, j`` (their SHARED g-columns);
    each is a facet vector on the King hyperplane ``theta . dim B = 0`` (Plan-45 fan idiom)."""
    pi = eg.vertices[i]["pair"]
    pj = eg.vertices[j]["pair"]
    return [list(c) for c in sorted(pi.g_key() & pj.g_key())]


def _iso(X, Y):
    """A cheap-prefiltered ``is_isomorphic`` (dim + dim-vector before the Hom certificate)."""
    from quiverlab.modules.hom import is_isomorphic
    return (X.dim == Y.dim and X.dimension_vector() == Y.dimension_vector()
            and is_isomorphic(X, Y))


def _wall_dict(wall, wid, facets, verts, eg, n):
    """A JSON-ready wall record from a :class:`Wall` + its grouped exchange-edge facets.
    For ``n == 3`` the drawing rays are the DISTINCT shared g-vectors of the grouped edges
    (each lies IN D(B) -- self-certified by the caller), L1-projected."""
    from quiverlab.tautilting.stability import _l1_project
    d = {"id": wid,
         "brick_dimvec": {verts[k]: int(wall.brick_dimvec[k]) for k in range(n)},
         "brick_name": wall.brick_name,
         "equality": [int(x) for x in wall.equality],
         "inequalities": [list(map(int, ineq)) for ineq in wall.inequalities],
         "is_full_hyperplane": wall.is_full_hyperplane,
         "codim": wall.codim,
         "facets": [list(f) for f in facets]}
    if n <= 2:
        d["rays"] = None if wall.rays is None else [list(r) for r in wall.rays]
    else:   # n == 3: rays = distinct shared g-vectors of the grouped facet edges
        seen = set()
        rays = []
        for (i, j) in facets:
            for g in _shared_gvectors(eg, i, j):
                key = tuple(g)
                if key not in seen:
                    seen.add(key)
                    rays.append(g)
        d["rays"] = [[str(Fraction(x)) for x in g] for g in rays]
        proj = [_l1_project(g) for g in rays]
        d["rays_l1"] = [p[0] for p in proj]
        d["faces"] = [p[1] for p in proj]
        d["net2d"] = [p[2] for p in proj]
    return d


def _build_walls_complete(A, eg, brs, verts, n, budget):
    """One wall per canonical brick ISO-CLASS (``torsion.bricks(A)``). Each edge is assigned
    to the wall of the brick it is ISO to -- NEVER grouped by dim-vector (BST Rem 3.19:
    kZ2/rad^2 carries two non-isomorphic (1,1)-bricks P1, P2 on opposite half-rays of one
    hyperplane; a dim-vector key would merge them into 3 walls, the truth is 4). The shipped
    :func:`torsion._edge_brick` disambiguates same-dim-vector bricks by torsion-class
    membership -- reused here so facets <-> edges stays a bijection."""
    from quiverlab.tautilting.torsion import _edge_brick, _torsion_universe
    universe = _torsion_universe(A, budget=budget)
    edge_bricks = {}
    for (i, j) in eg.arrows:
        dv = eg.arrows[(i, j)]["brick"]
        edge_bricks[(i, j)] = _edge_brick(
            A, eg.vertices[i]["pair"], eg.vertices[j]["pair"], dv, universe)
    walls = []
    for wid, B in enumerate(brs):
        wall = wall_of_brick(B, budget=budget)
        facets = sorted((i, j) for (i, j), EB in edge_bricks.items()
                        if EB is not None and _iso(EB, B))
        walls.append(_wall_dict(wall, wid, facets, verts, eg, n))
    return walls


def _build_walls_bounded(A, eg, verts, n):
    """The honest BOUNDED-REGION walls for a tau-tilting-INFINITE (budget-capped) algebra:
    the discovered exchange-edge facets grouped by their King wall-normal (the hyperplane),
    each flagged ``partial`` -- the brick module + full D(B) inequality system are NOT
    enumerated (that needs the COMPLETE torsion universe, which does not exist here). No
    count is claimed (DIJ: brick-finite <=> tau-tilting-finite). Cheap by design: only the
    one exchange-graph BFS is spent, never the per-brick submodule/iso work."""
    groups = {}
    order = []
    for (i, j) in eg.arrows:
        dv = eg.arrows[(i, j)]["brick"]               # the King wall normal (free from BFS)
        key = tuple(int(dv[v]) if dv else 0 for v in verts) if dv else None
        if key is None:
            continue
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append((i, j))
    walls = []
    for wid, key in enumerate(order):
        walls.append({
            "id": wid,
            "brick_dimvec": {verts[k]: key[k] for k in range(n)},
            "brick_name": None,
            "equality": list(key),
            "inequalities": [],
            "is_full_hyperplane": False,
            "codim": 1,
            "rays": None,
            "facets": [list(f) for f in sorted(groups[key])],
            "partial": True,
        })
    return walls


_TRUNC = ("A appears to be tau-tilting-infinite / brick-infinite (the exchange-graph BFS "
          "did not close within budget {N}); the region shown is the sub-fan explored from "
          "(A,0) -- each chamber and wall is exact, but the structure is NOT complete and no "
          "count is claimed (DIJ: brick-finite <=> tau-tilting-finite).")


def wall_chamber_structure(A, *, budget=512):
    """The full wall-and-chamber payload of ``A`` via bricks (Plan 63 / R25). Chambers =
    g-vector cones of the support tau-tilting pairs (the Plan-45 exchange graph); walls = one
    per brick ISO-CLASS, each carrying its exact ``D(B)`` inequality system + (rank <= 3)
    drawing rays; the chamber<->wall adjacency (exchange edges grouped by brick); the four
    counts; and ``render in {fan2d, fan3d, table}``. Certified COMPLETE iff ``A`` is
    brick-finite <=> tau-tilting-finite (DIJ, decided by the BFS closing); otherwise a BOUNDED
    region with ``complete=False``, ``status="budget"``, a ``truncation`` note, and NO counts.
    The brick / is_isomorphic char caveat propagates loudly (char 0 / char > dim; QQ default).
    """
    from quiverlab.tautilting.mutation import exchange_graph
    verts = list(A.quiver.vertices)
    n = len(verts)
    render = "fan2d" if n == 2 else ("fan3d" if n == 3 else "table")
    eg = exchange_graph(A, budget_pairs=budget)
    out = {
        "kind": "wall_chamber",
        "n": n,
        "complete": eg.is_complete,
        "status": eg.status,
        "render": render,
        "references": list(_CITATIONS),
    }
    chambers = _build_chambers(eg, verts)
    out["chambers"] = chambers
    out["num_chambers"] = len(chambers)
    if not eg.is_complete:
        walls = _build_walls_bounded(A, eg, verts, n)
        out["walls"] = walls
        out["num_walls"] = len(walls)
        out["counts"] = None
        out["green_count"] = None
        out["truncation"] = _TRUNC.format(N=budget)
        return out
    from quiverlab.tautilting.green import maximal_green_sequences
    from quiverlab.tautilting.torsion import bricks as torsion_bricks
    brs = torsion_bricks(A, budget=budget)
    walls = _build_walls_complete(A, eg, brs, verts, n, budget)
    out["walls"] = walls
    out["num_walls"] = len(walls)
    out["counts"] = {"chambers": len(eg.vertices), "walls": len(walls),
                     "bricks": len(brs), "s_tau_tilt": len(eg.vertices)}
    out["green_count"] = maximal_green_sequences(A, cap=budget)["count"]
    out["truncation"] = None
    return out


def _is_green_path(eg, orient, seq):
    """True iff ``seq`` (a list of pair-ids) is a monotone DOWNWARD chamber path from the
    source ``(A,0)`` to the sink ``(0,A)`` -- a maximal green sequence read as a green path
    through the wall-and-chamber structure (BST/Plan-45). Falsifiable: rejects a path that
    does not start at ``(A,0)``, does not end at ``(0,A)``, or crosses any wall upward."""
    if not seq:
        return False
    verts = list(eg.vertices[seq[0]]["pair"].algebra.quiver.vertices)
    if not eg.vertices[seq[0]]["is_initial"]:
        return False
    last = eg.vertices[seq[-1]]
    if last["summand_dimvecs"] or set(last["support"]) != set(verts):
        return False                                  # the sink must be the terminal (0, A)
    for t in range(len(seq) - 1):
        a, b = seq[t], seq[t + 1]
        e = (min(a, b), max(a, b))
        if e not in eg.arrows:
            return False                              # not an exchange edge
        # hasse_orientation keys on (min,max): "down" => min is the higher (downward source);
        # "up" => max is. The green step must go DOWNWARD out of seq[t].
        source = e[0] if orient[e] == "down" else e[1]
        if source != a:
            return False                              # the step goes UP, not down
    return True

