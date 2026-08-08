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
        # the only stability space is R; theta . (d,) = 0 with d >= 1 forces theta = 0.
        rays = ()
    else:
        rays = None                                   # n >= 3: geometry via grouped facets
    return Wall(brick_dimvec=dv, brick_name=name, equality=dv, inequalities=ineqs,
                is_full_hyperplane=is_full, codim=1, rays=rays)
