"""Minimal left/right add(M)-approximations (Plan 44 / C7).

A right add(M)-approximation ``f: X -> C`` (``X in add M``) means every ``M -> C``
factors through ``f`` (``Hom_A(M, -)`` surjectivity). It is right MINIMAL iff ``X`` has
no summand mapping to ``0``; the minimal one is unique, with ``X = M^k`` where ``k`` = the
number of generators of ``Hom_A(M, C)`` as a right ``End_A(M)``-module (= ``dim Hom(M,C) /
Hom(M,C)*rad End(M)``, Nakayama). Dually a left add(M)-approximation ``g: C -> M^k`` is a
minimal LEFT ``End_A(M)``-module generating set of ``Hom_A(C, M)`` stacked into a map out
of ``C``. Float-free, exact linear algebra over the Domain (ASS I.2, I.5).

The minimality selector reads ``rad(End_A(M))`` via ``endomorphism._radical_endos`` -- the
trace-form radical, rigorous over ``char 0`` or ``char > dim M`` -- so both constructors
inherit that loud char refusal off scope (never a silently inflated ``k``).
"""
from __future__ import annotations

from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.endomorphism import _radical_endos
from quiverlab.modules.hom import _assert_comparable
from quiverlab.modules.morphism import ModuleHom, direct_sum, hom_basis


def _min_add_generators(M, C):
    """The sublist of ``hom_basis(M, C)`` that minimally generates ``Hom_A(M, C)`` as a
    RIGHT ``End_A(M)``-module: a basis modulo ``Hom(M, C)*rad End(M)`` (the right action
    ``g . rho = rho.then(g)``). Inherits the trace-form char refusal via ``_radical_endos``."""
    homs = hom_basis(M, C)                                   # g: M -> C
    if not homs:
        return homs
    dom = C.domain

    def vec(g):                                             # column-major flatten (C.dim x M.dim)
        return [g.matrix[i][j] for j in range(M.dim) for i in range(C.dim)]

    rad = _radical_endos(M)                                 # list[ModuleHom], may be []
    sub = [vec(rho.then(g)) for g in homs for rho in rad]   # in Hom(M,C)*rad End(M)
    cand = [vec(g) for g in homs]
    keep = lm.independent_modulo(cand, sub, dom)            # minimal (mod-radical) generators
    return [homs[i] for i in keep]


def _min_add_generators_left(M, C):
    """The sublist of ``hom_basis(C, M)`` that minimally generates ``Hom_A(C, M)`` as a
    LEFT ``End_A(M)``-module: a basis modulo ``rad End(M)*Hom(C, M)`` (the left action
    ``rho . h = h.then(rho)``)."""
    homs = hom_basis(C, M)                                   # h: C -> M
    if not homs:
        return homs
    dom = C.domain

    def vec(h):                                            # column-major flatten (M.dim x C.dim)
        return [h.matrix[i][j] for j in range(C.dim) for i in range(M.dim)]

    rad = _radical_endos(M)
    sub = [vec(h.then(rho)) for h in homs for rho in rad]   # rho . h = h.then(rho)
    cand = [vec(h) for h in homs]
    keep = lm.independent_modulo(cand, sub, dom)
    return [homs[i] for i in keep]


def right_add_approximation(M, C):
    """The minimal right ``add(M)``-approximation ``f: M^k ->> C``: every ``M -> C``
    factors through ``f``, ``k`` minimal. ``f`` is ``[g_0 | ... | g_{k-1}]`` on the ``k``
    blocks of ``M^k``, the ``g_i`` a minimal right ``End_A(M)``-module generating set of
    ``Hom_A(M, C)``. ``check=True`` at construction re-certifies ``f`` IS an A-module map."""
    _assert_comparable(M, C, "add-approximation")
    dom = C.domain
    gens = _min_add_generators(M, C)
    k = len(gens)
    if k == 0:                                             # Hom(M, C) = 0: source is 0
        from quiverlab.modules.duality import _zero_module
        Z = _zero_module(M.algebra, side=M.side)
        return ModuleHom(Z, C, lm.zeros(C.dim, 0, dom), check=False)
    D, _incls, _projs = direct_sum(*([M] * k))             # M^k, consecutive M.dim blocks
    fmat = lm.zeros(C.dim, D.dim, dom)
    for i, g in enumerate(gens):
        off = i * M.dim
        for r in range(C.dim):
            gr = g.matrix[r]
            for j in range(M.dim):
                fmat[r][off + j] = gr[j]
    return ModuleHom(D, C, fmat, check=True)


def left_add_approximation(M, C):
    """The minimal left ``add(M)``-approximation ``g: C >-> M^k``: every ``C -> M``
    factors through ``g``, ``k`` minimal. Built DIRECTLY (not by dualizing) so ``g.src``
    IS the passed ``C``: ``g`` stacks the ``h_i: C -> M`` (a minimal left ``End_A(M)``-module
    generating set of ``Hom_A(C, M)``) into the ``k`` blocks of ``M^k``."""
    _assert_comparable(M, C, "add-approximation")
    dom = C.domain
    gens = _min_add_generators_left(M, C)
    k = len(gens)
    if k == 0:                                             # Hom(C, M) = 0: target is 0
        from quiverlab.modules.duality import _zero_module
        Z = _zero_module(M.algebra, side=M.side)
        return ModuleHom(C, Z, lm.zeros(0, C.dim, dom), check=False)
    D, _incls, _projs = direct_sum(*([M] * k))
    gmat = lm.zeros(D.dim, C.dim, dom)
    for i, h in enumerate(gens):
        off = i * M.dim
        for r in range(M.dim):
            hr = h.matrix[r]
            for j in range(C.dim):
                gmat[off + r][j] = hr[j]
    return ModuleHom(C, D, gmat, check=True)
