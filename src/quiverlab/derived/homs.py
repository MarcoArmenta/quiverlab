"""Reified hyper-Hom: a basis of H^n(Hom^.(X, Y)) as chain maps X -> Y[n]
(Plan 43 / derived category). Thin accessor over P39's private Hom-total-complex
internals (``_hom_total_blocks`` / ``_delta_total`` / ``_combine_homs``) -- never
recompute the Hom complex a second way (the P41 basis-mismatch discipline).

Why a cocycle IS a chain map ``X -> Y[n]`` (Weibel 2.7.4, convention (*) in
``complexes.py``): ``Hom^n(X, Y) = (+)_p Hom(X_p, Y_{p-n})``; a component
``f_p: X_p -> Y_{p-n}`` lands in ``Y.shift(n).term(p) = Y.term(p-n)``. The cocycle
condition ``d^Y f_p = (-1)^n f_{p-1} d^X`` is exactly the chain-map square for
``X -> Y[n]`` (the shifted complex ``Y[n]`` carries ``(-1)^n d^Y``). So
``ChainMap(X, Y.shift(n), comps, check=True)`` PASSES on any cocycle and FAILS on a
non-cocycle -- the reification is self-certifying.

Canonical (byte-reproducible) coset representatives: each reified class is a
deterministic cocycle reduced modulo the coboundary span (see ``_reduce_mod_span``),
so the same class is produced regardless of the incidental cocycle-basis order --
the CS/Plan-17 free-variables-zero canonicalisation, adapted to a column span.
"""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.fields import linalg
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.complexes import (ChainComplex, ChainMap, hyper_hom_dims,
                                         _hom_total_blocks, _delta_total,
                                         _combine_homs)


def _reduce_mod_span(r, span_cols, dom):
    """Canonical representative of the coset ``r + span(span_cols)``: reduce ``r`` by
    the RREF of the span (the span vectors stacked as ROWS -- row operations preserve
    the row space, so every RREF row is itself a span vector) until ``r`` has a zero
    coordinate at every pivot position of the span. Only span vectors are ever
    subtracted, so a cocycle stays a cocycle and its homology class is unchanged (the
    span here is the coboundary span, and a coboundary is a cocycle).

    NB deviation from the plan's ``reduce_mod_nullspace(r, transpose(B))``: that reduces
    modulo the *left*-nullspace of the coboundary matrix ``B`` (the orthogonal
    complement of the coboundary span), which would move ``r`` OUT of ``ker delta^n``
    and fail the ChainMap self-cert. This reduces modulo the span itself -- correct AND
    canonical."""
    if not span_cols:
        return list(r)
    R, piv = linalg.rref([list(c) for c in span_cols], dom)   # rows span the coboundaries
    y = list(r)
    for row, pc in zip(R, piv):
        c = y[pc]
        if dom.is_zero(c):
            continue
        y = [dom.sub(a, dom.mul(c, b)) for a, b in zip(y, row)]
    return y


def _coset_reps(cocycle_cols, coboundary_cols, dom):
    """Canonical representatives of ``ker/im``: a maximal subset of ``cocycle_cols``
    independent modulo ``coboundary_cols`` (its size = ``dim H^n``), each reduced
    modulo the coboundary span for byte-stability."""
    idx = lm.independent_modulo(cocycle_cols, coboundary_cols, dom)
    reps = [list(cocycle_cols[i]) for i in idx]
    if coboundary_cols:
        reps = [_reduce_mod_span(r, coboundary_cols, dom) for r in reps]
    return reps


def hyper_hom_basis(X, Y, n):
    """A basis of ``H^n(Hom^.(X, Y))`` reified as chain maps ``X -> Y.shift(n)``. For
    ``X`` PERFECT this is ``Hom_{D^b(mod A)}(X, Y[n])`` (P39 header). Each returned map
    is built from a canonical coset representative of
    ``ker(delta^n) / im(delta^{n-1})`` and is a genuine chain map
    (``ChainMap(..., check=True)`` is the reification self-cert). SELF-CERT:
    ``len(result) == hyper_hom_dims(X, Y, n, n)[n]`` (asserted here). Raises loudly if
    ``X`` is not certified perfect (use ``projective_model(X, ...)`` first)."""
    if not X.is_perfect():
        raise QuiverlabError(
            "hyper_hom_basis: X must be a perfect complex; resolve it with "
            "projective_model(X, ...) first (P39).")
    dom = X.domain
    blocks, cdim = _hom_total_blocks(X, Y, n, dom)
    dn, _s, _t = _delta_total(X, Y, n, dom)             # delta^n : Hom^n -> Hom^{n+1}
    dn1, _s1, _t1 = _delta_total(X, Y, n - 1, dom)      # delta^{n-1}: Hom^{n-1}->Hom^n
    if dn and dn[0]:
        cocycles = lm.kernel_columns(dn, dom)
    else:                                               # delta^n = 0 => every cochain
        ident = lm.identity(cdim, dom)
        cocycles = [lm.col(ident, j) for j in range(cdim)]
    cobounds = ([lm.col(dn1, j) for j in range(len(dn1[0]))]
                if (dn1 and dn1[0]) else [])
    reps = _coset_reps(cocycles, cobounds, dom)
    # Y.shift(0) is a structural copy of Y; use Y itself at n == 0 so the reified
    # degree-0 classes share Y as their target object -- ChainMap.then (which matches
    # the middle complex by identity) then composes End(T) classes (Task 3).
    Yn = Y if n == 0 else Y.shift(n)
    maps = []
    for rep in reps:
        comps = {}
        for b in blocks:
            coeffs = rep[b["offset"]: b["offset"] + b["count"]]
            comps[b["p"]] = _combine_homs(b["homs"], coeffs,
                                          X.term(b["p"]).dim, Y.term(b["q"]).dim, dom)
        maps.append(ChainMap(X, Yn, comps, check=True))  # cocycle => valid chain map
    # self-cert: the homotopy-quotient dimension equals P39's rank-formula dimension.
    if len(maps) != hyper_hom_dims(X, Y, n, n)[n]:
        raise QuiverlabError(
            "hyper_hom_basis: reified class count != hyper_hom_dims (coset/rank "
            "mismatch -- basis reification bug)")
    return maps
