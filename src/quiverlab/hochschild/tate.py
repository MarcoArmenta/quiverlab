"""Tate-Hochschild (singular Hochschild) cohomology -- the public surface (Plan 76 / R3).

`HHhat^m(A) = Extbar^m_{A^e}(A, A)` is defined for EVERY integer m, positive and
negative, as the cohomology of `Hom_{A^e}(T, A)` for a complete resolution `T` of the
diagonal bimodule (`engine/complete_resolution.py` builds `T`).  Above the Gorenstein
dimension `d` of `A^e` it coincides with ordinary Hochschild cohomology
(Bergh-Jorgensen, JNCG 7 (2013)):

    HHhat^m(A) = HH^m(A)   for every m >= d + 1 =: agrees_from,

and `d = 0` when `A` is self-injective, so agreement starts in degree 1.  **Degree 0 is
NOT ordinary**: `HHhat^0` is the stable centre, a proper quotient of `Z(A) = HH^0`
(k[x]/(x^n) over a field of characteristic not dividing n: `HHhat^0 = n-1` against
`HH^0 = n`).  The negative degrees are the genuinely new object.

THREE ROUTES, each honest about what it does not know:

  * `native` -- the GF(p) complete-resolution splice.  Knows EVERY degree, including
    `HHhat^0` and all negatives.  Scope: **self-injective** `A` over `GF(p)` (then
    `A^e` is self-injective, the `A^e`-dual of a projective is projective and the
    splice is well defined) -- periodic or not.
  * `duality` -- the Bergh-Jorgensen duality for SYMMETRIC `A` over any exact Domain,
    `dim HHhat^{-j} = dim HHhat^{j-1}`, composed with the threshold.  It therefore
    supplies `m >= 1` (`= HH^m`) and `m <= -2` (`= HH^{|m|-1}`).  Degrees `0` and `-1`
    are a CLOSED 2-cycle -- the duality links them to each other and to no
    threshold-known positive -- so both are returned `None`, never `HH^0`.
  * `positive` -- the threshold alone, any Gorenstein `A` over any Domain: `m >=
    agrees_from` only, `neg_dims = None`, degree 0 `None`.  A truthful partial answer.

`engine="auto"` takes `native` when it applies, else `duality` for symmetric input,
else `positive`.  An explicit `engine=` refuses loudly when its hypothesis fails.

DEFERRED (loud): a Gorenstein but NOT self-injective algebra.  `D_{A^e}(P_n)` is
projective only when `A^e` is self-injective, so the dual-splice cannot build its
complete resolution; the construction that can is a periodicity-extension splice
(Usui), the named follow-up.  The eventual-periodicity certificate `tate_periodicity`
is computed for such an algebra regardless -- it is a property of the algebra, not of
the route.
"""
from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.fields.primefield import PrimeField

__all__ = ["TateHochschild", "TatePeriodicity", "nakayama_permutation",
           "tate_hochschild", "tate_hochschild_block", "tate_periodicity"]

_REFERENCES = ("bergh_jorgensen_tate", "wang_singular_hh", "keller_singular_hh",
               "usui_tate_periodic")

_ENGINES = ("auto", "native", "duality", "positive")


@dataclass(frozen=True)
class TateHochschild:
    """A frozen Tate-Hochschild report.  A data record, not a predicate (no __bool__).

    `pos_dims[m]` is `dim HHhat^m` for `0 <= m <= top`, or `None` where the ROUTE does
    not know it -- `pos_dims[0]` is filled on `native` only, because `HHhat^0` is not
    `HH^0`.  `neg_dims[j-1]` is `dim HHhat^{-j}` for `1 <= j <= top`, `None` entrywise
    where unknown, and `None` as a whole on the `positive` route.  `ordinary_pos` is
    always the ORDINARY `HH^0..HH^top` -- the agreement anchor; its degree 0 is `HH^0`
    and must never be read as `HHhat^0`.
    """
    top: int
    pos_dims: tuple
    neg_dims: tuple | None
    hat_hh0: int | None
    agrees_from: int | None
    ordinary_pos: tuple
    period: int | None
    periodicity_degree: int | None
    scope: str
    engine: str
    cup: dict | None
    references: tuple
    note: str | None


@dataclass(frozen=True)
class TatePeriodicity:
    """Eventual periodicity of `A` over `A^e` (Usui's hypothesis for the invertible
    homogeneous element of the Tate ring)."""
    period: int | None
    eventually_periodic: bool
    invertible_element_degree: int | None
    complexity: int | None
    note: str


def nakayama_permutation(A):
    """{vertex: vertex} with `soc(e_v A) = S_{pi(v)}` -- re-exported from the engine."""
    from quiverlab.engine.complete_resolution import nakayama_permutation as _np
    return _np(A)


# ---------------------------------------------------------------------------
# scope gates
# ---------------------------------------------------------------------------
def _is_selfinjective(A):
    try:
        return bool(A.is_selfinjective())
    except QuiverlabError:
        return False


def _is_symmetric(A):
    try:
        return bool(A.is_symmetric())
    except QuiverlabError:
        return False


def _require_gorenstein(A):
    """`HHhat^*` exists only when `A^e` is Gorenstein.  Self-injective guarantees it;
    otherwise ask for the Gorenstein certificate and refuse loudly without one."""
    if _is_selfinjective(A):
        return True
    try:
        gor = bool(A.is_gorenstein())
    except QuiverlabError:
        gor = False
    if not gor:
        raise QuiverlabError(
            "Tate-Hochschild cohomology needs a complete resolution over A^e, which "
            "requires A^e to be (Iwanaga-)Gorenstein; this algebra is neither "
            "self-injective nor certified Gorenstein",
            hint="the eventual-periodicity certificate tate_periodicity(A) is still "
                 "available -- it is a property of the algebra, not of a route")
    return False


def _deferred_refusal():
    return QuiverlabError(
        "Tate-Hochschild cohomology of a Gorenstein but NOT self-injective algebra is "
        "DEFERRED: the D_{A^e}(P_n)-dual splice needs A^e self-injective for the dual "
        "half to consist of projectives, which fails here",
        hint="the periodicity-extension complete resolution (Usui) is the named "
             "follow-up task; tate_periodicity(A) is computed regardless")


def _prime_of(A):
    dom = A.domain
    return dom.p if isinstance(dom, PrimeField) else None


def _choose_engine(A, engine, selfinj):
    if engine not in _ENGINES:
        raise QuiverlabError(
            f"unknown Tate-Hochschild engine {engine!r}; valid engines are "
            f"{', '.join(_ENGINES)}")
    p = _prime_of(A)
    if engine == "auto":
        if selfinj and p is not None:
            return "native"
        if _is_symmetric(A):
            return "duality"
        return "positive"
    if engine == "native":
        if p is None:
            raise QuiverlabError(
                "the native Tate-Hochschild route is the GF(p) complete-resolution "
                f"splice; this algebra is over {A.domain.name}",
                hint="use engine='duality' (symmetric, any exact Domain) or "
                     "engine='positive' (the threshold alone), or build the algebra "
                     "over GF(p)")
        if not selfinj:
            raise _deferred_refusal()
    if engine == "duality" and not _is_symmetric(A):
        raise QuiverlabError(
            "the duality route needs a SYMMETRIC algebra (nu^2 = id), which is what "
            "makes dim HHhat^{-j} = dim HHhat^{j-1}; this algebra is not symmetric",
            hint="use engine='native' over GF(p) for the full ring, or "
                 "engine='positive' for the threshold degrees only")
    return engine


# ---------------------------------------------------------------------------
# the public surface
# ---------------------------------------------------------------------------
def tate_hochschild(A, top, *, engine="auto", max_cells=4_000_000):
    """`HHhat^m(A)` for `-top <= m <= top`.  See the module docstring for the routes."""
    top = int(top)
    if top < 0:
        raise QuiverlabError("the Tate-Hochschild top degree must be >= 0")
    selfinj = _require_gorenstein(A)
    route = _choose_engine(A, engine, selfinj)

    ordinary = tuple(A.hochschild_cohomology(top, max_cells=max_cells).dims)
    agrees_from = _agrees_from(A, selfinj)
    scope = "self_injective" if selfinj else "positive_only"

    if route == "native":
        return _native(A, top, ordinary, agrees_from, scope, max_cells)
    if route == "duality":
        return _duality(A, top, ordinary, agrees_from, scope)
    return _positive(A, top, ordinary, agrees_from, scope)


def _agrees_from(A, selfinj):
    """`d + 1` with `d` the Gorenstein dimension of `A^e`.

    Self-injective `A` makes `A^e` self-injective, so `d = 0` and the threshold is 1 --
    a theorem, not a computation.  Off that scope the threshold would need
    `gorenstein_dimension(A^e)`; rather than guess it we return `None` and say so in
    the note, and the `positive` route then claims no agreement degree at all."""
    return 1 if selfinj else None


def _native(A, top, ordinary, agrees_from, scope, max_cells):
    from quiverlab.engine.complete_resolution import tate_cohomology_dims
    p = _prime_of(A)
    dims = tate_cohomology_dims(A, top, p)
    pos = tuple(int(dims[m]) for m in range(top + 1))
    neg = tuple(int(dims[-j]) for j in range(1, top + 1))
    cert = tate_periodicity(A, p=p)
    note = ("the native route builds the complete resolution, so every degree -- "
            "including HHhat^0, which is the stable centre and NOT HH^0 -- is "
            "engine-computed")
    return TateHochschild(
        top=top, pos_dims=pos, neg_dims=neg, hat_hh0=pos[0] if pos else None,
        agrees_from=agrees_from, ordinary_pos=ordinary, period=cert.period,
        periodicity_degree=cert.invertible_element_degree, scope=scope,
        engine="native", cup=None, references=_REFERENCES, note=note)


def _duality(A, top, ordinary, agrees_from, scope):
    """Symmetric: `HHhat^m = HH^m` for `m >= 1` and `HHhat^{-j} = HH^{j-1}` for `j >= 2`.

    Degrees 0 and -1 form a closed 2-cycle under the reflection `n <-> -(n+1)` -- the
    duality relates them to EACH OTHER and to nothing the threshold knows -- so both
    stay `None`.  Returning `HH^0` there would be wrong, not merely imprecise."""
    pos = (None,) + tuple(ordinary[1:top + 1])
    neg = (None,) + tuple(ordinary[j - 1] for j in range(2, top + 1))
    note = ("the duality route supplies degrees m >= 1 (= HH^m) and m <= -2 "
            "(= HH^{|m|-1}); degrees 0 and -1 are a closed 2-cycle it cannot reach, so "
            "they are absent -- HHhat^0 is NOT HH^0 and is native-only")
    return TateHochschild(
        top=top, pos_dims=pos, neg_dims=neg, hat_hh0=None, agrees_from=agrees_from,
        ordinary_pos=ordinary, period=None, periodicity_degree=None, scope=scope,
        engine="duality", cup=None, references=_REFERENCES, note=note)


def _positive(A, top, ordinary, agrees_from, scope):
    lo = agrees_from if agrees_from is not None else None
    if lo is None:
        pos = tuple([None] * (top + 1))
        note = ("the Gorenstein dimension of A^e was not certified, so no agreement "
                "degree can be named and no Tate value is claimed; only the ordinary "
                "HH^* anchor is reported")
    else:
        pos = tuple([None] * min(lo, top + 1)
                    + [ordinary[m] for m in range(lo, top + 1)])
        note = (f"the positive route reports only degrees m >= {lo}, where the "
                f"Bergh-Jorgensen threshold makes HHhat^m = HH^m; degree 0 is absent "
                f"(HHhat^0 is the stable centre, NOT HH^0) and no negative degree is "
                f"claimed")
    return TateHochschild(
        top=top, pos_dims=pos, neg_dims=None, hat_hh0=None, agrees_from=agrees_from,
        ordinary_pos=ordinary, period=None, periodicity_degree=None, scope=scope,
        engine="positive", cup=None, references=_REFERENCES, note=note)


# ---------------------------------------------------------------------------
# eventual periodicity
# ---------------------------------------------------------------------------
def tate_periodicity(A, *, max_period=12, p=32003, bound=None):
    """Is `A` eventually periodic over `A^e`, and with what period?

    The cheap necessary signal is the `A^e`-Betti growth: complexity 1 means bounded
    ranks (a periodic candidate), complexity >= 2 rules periodicity out.  The period
    itself is read off the minimal `A^e`-resolution's ranks and then CERTIFIED by
    comparing the syzygy differentials, never guessed from the ranks alone; when no
    period is confirmed within `max_period` the answer is an honest `None`.
    """
    prime = _prime_of(A) or int(p)
    try:
        cx = int(A.complexity(max(2 * max_period, 8)))
    except QuiverlabError:
        cx = None
    if cx is not None and cx != 1:
        return TatePeriodicity(
            period=None, eventually_periodic=False, invertible_element_degree=None,
            complexity=cx,
            note=("the A^e-Betti numbers grow (complexity %s), so A is not eventually "
                  "periodic over A^e; the Tate ring is still served as a bounded "
                  "window" % cx))
    period = _period_from_resolution(A, max_period, prime)
    if period is None:
        return TatePeriodicity(
            period=None, eventually_periodic=False, invertible_element_degree=None,
            complexity=cx,
            note=("no period <= %d was confirmed on the minimal A^e-resolution; the "
                  "answer is absent, not negative" % max_period))
    return TatePeriodicity(
        period=period, eventually_periodic=True, invertible_element_degree=period,
        complexity=cx,
        note=("Omega^{n+%d}_{A^e}(A) = Omega^n_{A^e}(A) on the minimal A^e-resolution; "
              "for Gorenstein A this is Usui's criterion for an invertible homogeneous "
              "element of degree %d in the Tate-Hochschild ring" % (period, period)))


def _period_from_resolution(A, max_period, p):
    """Least `q <= max_period` with the minimal `A^e`-resolution `q`-periodic.

    Certified on the DIFFERENTIALS (the syzygy data), not on the Betti numbers: equal
    ranks are necessary and famously not sufficient."""
    import numpy as np

    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.complete_resolution import _vertex_labels
    from quiverlab.engine.resolutions_minimal import minimal_resolution
    labels = _vertex_labels(A)
    Aq = A if len(labels) > 1 else A.unit_adapted()
    depth = 2 * max_period + 2
    rks, cols, eng, truncated = minimal_resolution(to_engine(Aq), depth, p)
    if truncated is not None:
        return None
    tags = eng.corner_tags

    def same(n, k):
        if rks.get(n, 0) != rks.get(k, 0) or rks.get(n - 1, 0) != rks.get(k - 1, 0):
            return False
        if tags is not None and (tags.get(n) != tags.get(k)
                                 or tags.get(n - 1) != tags.get(k - 1)):
            return False
        a, b = cols.get(n) or [], cols.get(k) or []
        return len(a) == len(b) and all(
            np.array_equal(np.asarray(x) % p, np.asarray(y) % p) for x, y in zip(a, b))

    for q in range(1, max_period + 1):
        start = 1
        if all(same(n, n + q) for n in range(start, start + q + 1)):
            return q
    return None


# ---------------------------------------------------------------------------
# runner block
# ---------------------------------------------------------------------------
def tate_hochschild_block(A, top, *, engine="auto", max_cells=4_000_000):
    """The JSON block both runners serve (server tier and the Pyodide twin)."""
    r = tate_hochschild(A, top, engine=engine, max_cells=max_cells)
    return {
        "kind": "tate_hochschild",
        "top": r.top,
        "pos_dims": list(r.pos_dims),
        "neg_dims": None if r.neg_dims is None else list(r.neg_dims),
        "hat_hh0": r.hat_hh0,
        "agrees_from": r.agrees_from,
        "ordinary_pos": list(r.ordinary_pos),
        "period": r.period,
        "periodicity_degree": r.periodicity_degree,
        "scope": r.scope,
        "engine": r.engine,
        "references": list(r.references),
        "note": r.note,
    }
