# SPDX-License-Identifier: MIT
"""Generalized Koszulity recognizers beyond the quadratic case (Plan 77 / R36).

Plan 27 owns QUADRATIC Koszulity (``modules/koszul.py``: ``is_quadratic``,
``g_quadratic_certificate`` = Priddy PBW, ``quadratic_dual``, ``froberg_obstruction``,
and the three-valued ``YonedaPresentation.koszul``).  This module owns everything
past it, as a thin recognizer layer over the SHIPPED Plan-05 minimal projective
resolutions of the simples and the Plan-27 Yoneda engine -- no new homological
engine is built here:

* :func:`generation_degrees` -- the headline primitive: the internal (path-length)
  generation degrees ``l_i(n)`` of the minimal resolution of ``S_i``.  Plan 27 knows
  only the HOMOLOGICAL degrees of the Yoneda generators, and those do not see ``N``
  at all (``k[x]/x^3``, ``k[x]/x^4``, ``k[x]/x^5`` all have ``generators_by_degree ==
  {1: 1, 2: 1}``); the internal degrees are what distinguishes them.
* :func:`n_koszul_certificate` -- Berger's 2-N alternation on an N-homogeneous
  algebra (``N = 2`` defers to Plan 27 verbatim).
* :func:`k2_certificate` -- Cassidy-Shelton K2 through an explicit certified window.
* :func:`almost_koszul_certificate` -- the Brenner-Butler-King ``(p,q)`` classifier
  for the algebras a Koszul route refuses (the Dynkin preprojectives).
* :func:`multi_koszul_certificate` -- Herscovich's connected-graded notion, scoped.
* :func:`koszul_profile` / :func:`koszul_profile_block` -- the assembled record and
  the no-code block shared byte-identically by both runners.

Everything runs over any exact ``Domain`` (integer path lengths are field-free
combinatorics; the only linear algebra is inherited from Plan 05/27).

References: ``berger_nonquadratic`` (Berger, J. Algebra 239 (2001) 705-734),
``cassidy_shelton`` (Math. Z. 260 (2008) 93-114), ``brenner_butler_king``
(Algebr. Represent. Theory 5 (2002) 331-368), ``herscovich_multikoszul``
(arXiv:1305.1678), ``green_marcos_martinezvilla_zhang`` (JPAA 193 (2004) 141-162,
the semisimple-base foundation), ``chouhy_degenerations`` (context).
"""
from quiverlab.errors import DepthLimitError, QuiverlabError
from quiverlab.modules import koszul as _koszul
from quiverlab.modules.builders import projective, simple
from quiverlab.modules.resolution import minimal_resolution

__all__ = [
    "generation_degrees", "n_homogeneous_degree", "n_koszul_certificate",
    "k2_certificate", "almost_koszul_certificate", "multi_koszul_certificate",
    "koszul_profile", "koszul_profile_block",
]

_NOT_GRADED = ("not length-graded: internal degrees are undefined (the ideal I is "
               "inhomogeneous, so the minimal resolution is not graded)")


class _NonPure(Exception):
    """A generator column of ``d_n`` mixes internal degrees.

    DEFENSIVE ONLY.  On a length-graded ``A`` the minimal resolution is graded and
    every column of a differential is homogeneous, so this cannot fire; it guards a
    bug or a mis-gated non-graded input rather than validating the extraction.  The
    REAL correctness anchor for :func:`generation_degrees` is Berger's closed form
    ``delta(n)``, an EXTERNAL ground truth -- any independent re-derivation of the
    internal degrees would use the same single method (predecessor degree + the
    differential's path length), so "another extractor agrees" would prove nothing.
    """

    def __init__(self, degree, degrees):
        super().__init__(f"P_{degree} generator column mixes internal degrees "
                         f"{sorted(degrees)}")
        self.degree = degree
        self.degrees = sorted(degrees)


def _require_length_graded(A, what):
    _koszul._require_presentation(A, what)
    if not _koszul._is_length_graded(A):
        raise QuiverlabError(
            f"{what} requires a length-graded algebra -- {_NOT_GRADED}",
            hint="every defining relation must be homogeneous (all terms of one "
                 "path length); try k2_certificate, which needs only the "
                 "homological Yoneda generation degrees")


def _label_length(label, arrows, idempotents):
    """Path length of an algebra/projective basis label.

    The label grammar is ``e_<v>`` for a trivial path and a ``*``-joined arrow word
    otherwise.  The length is counted against the quiver's ARROW SET rather than by
    counting ``*`` separators, so an arrow whose name contains ``*`` or looks like an
    idempotent cannot silently corrupt a degree.
    """
    tokens = label.split("*")
    if all(t in arrows for t in tokens):
        return len(tokens)
    if label in idempotents:
        return 0
    raise QuiverlabError(
        f"unreadable path label {label!r}: not a vertex idempotent and not a "
        "'*'-joined word of quiver arrows",
        hint="the internal-degree extractor reads path lengths off the basis "
             "labels; a changed label grammar must be reflected here")


def _grammar(A):
    return set(A.quiver.arrows), {f"e_{v}" for v in A.quiver.vertices}


def _projective_basis_lengths(A, v, cache):
    """``(lengths, dim)`` for ``projective(A, v)``: the path length of each basis
    label, in the module's basis order.  Memoized per vertex in ``cache``."""
    if v not in cache:
        arrows, idems = _grammar(A)
        P = projective(A, v)
        cache[v] = ([_label_length(lab, arrows, idems) for lab in P._pv_basis_labels],
                    P.dim)
    return cache[v]


def _offsets(A, vertices, cache):
    """Block offsets ``(start, dim, vertex)`` of the summands of a resolution term."""
    out, start = [], 0
    for v in vertices:
        _, d = _projective_basis_lengths(A, v, cache)
        out.append((start, d, v))
        start += d
    return out


def _generation_degrees(A, i, length, max_term_dim, cache):
    """Worker for :func:`generation_degrees`; ``cache`` is shared across simples."""
    dom = A.domain
    terms, dmats = minimal_resolution(simple(A, i), length,
                                      max_term_dim=max_term_dim)
    out = []
    prev_offsets, prev_degrees = None, None
    for n, term in enumerate(terms):
        if not term.vertices:
            out.append([])                       # the resolution terminated here
            break
        offsets = _offsets(A, term.vertices, cache)
        if n == 0:
            # P_0 covers S_i in internal degree 0.
            degrees = [0] * len(offsets)
        else:
            dn = dmats[n]                        # rows index P_{n-1}, cols index P_n
            degrees = []
            for (start, _dim, _v) in offsets:
                column = start                   # the summand's canonical generator
                seen = set()
                for r in range(len(dn)):
                    if dom.is_zero(dn[r][column]):
                        continue
                    for k, (pstart, pdim, pv) in enumerate(prev_offsets):
                        if pstart <= r < pstart + pdim:
                            lengths, _ = _projective_basis_lengths(A, pv, cache)
                            seen.add(prev_degrees[k] + lengths[r - pstart])
                            break
                if len(seen) != 1:
                    raise _NonPure(n, seen)
                degrees.append(seen.pop())
        out.append(sorted(set(degrees)))
        prev_offsets, prev_degrees = offsets, degrees
    return out


def generation_degrees(A, i, length, *, max_term_dim=200_000):
    """Internal (path-length) generation degrees of the minimal resolution of ``S_i``.

    Returns one entry per homological degree ``n = 0, 1, ...``: the SORTED SET of
    internal degrees in which the summand generators of ``P_n`` sit.  A terminated
    resolution contributes a final empty entry (``P_n = 0``).

    For a length-graded ``A`` the resolution is graded, so every differential is
    homogeneous and the internal degree of a ``P_n`` summand generator is well
    defined: ``P_0`` covers ``S_i`` in degree ``0``, and for ``n >= 1`` a summand's
    generator column has all its nonzero entries in one internal degree, namely the
    predecessor generator's degree plus the path length of the entry's basis vector.

    Raises ``QuiverlabError`` when ``A`` is not length-graded (internal degrees are
    then undefined) and ``DepthLimitError`` when a syzygy overshoots ``max_term_dim``
    (the caller reports the truncated window honestly).
    """
    _require_length_graded(A, "generation_degrees")
    return _generation_degrees(A, i, length, max_term_dim, {})


def n_homogeneous_degree(A):
    """The single relation length ``N >= 2`` when ``A`` is N-homogeneous, else ``None``.

    N-homogeneous means every defining relation is homogeneous of the SAME path
    length ``N``.  A HEREDITARY ``kQ`` (no relations, including the semisimple case)
    returns ``2``: it is the quadratic algebra with an empty relation space, matching
    ``koszul.is_quadratic``'s vacuous-True convention, so callers get the ``N = 2``
    defer-to-Plan-27 branch with no separate sentinel.
    """
    _koszul._require_presentation(A, "n_homogeneous_degree")
    relations = A.relations or []
    if not relations:
        return 2
    degrees = set()
    for rel in relations:
        if rel.min_length != rel.max_length:
            return None                          # an inhomogeneous relation
        degrees.add(rel.min_length)
    return degrees.pop() if len(degrees) == 1 else None


def berger_degree(n, N):
    """Berger's pure-resolution internal degree ``delta(n)`` for an N-Koszul algebra:
    ``(N/2)*n`` for ``n`` even, ``(N/2)*(n-1) + 1`` for ``n`` odd -- written in exact
    integer form so no division is ever inexact."""
    return (n // 2) * N if n % 2 == 0 else (n // 2) * N + 1


def _yoneda(A, top, yoneda):
    """The Plan-27 ``YonedaPresentation``, computed once and threaded by the profile."""
    return yoneda if yoneda is not None else A.ext_algebra(top)


def _all_generation_degrees(A, length, max_term_dim=200_000):
    """``{vertex: generation_degrees(...)}`` for every simple, sharing one cache."""
    _require_length_graded(A, "generation_degrees")
    cache = {}
    return {v: _generation_degrees(A, v, length, max_term_dim, cache)
            for v in A.quiver.vertices}


# ---------------------------------------------------------------------------
# K2 (Cassidy-Shelton)
# ---------------------------------------------------------------------------
def k2_certificate(A, top=8, *, yoneda=None):
    """Cassidy-Shelton K2: is ``E(A) = Ext(k,k)`` generated as an algebra in
    cohomological degrees 1 and 2?

    Decided from the shipped Plan-27 ``generators_by_degree`` through the EXPLICIT
    certified window ``W = certified_through_degree``.  Three-valued:

    * ``False`` -- a genuine Yoneda generator sits in degree ``>= 3`` INSIDE ``W``.
      Decisive regardless of completeness (the obstruction is exhibited).
    * ``True`` -- the generators sit in degrees ``{1, 2}`` AND the window is complete
      (``gl.dim A`` finite and exact, so ``E(A)`` is finite-dimensional and fully
      computed).
    * ``None`` -- the generators sit in degrees ``{1, 2}`` through ``W`` but the
      window is NOT complete: honest inconclusive, a generator could still appear
      beyond ``W``.  Never reported as an unconditional ``True``.

    K2 needs only HOMOLOGICAL data, so unlike the N-Koszul / almost-Koszul
    recognizers it does not require ``A`` to be length-graded.
    """
    Y = _yoneda(A, top, yoneda)
    degrees = sorted(Y.generators_by_degree)
    window = Y.certified_through_degree
    complete = Y.is_finite_dimensional is True
    beyond = [d for d in degrees if d >= 3]
    if beyond:
        verdict = False
        reason = ("a Yoneda generator sits in cohomological degree "
                  + str(beyond[0]) + ", inside the certified window "
                  + str(window) + ": E(A) is NOT generated in degrees 1 and 2")
    elif complete:
        verdict = True
        reason = ("gl.dim A is finite and exact, so E(A) is finite-dimensional and "
                  "completely computed through degree " + str(window)
                  + "; its generators sit in degrees "
                  + (str(degrees) if degrees else "{} (semisimple: vacuously K2)"))
    else:
        verdict = None
        reason = ("K2 holds through degree " + str(window)
                  + ", but gl.dim A is not exact-finite, so a Yoneda generator "
                    "could still appear beyond the window: inconclusive, not a "
                    "claim of K2")
    return {"verdict": verdict, "window": window, "complete": complete,
            "generator_degrees": degrees, "reason": reason}


# ---------------------------------------------------------------------------
# N-Koszul (Berger)
# ---------------------------------------------------------------------------
def n_koszul_certificate(A, top=8, *, yoneda=None):
    """Berger's N-Koszul recognizer: is the minimal resolution of every simple PURE,
    with ``P_n`` generated in the single internal degree ``delta(n)``?

    ``N = 2`` DEFERS TO PLAN 27 verbatim (``ext_algebra(top).koszul``) -- Plan 77
    adds no second opinion in the quadratic case, and the hereditary ``kQ`` lands in
    this branch too (``n_homogeneous_degree`` returns 2 vacuously).

    For ``N >= 3`` the verdict is cross-checked against Berger's Ext-algebra
    characterization (N-Koszul iff the Yoneda algebra is generated in degrees 0, 1, 2
    -- i.e. iff K2), which is a genuinely independent reading of the same algebra:
    the purity side reads INTERNAL degrees off the minimal resolution, the K2 side
    reads HOMOLOGICAL generator degrees off the Yoneda engine.  ``k2_agrees`` records
    whether they agree; a disagreement is surfaced, never smoothed over.
    """
    Y = _yoneda(A, top, yoneda)
    window = Y.certified_through_degree
    complete = Y.is_finite_dimensional is True
    out = {"n_homogeneous": None, "verdict": None, "window": window,
           "complete": complete, "berger_expected": [], "internal_degrees": {},
           "internal_degree_sets": {}, "pure": None, "impurity": None,
           "k2_agrees": None, "reason": ""}

    N = n_homogeneous_degree(A)
    out["n_homogeneous"] = N
    if N is None:
        out["reason"] = ("relations are not homogeneous of a single degree, so A is "
                         "not N-homogeneous and Berger's N-Koszul property does not "
                         "apply; see the K2 and multi-Koszul fields")
        return out
    out["berger_expected"] = [berger_degree(n, N) for n in range(window + 1)]

    if N == 2:
        out["verdict"] = Y.koszul
        out["reason"] = ("N=2: quadratic Koszulity is Plan 27's "
                         "g_quadratic_certificate/Froberg verdict, reported verbatim "
                         "(a hereditary kQ is quadratic-trivially); Plan 77 adds no "
                         "second opinion in the quadratic case")
        return out

    if not _koszul._is_length_graded(A):
        # Unreachable for N >= 3 (a single relation degree IS homogeneous), kept as
        # an explicit honest branch rather than an implicit assumption.
        out["reason"] = _NOT_GRADED
        return out

    try:
        per_simple = _all_generation_degrees(A, max(window, 1))
    except DepthLimitError as exc:
        out["reason"] = ("the minimal resolution could not be built through the "
                         "window: " + str(exc))
        return out
    except _NonPure as exc:
        out["pure"] = False
        out["verdict"] = False
        out["impurity"] = {"degree": exc.degree, "internal_degrees": exc.degrees}
        out["reason"] = ("a generator column of d_" + str(exc.degree) + " mixes "
                         "internal degrees " + str(exc.degrees)
                         + ": the resolution is not graded-pure")
        return out

    out["internal_degree_sets"] = {v: [list(s) for s in gd]
                                   for v, gd in per_simple.items()}
    flat, pure, impurity, mismatch = {}, True, None, None
    for v, gd in per_simple.items():
        row = []
        for n, s in enumerate(gd):
            if not s:
                break                            # the resolution terminated
            if len(s) != 1:
                pure = False
                impurity = {"simple": v, "degree": n, "internal_degrees": s}
                break
            row.append(s[0])
            if s[0] != berger_degree(n, N) and mismatch is None:
                mismatch = {"simple": v, "degree": n, "internal_degree": s[0],
                            "berger": berger_degree(n, N)}
        flat[v] = row
        if impurity is not None:
            break
    out["internal_degrees"] = flat
    out["pure"] = pure and mismatch is None
    out["impurity"] = impurity

    k2 = k2_certificate(A, top, yoneda=Y)
    k2_holds = k2["verdict"] is not False        # no degree->=3 generator inside W
    out["k2_agrees"] = (k2_holds == out["pure"])

    if impurity is not None:
        out["verdict"] = False
        out["reason"] = ("P_" + str(impurity["degree"]) + " of S_"
                         + str(impurity["simple"]) + " is generated in multiple "
                         "internal degrees " + str(impurity["internal_degrees"])
                         + " rather than the single Berger degree delta("
                         + str(impurity["degree"]) + ") = "
                         + str(berger_degree(impurity["degree"], N))
                         + ": the resolution is not pure")
    elif mismatch is not None:
        out["verdict"] = False
        out["reason"] = ("P_" + str(mismatch["degree"]) + " of S_"
                         + str(mismatch["simple"]) + " is generated in internal "
                         "degree " + str(mismatch["internal_degree"])
                         + ", not Berger's delta(" + str(mismatch["degree"]) + ") = "
                         + str(mismatch["berger"]) + " for N=" + str(N)
                         + ": the 2-N alternation breaks")
    elif not out["k2_agrees"]:
        # Purity holds through the window yet the Yoneda algebra has a generator in
        # degree >= 3 -- Berger's N>=3 equivalence says these cannot both hold. Report
        # the disagreement instead of picking a side.
        out["verdict"] = False
        out["reason"] = ("the two N-Koszul certificates DISAGREE through the window: "
                         "the resolution is pure with Berger's degrees, but "
                         + k2["reason"] + " -- Berger's N>=3 characterization "
                         "(N-Koszul iff E(A) is generated in degrees 0,1,2) makes "
                         "this a genuine inconsistency, reported rather than resolved")
    elif complete:
        out["verdict"] = True
        out["reason"] = ("every simple has a pure minimal resolution with Berger's "
                         "internal degrees delta(n) for N=" + str(N)
                         + ", and gl.dim A is finite and exact, so the pattern is "
                           "complete: A is N-Koszul")
    elif k2["verdict"] is True:
        out["verdict"] = True
        out["reason"] = ("the resolution is pure with Berger's degrees and E(A) is "
                         "completely generated in degrees 0,1,2 (Berger's N>=3 "
                         "characterization): A is N-Koszul")
    else:
        out["verdict"] = None
        out["reason"] = ("the Berger pattern delta(n) for N=" + str(N)
                         + " holds through degree " + str(window)
                         + ", but gl.dim A is not exact-finite and K2 is only "
                           "window-certified, so N-Koszulity is certified THROUGH "
                           "THE WINDOW only, never claimed unconditionally")
    return out


# ---------------------------------------------------------------------------
# (p,q)-almost-Koszul (Brenner-Butler-King)
# ---------------------------------------------------------------------------
def _top_degree(A):
    """``p`` = the top graded degree of ``A``: the longest path in its own basis.

    ``A.basis_labels`` IS the irreducible-path basis, so for a length-graded ``A``
    the concentration ``A = A_0 (+) ... (+) A_p`` is read off directly and no Groebner
    completion is involved.  (``koszul._algebra_graded_matrices`` would re-complete the
    reduction system under its default ``degree_bound=8``, which is too small for
    ``Pi(A5)``/``Pi(D4)`` -- it raises ``AdmissibilityError`` there.  The basis-label
    route is total on every finite-dimensional ``A``.)"""
    _koszul._require_presentation(A, "_top_degree")
    arrows, idems = _grammar(A)
    return max([0] + [_label_length(lab, arrows, idems) for lab in A.basis_labels])


def _first_jump(degrees):
    """``(n_star, e)`` -- the first homological degree whose internal degree exceeds it,
    and that internal degree.  ``None`` when the strand stays linear (``l(n) == n``)
    for every computed term: a Koszul resolution has no break.  Raises ``ValueError``
    when the breaking term is generated in SEVERAL internal degrees (``e`` ambiguous).
    """
    for n, s in enumerate(degrees):
        if not s:
            return None                          # the resolution terminated, linear
        if s == [n]:
            continue
        if len(s) != 1:
            raise ValueError("P_" + str(n) + " is generated in several internal "
                             "degrees " + str(s) + ", so the break degree e is "
                             "ambiguous")
        return n, s[0]
    return None


def almost_koszul_certificate(A, top=8, *, yoneda=None):
    """Brenner-Butler-King ``(p,q)``-almost-Koszul, as a SIGNATURE recognizer.

    BBK: ``A`` is concentrated in degrees ``0..p`` and a linear complex of projectives
    resolves the simple up to an error given by the degree ``p+q`` part.  Pointwise on
    the minimal resolution that reads: ``l_i(n) = n`` (a linear strand) for ``n <= q``,
    then the first break at ``n = q+1`` to internal degree ``e = p+q``.  The engine
    recovers ``p`` as the top degree of ``A`` and the break ``(n*, e)`` from the
    internal degrees, giving TWO readings of BBK's ``q`` -- ``q = e - p`` and
    ``q = n* - 1`` -- which must AGREE (they do on every Dynkin preprojective and on
    every ``k[x]/x^N``); a disagreement is reported, never averaged away.

    ``verdict=True`` requires: a linear prefix, a first break with ``(n*, e)`` agreeing
    across the simples that break, the two ``q`` readings agreeing, and ``q >= 2``.
    ``verdict=None`` when there is NO break (the resolution is linear -- the algebra is
    genuinely Koszul, e.g. ``Pi(A2)`` or ``k[x]/x^2``, and a Koszul algebra needs no
    "error"), when ``q == 1`` (BBK's degenerate Koszul-type boundary -- the whole
    ``k[x]/x^N`` family, which is N-Koszul and is caught by
    :func:`n_koszul_certificate`), or when the data disagree across simples.

    **Scoped OUT, stated:** this does NOT verify BBK's full definition -- neither the
    existence of the finite linear complex of projectives with a single error term nor
    the ``2(h-1)`` periodicity.  What is reported is the ``(p,q)`` label plus the
    single-break structure, which is what BBK's definition checks pointwise, plus the
    WINDOW-OBSERVED spacing of the internal-degree jumps -- which is NOT BBK's bimodule
    period ``2(h-1)`` and is never claimed to be (for ``Pi(A3)``, ``h = 4`` and BBK's
    period is 6, while the observed jump spacing is 3).
    """
    Y = _yoneda(A, top, yoneda)
    window = Y.certified_through_degree
    obstruction = Y.koszul_obstruction
    out = {"verdict": None, "p": None, "q": None, "break_hom_degree": None,
           "break_internal_degree": None, "linear_steps": None,
           "jump_degrees": [], "jump_spacing": None, "window": window,
           "complete": Y.is_finite_dimensional is True,
           "seam_obstruction_degree": obstruction[0] if obstruction else None,
           "reason": ""}

    if not _koszul._is_length_graded(A):
        out["reason"] = _NOT_GRADED
        return out
    p = _top_degree(A)
    out["p"] = p
    try:
        per_simple = _all_generation_degrees(A, max(window, 1))
    except DepthLimitError as exc:
        out["reason"] = ("the minimal resolution could not be built through the "
                         "window: " + str(exc))
        return out
    except _NonPure as exc:
        out["reason"] = ("a generator column of d_" + str(exc.degree) + " mixes "
                         "internal degrees " + str(exc.degrees)
                         + ": the break degree is not well defined")
        return out

    breaks = {}
    for v, gd in per_simple.items():
        try:
            jump = _first_jump(gd)
        except ValueError as exc:
            out["reason"] = "S_" + str(v) + ": " + str(exc)
            return out
        if jump is not None:
            breaks[v] = jump

    if not breaks:
        out["reason"] = ("the minimal resolution of every simple stays linear "
                         "(l(n) = n) through degree " + str(window)
                         + ": there is no almost-Koszul error term -- the algebra is "
                           "Koszul in this range, and a Koszul algebra needs none "
                           "(see the quadratic and K2 verdicts)")
        return out

    n_stars = {n for (n, _e) in breaks.values()}
    es = {e for (_n, e) in breaks.values()}
    if len(n_stars) != 1 or len(es) != 1:
        out["reason"] = ("the simples break in different places -- (n*, e) = "
                         + str(sorted(breaks.items()))
                         + " -- so there is no single (p,q) label")
        return out
    witness = sorted(breaks)[0]
    n_star, e = breaks[witness]
    q = e - p
    out.update({"break_hom_degree": n_star, "break_internal_degree": e,
                "q": q, "linear_steps": n_star - 1})

    # The window-observed jump spacing: the homological degrees at which the excess
    # l(n) - n increases. Recorded; never claimed to be BBK's period 2(h-1).
    jumps, previous = [], 0
    for n, s in enumerate(per_simple[witness]):
        if not s or len(s) != 1:
            break
        excess = s[0] - n
        if excess > previous:
            jumps.append(n)
        previous = excess
    out["jump_degrees"] = jumps
    spacings = {b - a for a, b in zip(jumps, jumps[1:])}
    out["jump_spacing"] = spacings.pop() if len(spacings) == 1 else None

    if q != n_star - 1:
        out["reason"] = ("the two readings of BBK's q disagree: e - p = " + str(q)
                         + " but n* - 1 = " + str(n_star - 1)
                         + "; BBK's definition puts the error of internal degree p+q "
                           "at homological step q+1, so this is not the "
                           "almost-Koszul signature")
        return out
    if q < 2:
        out["reason"] = ("q = e - p = " + str(q) + " < 2: the error sits at internal "
                         "degree p+1, BBK's degenerate Koszul-type boundary, not a "
                         "genuine almost-Koszul break (an N-homogeneous algebra here "
                         "is N-Koszul -- see the n_koszul field)")
        return out
    out["verdict"] = True
    out["reason"] = ("A is concentrated in degrees 0.." + str(p)
                     + " and the minimal resolution of every simple is linear for "
                     + str(n_star - 1) + " steps and then breaks once, at homological "
                       "degree " + str(n_star) + " to internal degree " + str(e)
                     + " = p + q: the (p, q) = (" + str(p) + ", " + str(q)
                     + ")-almost-Koszul signature (Brenner-Butler-King). The full BBK "
                       "definition -- the finite linear complex with a single error "
                       "term and the 2(h-1) periodicity -- is NOT verified here")
    return out
