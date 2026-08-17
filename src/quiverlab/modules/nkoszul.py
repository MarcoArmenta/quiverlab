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
