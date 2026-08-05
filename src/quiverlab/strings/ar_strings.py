"""String-module AR translate by hooks/cohooks (Butler-Ringel 1987).

The AR translate of a string module is again a string module, described
combinatorially by Butler-Ringel via HOOKS and COHOOKS at the two ends of the
string:

  - *adding a cohook* at an end = append a DIRECT arrow, then as many INVERSE
    arrows as the string conditions allow (a maximal inverse arm);
  - *adding a hook* at an end = append an INVERSE arrow, then a maximal DIRECT arm;
  - deleting a hook / cohook are the inverse operations.

``tau M(w)`` performs, at each end, "delete a hook if the end terminates in one,
else add a cohook"; ``tau^- M(w)`` is the dual ("delete a cohook else add a hook").
A projective string has ``tau = 0``; an injective string has ``tau^- = 0``.

CONVENTION / SCOPE (arbitrated, not assumed -- the house cup-sign precedent): the
hook/cohook *side* and the delete-vs-add choice at each end are ARBITRATED against
the trusted Plan-23/41 engine translate (``Module.tau``) on EVERY call. The
closed-form arm rule (``_combinatorial_tau``) reproduces the engine translate for
the string chains this subsystem targets (linear A_n, gentle A_n/rad^2 and their
relatives -- verified degreewise); at a few AR-quiver BOUNDARY strings (the
translate of a non-injective projective sitting at a source/sink, a structural
"jump" rather than a local arm move) the closed-form arm operation does not build
the correct string, so the result is COMPLETED from the trusted engine translate by
identifying its string. The engine translate is the arbiter throughout -- we never
return a walk whose module is not ``is_isomorphic`` to it, and never a guessed
convention. Butler-Ringel 1987, Communications in Algebra 15."""
from __future__ import annotations

from quiverlab.errors import QuiverlabError
from quiverlab.strings.modules import string_module
from quiverlab.strings.walks import (_canonical, _is_trivial, _reverse_inverse,
                                     enumerate_strings, is_valid_walk,
                                     letter_source, letter_target)


# ---------------------------------------------------------------------------
# Hook / cohook primitives (right end; the left end is the mirror)
# ---------------------------------------------------------------------------


def _right_vertex(A, w):
    return w[0][1] if _is_trivial(w) else letter_target(A.quiver, w[-1])


def _max_arm(A, w, sign):
    """Append letters of the given ``sign`` (+1 direct, -1 inverse) greedily while the
    walk stays valid (a maximal arm)."""
    Q = A.quiver
    while True:
        e = _right_vertex(A, w)
        nxt = None
        for a in Q.arrows:
            ell = (a, sign)
            if letter_source(Q, ell) == e:
                cand = (ell,) if _is_trivial(w) else w + (ell,)
                if is_valid_walk(A, cand):
                    nxt = cand
                    break
        if nxt is None:
            return w
        w = nxt


def _add_cohook_right(A, w):
    """Append a direct arrow then a maximal inverse arm."""
    Q = A.quiver
    e = _right_vertex(A, w)
    for a in Q.arrows:
        ell = (a, +1)
        if letter_source(Q, ell) == e:
            w1 = (ell,) if _is_trivial(w) else w + (ell,)
            if is_valid_walk(A, w1):
                return _max_arm(A, w1, -1)
    return None


def _add_hook_right(A, w):
    """Append an inverse arrow then a maximal direct arm."""
    Q = A.quiver
    e = _right_vertex(A, w)
    for a in Q.arrows:
        ell = (a, -1)
        if letter_source(Q, ell) == e:
            w1 = (ell,) if _is_trivial(w) else w + (ell,)
            if is_valid_walk(A, w1):
                return _max_arm(A, w1, +1)
    return None


def _collapse(A, w, removed):
    return w if w else ((None, letter_source(A.quiver, removed)),)


def _delete_hook_right(A, w):
    """Remove a trailing maximal direct arm then one inverse letter (inverse of
    adding a hook). ``None`` if the end does not terminate in a hook."""
    if _is_trivial(w):
        return None
    w = tuple(w)
    while w and w[-1][1] > 0:
        w = w[:-1]
    if not w or w[-1][1] != -1:
        return None
    removed = w[-1]
    return _collapse(A, w[:-1], removed)


def _delete_cohook_right(A, w):
    """Remove a trailing maximal inverse arm then one direct letter (inverse of
    adding a cohook). ``None`` if the end does not terminate in a cohook."""
    if _is_trivial(w):
        return None
    w = tuple(w)
    while w and w[-1][1] < 0:
        w = w[:-1]
    if not w or w[-1][1] != 1:
        return None
    removed = w[-1]
    return _collapse(A, w[:-1], removed)


def _mirror(op, A, w):
    """Apply a right-end operation to the LEFT end via reverse-inversion."""
    rw = w if _is_trivial(w) else _reverse_inverse(w)
    r = op(A, rw)
    if r is None:
        return None
    return r if _is_trivial(r) else _reverse_inverse(r)


def _ends_with_hook(w):
    if _is_trivial(w):
        return False
    ww = tuple(w)
    while ww and ww[-1][1] > 0:
        ww = ww[:-1]
    return bool(ww) and ww[-1][1] < 0


def _ends_with_cohook(w):
    if _is_trivial(w):
        return False
    ww = tuple(w)
    while ww and ww[-1][1] < 0:
        ww = ww[:-1]
    return bool(ww) and ww[-1][1] > 0


def _tau_right(A, w):
    if not _is_trivial(w) and _ends_with_hook(w):
        return _delete_hook_right(A, w)
    return _add_cohook_right(A, w)


def _taum_right(A, w):
    if not _is_trivial(w) and _ends_with_cohook(w):
        return _delete_cohook_right(A, w)
    return _add_hook_right(A, w)


def _combinatorial_tau(A, w):
    r = _tau_right(A, w)
    if r is None:
        return None
    return _mirror(_tau_right, A, r)


def _combinatorial_tau_minus(A, w):
    r = _taum_right(A, w)
    if r is None:
        return None
    return _mirror(_taum_right, A, r)


# ---------------------------------------------------------------------------
# Module <-> walk identification (engine-arbitrated)
# ---------------------------------------------------------------------------


def _matches(A, walk, T):
    if walk is None:
        return False
    try:
        M = string_module(A, walk)
    except QuiverlabError:
        return False
    from quiverlab.modules.hom import is_isomorphic
    return (M.dim == T.dim and M.dimension_vector() == T.dimension_vector()
            and is_isomorphic(M, T))


def _walk_of_module(A, T, prefer=None):
    """Return a canonical walk ``w`` with ``string_module(w)`` isomorphic to ``T``.

    Prefers the combinatorial candidate ``prefer`` when it verifies against ``T``;
    otherwise identifies the string by searching the census (the AR translate of a
    string is a string, so it is in the census of a rep-finite string algebra).
    Loud if no string matches (a rep-infinite band region is out of scope)."""
    if _matches(A, prefer, T):
        return _canonical(prefer)
    from quiverlab.modules.hom import is_isomorphic
    ml = max(8, 2 * T.dim)
    cen = enumerate_strings(A, max_length=ml, budget=1 << 16)
    dvT = T.dimension_vector()
    for w in cen.walks:
        M = string_module(A, w)
        if M.dim == T.dim and M.dimension_vector() == dvT and is_isomorphic(M, T):
            return w
    raise QuiverlabError(
        "strings: could not express the AR translate as a string module",
        hint="the translate is outside the enumerated string census (a rep-infinite "
             "band region is out of scope for the combinatorial string translate)")


# ---------------------------------------------------------------------------
# Public surface
# ---------------------------------------------------------------------------


def string_tau(A, walk):
    """The walk ``w'`` with ``string_module(w') ~ tau(string_module(walk))``, by the
    Butler-Ringel hook/cohook combinatorics -- arbitrated against the trusted engine
    translate. ``None`` when ``M(walk)`` is projective (``tau = 0``).

    HONESTY NOTE (devil's-advocate round, 2026-08-05): this function COMPUTES
    the engine translate ``Module.tau()`` on every call and then identifies a
    walk presentation for it (combinatorial rule first, census fallback), with
    ``is_isomorphic`` verification before anything is returned. It is neither
    faster than nor independent of the engine tau -- its value is the WALK
    (the combinatorial presentation the engine does not provide). It can raise
    where the engine succeeds (no walk found within budget); it can never
    return a wrong walk."""
    T = string_module(A, walk).tau()
    if T.dim == 0:
        return None
    return _walk_of_module(A, T, prefer=_combinatorial_tau(A, walk))


def string_tau_minus(A, walk):
    """The walk ``w'`` with ``string_module(w') ~ tau^-(string_module(walk))``.
    ``None`` when ``M(walk)`` is injective (``tau^- = 0``)."""
    T = string_module(A, walk).tau_minus()
    if T.dim == 0:
        return None
    return _walk_of_module(A, T, prefer=_combinatorial_tau_minus(A, walk))
