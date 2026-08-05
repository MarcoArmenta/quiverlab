"""Reduced walks / strings / bands over a string algebra (Plan 46 / C5).

A string algebra is ``kQ/I`` with ``I`` monomial admissible, at most 2 arrows in
and 2 out at each vertex, and the special-biserial branching condition (for each
arrow at most one nonzero and one zero length-2 continuation on each side). Its
indecomposable modules are the STRING modules (from reduced walks avoiding ``I``)
and the BAND modules (from cyclic reduced walks); Butler-Ringel 1987. Everything
here is exact combinatorics on the reduction system -- float-free, field-agnostic.

A **letter** is ``(arrow_name: str, direction: +1|-1)``: ``+1`` a direct arrow,
``-1`` a formal inverse. A **walk** is a tuple of letters read left-to-right in the
WALK sense (letter ``i`` ends where letter ``i+1`` starts; inverses reverse
endpoints). A **trivial walk** at a vertex ``v`` is written ``((None, v),)`` -- the
one-vertex string whose module is the simple ``S_v`` (a bare ``()`` cannot name its
vertex, per the Task-1 bookkeeping note).

Canonicalisation: a walk ``w`` and its formal inverse ``w^{-1}`` give isomorphic
modules, so a census keeps the lexicographically smaller of the two (``_canonical``).
A trivial walk is its own inverse.
"""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.invariants.recognizers import _len2_in_ideal, is_string
from quiverlab.resolutions_cs.build import reduction_system_of

# ---------------------------------------------------------------------------
# Letter algebra
# ---------------------------------------------------------------------------


def invert(ell):
    """The formal inverse of a letter. A trivial letter ``(None, v)`` is its own
    inverse (reversing a one-vertex walk yields itself)."""
    name, d = ell
    return ell if name is None else (name, -d)


def letter_source(Q, ell):
    name, d = ell
    return Q.source(name) if d > 0 else Q.target(name)


def letter_target(Q, ell):
    name, d = ell
    return Q.target(name) if d > 0 else Q.source(name)


def _is_trivial(walk):
    return len(walk) == 1 and walk[0][0] is None


# ---------------------------------------------------------------------------
# String signs sigma / epsilon (Butler-Ringel S1/S2/S3)
# ---------------------------------------------------------------------------

_SIGNS_CACHE = {}


def _require_string(A):
    if A.quiver is None or A.relations is None:
        raise QuiverlabError("strings: need a quiver-presented algebra",
                             hint="structure-constant algebras carry no walks")
    if not is_string(A):
        raise QuiverlabError(
            "strings: this algebra is not a string algebra",
            hint="strings/bands are defined for special-biserial monomial kQ/I "
                 "with <=2 arrows in/out per vertex (is_string(A) is False)")


def string_signs(A):
    """``(sigma, epsilon)``: ``Q_1 -> {+1,-1}``.

    (S1) arrows with equal SOURCE get distinct ``sigma``; (S2) arrows with equal
    TARGET get distinct ``epsilon``; (S3) for ``beta*gamma`` NOT in ``I`` with
    ``target(beta)=source(gamma)``, ``sigma(gamma) = -epsilon(beta)``. Built
    greedily from the (<=2) branches at each vertex and reconciled for (S3); a
    string algebra always admits such an assignment. Loud if ``A`` is not a string
    algebra, or the constraints are inconsistent. Memoised by ``id(A)``."""
    _require_string(A)
    key = id(A)
    hit = _SIGNS_CACHE.get(key)
    if hit is not None and hit[0] is A:
        return hit[1]
    Q = A.quiver
    rs = reduction_system_of(A)
    arrows = list(Q.arrows)
    sigma, epsilon = {}, {}
    for v in Q.vertices:
        outs = [a for a in arrows if Q.source(a) == v]
        for i, a in enumerate(outs):
            sigma[a] = 1 if i == 0 else -1
        ins = [a for a in arrows if Q.target(a) == v]
        for i, a in enumerate(ins):
            epsilon[a] = 1 if i == 0 else -1
    # (S3): reconcile via the nonzero composable pairs. Propagate; if a forced
    # value contradicts an assigned one even after a single class-flip, raise.
    for b in arrows:
        for g in arrows:
            if Q.target(b) == Q.source(g) and not _len2_in_ideal(rs, b, g):
                want = -epsilon[b]
                if sigma[g] != want:
                    _flip_sigma_class(Q, sigma, Q.source(g))
                    if sigma[g] != want:
                        raise QuiverlabError(
                            "strings: sigma/epsilon constraints are inconsistent",
                            hint=f"arrows {b!r}, {g!r} force a contradictory sign; "
                                 "the algebra may violate the string branch condition")
    result = (sigma, epsilon)
    _SIGNS_CACHE[key] = (A, result)
    return result


def _flip_sigma_class(Q, sigma, v):
    for a in Q.arrows:
        if Q.source(a) == v:
            sigma[a] = -sigma[a]


# ---------------------------------------------------------------------------
# Walk validity
# ---------------------------------------------------------------------------


def is_valid_walk(A, walk, rs=None):
    """Composable, reduced, relation-avoiding, and sign-consistent."""
    Q = A.quiver
    if not walk or _is_trivial(walk):
        return True                                   # trivial walk e_v is valid
    if rs is None:
        rs = reduction_system_of(A)
    for i in range(len(walk) - 1):
        if letter_target(Q, walk[i]) != letter_source(Q, walk[i + 1]):
            return False                              # not composable
        if walk[i + 1] == invert(walk[i]):
            return False                              # not reduced (backtrack)
        if not _pair_ok(rs, walk[i], walk[i + 1]):
            return False                              # hits a relation
    return _signs_consistent(A, walk)


def _pair_ok(rs, ell1, ell2):
    (n1, d1), (n2, d2) = ell1, ell2
    if d1 > 0 and d2 > 0:
        return not _len2_in_ideal(rs, n1, n2)         # direct a*b must be nonzero
    if d1 < 0 and d2 < 0:
        return not _len2_in_ideal(rs, n2, n1)         # inverse pair: reverse b*a
    return True                                       # a direct<->inverse turn is
                                                      # governed by sign-consistency


def _signs_consistent(A, walk):
    """A direct->inverse (peak) or inverse->direct (valley) turn must use two
    DISTINCT arrows on the shared side -- exactly the Butler-Ringel condition that
    the string module is well-defined. For special-biserial (<=2 in/out) this is
    equivalent to the two arrows being distinct, which the sigma/epsilon distinctness
    encodes."""
    sigma, epsilon = string_signs(A)
    for i in range(len(walk) - 1):
        (n1, d1), (n2, d2) = walk[i], walk[i + 1]
        if d1 > 0 and d2 < 0:                          # peak: ... a  b^-1 ...
            if n1 == n2 or epsilon[n1] == epsilon[n2]:
                return False
        elif d1 < 0 and d2 > 0:                        # valley: ... a^-1 b ...
            if n1 == n2 or sigma[n1] == sigma[n2]:
                return False
    return True


# ---------------------------------------------------------------------------
# Canonical forms
# ---------------------------------------------------------------------------


def _reverse_inverse(walk):
    return tuple(invert(e) for e in reversed(walk))


def _canonical(walk):
    """The lexicographically smaller of ``walk`` and its reverse-inverse. A trivial
    walk is its own reverse-inverse, so it is returned unchanged."""
    if _is_trivial(walk):
        return tuple(walk)
    return min(tuple(walk), _reverse_inverse(walk))


def _sort_key(walk):
    # None-safe key so trivial ((None, v),) and honest letters never compare
    # None against a str name.
    return tuple((nm is None, "" if nm is None else nm, d) for nm, d in walk)


def _rotations(walk):
    L = len(walk)
    return [walk[k:] + walk[:k] for k in range(L)]


def _is_proper_power(walk):
    L = len(walk)
    for d in range(1, L):
        if L % d == 0 and walk[:d] * (L // d) == tuple(walk):
            return True
    return False


def _canonical_cyclic(walk):
    """Rotation- and inversion-minimal representative of a cyclic walk."""
    candidates = []
    for rot in _rotations(tuple(walk)):
        candidates.append(rot)
        candidates.append(_reverse_inverse(rot))
    return min(candidates, key=lambda w: (_sort_key(w), w))


# ---------------------------------------------------------------------------
# Census
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StringCensus:
    walks: tuple
    status: str          # "complete" | "budget"
    max_length: int
    count: int
    has_bands: bool


def enumerate_strings(A, max_length=8, budget=4096):
    """Enumerate the string modules of a string algebra as canonical reduced walks.

    HONEST CONTRACT: ``status == "complete"`` ONLY when ``find_bands(A, max_length)``
    is empty (rep-finite) AND the depth-first search closed within ``budget`` without
    any walk hitting ``max_length`` -- then the list is ALL indecomposable string
    modules (every non-band indecomposable). Bands present, or a budget/length cut,
    => ``status == "budget"`` (a length-capped sample, never a "complete" claim).

    Every valid walk is reached as a right-extension prefix from its first letter
    (all prefixes of a valid walk are valid), so no left-growth is needed; the output
    is canonicalised up to inversion.
    """
    _require_string(A)
    rs = reduction_system_of(A)
    Q = A.quiver
    alphabet = [(a, +1) for a in Q.arrows] + [(a, -1) for a in Q.arrows]
    out = set()
    state = {"truncated": False}

    for v in Q.vertices:                               # trivial walks = the simples
        out.add(((None, v),))

    def dfs(walk):
        out.add(_canonical(walk))
        if len(out) > budget:
            state["truncated"] = True
            return
        children = [walk + (ell,) for ell in alphabet]
        children = [c for c in children if is_valid_walk(A, c, rs)]
        if len(walk) >= max_length:
            if children:
                state["truncated"] = True
            return
        for c in children:
            dfs(c)

    for ell in alphabet:
        dfs((ell,))

    bands = find_bands(A, max_length)
    walks = tuple(sorted(out, key=_sort_key))
    complete = (not bands) and (not state["truncated"])
    return StringCensus(walks, "complete" if complete else "budget",
                        max_length, len(walks), bool(bands))


def find_bands(A, max_length=8):
    """Canonical band walks (cyclic strings) of length <= ``max_length``.

    A band is a cyclic reduced walk of length >= 1, not a proper power, containing
    BOTH a direct and an inverse letter (a pure directed / pure inverse cycle is not
    a band), every rotation a valid walk. Non-empty => ``A`` is rep-INFINITE."""
    _require_string(A)
    rs = reduction_system_of(A)
    Q = A.quiver
    alphabet = [(a, +1) for a in Q.arrows] + [(a, -1) for a in Q.arrows]
    bands, seen = [], set()

    def closes(walk):
        return (letter_target(Q, walk[-1]) == letter_source(Q, walk[0])
                and walk[0] != invert(walk[-1])
                and is_valid_walk(A, walk + (walk[0],), rs))          # wrap pair ok

    def is_band(walk):
        if {d for _, d in walk} != {+1, -1}:            # pure directed/inverse cycle
            return False
        if _is_proper_power(tuple(walk)):
            return False
        return all(is_valid_walk(A, rot, rs) for rot in _rotations(tuple(walk)))

    def grow(walk):
        if len(walk) >= 1 and closes(walk) and is_band(walk):
            cw = _canonical_cyclic(walk)
            if cw not in seen:
                seen.add(cw)
                bands.append(cw)
        if len(walk) >= max_length:
            return
        end = letter_target(Q, walk[-1])
        for ell in alphabet:
            if letter_source(Q, ell) == end and is_valid_walk(A, walk + (ell,), rs):
                grow(walk + (ell,))

    for ell in alphabet:
        grow((ell,))
    return bands
