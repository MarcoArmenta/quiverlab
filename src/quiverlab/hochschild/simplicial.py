"""The order complex of a finite poset, and its exact simplicial (co)homology --
quiverlab's SECOND, independent route to ``HH^*`` of an incidence algebra (Plan 75 / R9).

**The theorem.** For a finite poset ``P`` and any field ``k`` (Cibils 1989, generalizing
Gerstenhaber-Schack 1983 from the face poset of a simplicial complex to an ARBITRARY
finite poset)::

    HH^n(kP)  =  H^n(Delta(P); k)      for all n >= 0,

as graded rings. ``Delta(P)`` is the order complex (nerve): its ``p``-simplices are the
CHAINS ``x_0 < x_1 < ... < x_p``.

**Why this is the fast path.** The general route computes over the enveloping algebra
``(kP)^e`` -- a resolution term of size ``O(m^2)`` for ``m = dim kP``, rebuilt per prime.
This route computes the *combinatorial* cochain complex of ``Delta(P)``, whose degree-``p``
dimension is just the number of ``p``-chains. Better still, the boundary matrices are
INTEGER matrices computed ONCE: Smith normal form over ZZ then yields ``H^n(Delta; k)`` for
EVERY field simultaneously by universal coefficients::

    dim_k H^n(Delta; k)  =  rank_ZZ H_n  +  t_n(char k)  +  t_{n-1}(char k),

where ``t_n(p)`` counts the invariant factors of ``H_n(Delta; ZZ)`` divisible by ``p``.
**Torsion in integral homology is exactly the source of characteristic-dependent
``HH^*``** -- the ``RP^2`` pin (``H_1(RP^2; ZZ) = Z/2``) has ``HH^2`` of dimension 1 over
``GF(2)`` but 0 over ``QQ``, and this engine explains WHY from one integer computation,
which the general engine (a separate run per field) cannot.

Everything stays in exact integer arithmetic until the final reduction into the field.
"""
from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import rank as _rank

_SIMPLICIAL_CITATIONS = ("cibils_incidence", "gerstenhaber_schack_simplicial")


class OrderComplex:
    """The order complex ``Delta(P)`` as chains-by-dimension plus integer boundary maps.

    ``chains[p]`` is the ordered list of ``p``-simplices, each a tuple ``(x_0, ..., x_p)``
    STRICTLY increasing in ``P``. ``face_vector()`` is ``(len(chains[0]), len(chains[1]),
    ...)`` -- note the standard convention: entry ``p`` counts ``p``-SIMPLICES, i.e. chains
    with ``p + 1`` elements.
    """

    def __init__(self, chains, poset=None):
        self.chains = [list(c) for c in chains]
        self.poset = poset
        self._index = [{tuple(s): i for i, s in enumerate(level)} for level in self.chains]

    # ------------------------------------------------------------------ build
    @classmethod
    def of(cls, poset):
        """Enumerate the chains of ``poset`` by dimension.

        A chain is a subset totally ordered by ``<=``; it is emitted in its unique
        increasing order. Built incrementally (extend each ``p``-chain by an element
        strictly above its top), which is why no subset is ever generated twice and the
        output of each level is deterministic.
        """
        elems = list(poset.elements)
        strict_up = {x: [y for y in elems if x != y and poset.leq(x, y)] for x in elems}
        levels = [[(x,) for x in elems]]
        while levels[-1]:
            nxt = []
            for ch in levels[-1]:
                for y in strict_up[ch[-1]]:
                    nxt.append(ch + (y,))
            if not nxt:
                break
            levels.append(nxt)
        return cls(levels, poset=poset)

    # ------------------------------------------------------------------ shape
    def face_vector(self):
        """``(#0-simplices, #1-simplices, ...)`` -- entry ``p`` counts chains of LENGTH
        ``p`` (i.e. ``p + 1`` elements)."""
        return tuple(len(level) for level in self.chains)

    @property
    def top_dimension(self):
        return len(self.chains) - 1

    def euler_characteristic(self):
        """``sum (-1)^p f_p`` -- an integer invariant, independent of any field."""
        return sum((-1) ** p * len(level) for p, level in enumerate(self.chains))

    # ------------------------------------------------------------- boundaries
    def boundary(self, k):
        """The INTEGER matrix of ``d_k : C_k -> C_{k-1}``, the alternating face map
        ``d(x_0..x_k) = sum_j (-1)^j (x_0..^x_j..x_k)``.

        Rows are indexed by ``(k-1)``-simplices, columns by ``k``-simplices (so it is a
        map of column vectors, matching the rest of quiverlab). Returns ``[]`` when either
        side is empty -- callers treat an empty matrix as rank 0.
        """
        if k <= 0 or k > self.top_dimension:
            return []
        rows, cols = len(self.chains[k - 1]), len(self.chains[k])
        if rows == 0 or cols == 0:
            return []
        M = [[0] * cols for _ in range(rows)]
        idx = self._index[k - 1]
        for c, simplex in enumerate(self.chains[k]):
            for j in range(len(simplex)):
                face = simplex[:j] + simplex[j + 1:]
                M[idx[face]][c] += (-1) ** j
        return M


# --------------------------------------------------------------------------- #
# field cohomology (the dimension answer)
# --------------------------------------------------------------------------- #
def _rank_over(M, dom):
    if not M or not M[0]:
        return 0
    return _rank([[dom.coerce(x) for x in row] for row in M], dom)


def simplicial_cohomology_dims(oc, top, domain):
    """``[dim H^0, ..., dim H^top]`` of the order complex over the exact ``domain``.

    Over a FIELD the cochain complex is the transpose of the chain complex, so cohomology
    and homology have equal dimensions degreewise; we compute the homology form

        dim H_n  =  dim C_n  -  rank d_n  -  rank d_{n+1}

    which needs only ranks (no transposes, no bases). Degrees past the top dimension of
    the complex are ``0`` and are reported as such rather than truncated, so the returned
    list always has length ``top + 1``.
    """
    ranks = {}

    def r(k):
        if k not in ranks:
            ranks[k] = _rank_over(oc.boundary(k), domain)
        return ranks[k]

    out = []
    for n in range(top + 1):
        cn = len(oc.chains[n]) if n <= oc.top_dimension else 0
        out.append(cn - r(n) - r(n + 1))
    return out


# --------------------------------------------------------------------------- #
# integral homology by Smith normal form (the characteristic certificate)
# --------------------------------------------------------------------------- #
def _smith_invariant_factors(M):
    """The nonzero invariant factors of an integer matrix, by Smith normal form.

    Plain exact-integer SNF: repeatedly pick the smallest nonzero pivot, clear its row and
    column by integer row/column operations, recurse. No fractions ever appear, so the
    result is exact for any size this engine will meet. Returns the positive invariant
    factors in order; its length is ``rank M``.
    """
    A = [row[:] for row in M]
    if not A or not A[0]:
        return []
    rows, cols = len(A), len(A[0])
    factors = []
    t = 0
    while t < rows and t < cols:
        # find a pivot: the nonzero entry of least absolute value in the active block
        piv = None
        for i in range(t, rows):
            for j in range(t, cols):
                if A[i][j] and (piv is None or abs(A[i][j]) < abs(A[piv[0]][piv[1]])):
                    piv = (i, j)
        if piv is None:
            break
        A[t], A[piv[0]] = A[piv[0]], A[t]
        for row in A:
            row[t], row[piv[1]] = row[piv[1]], row[t]
        # clear the pivot row and column; repeat until both are clear (the pivot may
        # shrink on the way, which is exactly what makes this terminate)
        while True:
            dirty = False
            for i in range(t + 1, rows):
                if A[i][t]:
                    q = A[i][t] // A[t][t]
                    for j in range(t, cols):
                        A[i][j] -= q * A[t][j]
                    if A[i][t]:
                        A[t], A[i] = A[i], A[t]
                        dirty = True
            for j in range(t + 1, cols):
                if A[t][j]:
                    q = A[t][j] // A[t][t]
                    for i in range(t, rows):
                        A[i][j] -= q * A[i][t]
                    if A[t][j]:
                        for row in A:
                            row[t], row[j] = row[j], row[t]
                        dirty = True
            if not dirty:
                break
        factors.append(abs(A[t][t]))
        t += 1
    return factors


def integral_homology(oc, top):
    """``[(free_rank, invariant_factors), ...]`` for ``H_0 .. H_top`` over ``ZZ``.

    ``H_n = ZZ^{free_rank} (+) ZZ/f_1 (+) ... `` where the ``f_i > 1`` are the torsion
    invariant factors. Computed from the Smith normal forms of the integer boundary maps:
    the invariant factors of ``d_{n+1}`` give the torsion of ``H_n``, and
    ``free_rank = dim C_n - rank d_n - rank d_{n+1}``.

    This is the CHARACTERISTIC CERTIFICATE: a factor divisible by ``p`` is precisely why
    ``H^*(Delta; GF(p))`` exceeds ``H^*(Delta; QQ)``.
    """
    out = []
    for n in range(top + 1):
        cn = len(oc.chains[n]) if n <= oc.top_dimension else 0
        f_n = _smith_invariant_factors(oc.boundary(n))
        f_n1 = _smith_invariant_factors(oc.boundary(n + 1))
        free = cn - len(f_n) - len(f_n1)
        torsion = tuple(sorted(f for f in f_n1 if f > 1))
        out.append((free, torsion))
    return out


def cohomology_dims_from_integral(ih, top, characteristic):
    """Universal coefficients: field cohomology dims from the INTEGRAL homology alone.

        dim_k H^n  =  free_n  +  t_n(p)  +  t_{n-1}(p),

    with ``t_n(p)`` the number of invariant factors of ``H_n`` divisible by ``p``
    (``t_n(0) = 0`` in characteristic 0). This is the payoff of the integer route: ONE
    Smith normal form answers EVERY characteristic.
    """
    def t(n):
        if n < 0 or n >= len(ih) or characteristic == 0:
            return 0
        return sum(1 for f in ih[n][1] if f % characteristic == 0)

    out = []
    for n in range(top + 1):
        free = ih[n][0] if n < len(ih) else 0
        out.append(free + t(n) + t(n - 1))
    return out


def order_complex_of(algebra):
    """The order complex of an algebra KNOWN to be an incidence algebra, via the
    ``_poset`` provenance stashed by :func:`quiverlab.families.IncidenceAlgebra`.

    Refuses LOUDLY otherwise. Recognizing an arbitrary ``kQ/I`` as an incidence algebra is
    a genuinely harder problem (up to isomorphism it is not even obviously decidable here),
    and guessing would silently apply a theorem whose hypothesis was never checked -- so
    the provenance is required, not inferred.
    """
    P = getattr(algebra, "_poset", None)
    if P is None:
        raise QuiverlabError(
            "this algebra carries no poset provenance, so the order-complex route is "
            "unavailable: HH^*(kP) = H^*(Delta(P)) needs to KNOW that A is the incidence "
            "algebra of a poset, and quiverlab never GUESSES that from a presentation",
            hint="build it with quiverlab.families.IncidenceAlgebra(covers, field=...), "
                 "or use the general engines (engine='auto'/'cs'/'bar')")
    return OrderComplex.of(P)
