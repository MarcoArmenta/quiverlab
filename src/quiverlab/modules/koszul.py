"""Koszulity, quadraticity, the quadratic dual A^!, and the Fröberg identity
(Plan 27, Tier 2).

This module holds the *combinatorial* Koszul surface that the Yoneda-algebra
engine (``modules/ext_algebra.py``) integrates into its three-valued verdict:

- ``is_quadratic`` / ``g_quadratic_certificate`` -- quadraticity of the ideal vs
  the Priddy/PBW certificate that A is Koszul,
- ``quadratic_dual`` -- the Koszul dual A^! = kQ^op/(R^perp) (Polishchuk-Positselski
  convention), a genuine ``Algebra`` when finite-dimensional and a lightweight
  :class:`QuadraticDual` record (which still knows its graded dimensions) otherwise,
- ``froberg_obstruction`` / ``dual_dims_crosscheck`` -- the numeric Koszulity
  falsifier ``P(t)*C_A(-t) = I`` and the ``E(A) ~= (A^!)^op`` cross-check.

Everything here is *field-free where the field does not enter*: graded path counts
are integers (they count irreducible monomials by length), and the only linear
algebra over the domain is the orthogonal complement R^perp, done exactly.

**The locked convention (Plan 27).** For Koszul A the Yoneda algebra is
``E(A) = Ext_A^*(A/J, A/J) ~= (A^!)^op`` -- the quadratic dual lives on Q^op and the
op is a real consequence of (right modules) + (left-to-right product). The pairing
uses the repo's left-to-right law: for arrows ``alpha: i->j``, ``beta: j->k`` the
length-2 Q-path ``alpha*beta`` pairs with the length-2 Q^op-path ``beta^o*alpha^o``
(the reversed word -- exactly ``modules/opposite.py``'s reversal). In the two
mutually-reversed bases the pairing Gram matrix is the identity, so R^perp is the
ordinary orthogonal complement computed corner by corner, then relabelled on Q^op.

References: Priddy 1970 (``priddy`` -- G-quadratic/PBW ⇒ Koszul), Fröberg 1999
(``froberg_koszul`` -- the Hilbert-series criterion), Polishchuk-Positselski 2005
(``polishchuk_positselski`` -- quadratic-dual conventions).
"""
from dataclasses import dataclass

from quiverlab.combinat.quiver import Quiver
from quiverlab.combinat.relations import Relation
from quiverlab.core.algebra import Algebra
from quiverlab.errors import (
    AdmissibilityError, FieldError, NotFiniteDimensionalError, QuiverlabError,
)
from quiverlab.fields.linalg import nullspace
from quiverlab.groebner.certificate import certified_irreducibles, check_degree_bound
from quiverlab.groebner.complete import complete
from quiverlab.groebner.order import path_order
from quiverlab.groebner.reduction import reduce_comb, rule_from_comb


# ---------------------------------------------------------------------------
# Guards
# ---------------------------------------------------------------------------
def _require_presentation(A, what):
    """Every Koszul routine needs the bound-quiver presentation kQ/I. A
    structure-constant ``Algebra`` carries no quiver/path basis, so we refuse with
    the same QuiverlabError shape the rest of the module layer uses."""
    if A.quiver is None or A.basis_labels is None:
        raise QuiverlabError(
            f"{what} needs the quiver presentation kQ/I",
            hint="construct the algebra via Quiver.algebra(...); structure-constant "
                 "algebras carry no relations to test for quadraticity")


# ---------------------------------------------------------------------------
# Quadraticity and the G-quadratic (Koszul) certificate
# ---------------------------------------------------------------------------
def is_quadratic(A) -> bool:
    """True iff the defining ideal I of ``A = kQ/I`` is generated in degree 2.

    Decided from the **presentation** ``A.relations`` (assumed minimal, as an
    admissible ``Quiver.algebra`` presentation is): A is quadratic iff every
    defining relation is homogeneous of path-length exactly 2 (all its terms are
    length-2 paths). We read this off the presentation rather than the completed
    reduction system because "quadratic" is a property of the *generators* of I,
    whereas "G-quadratic" (``g_quadratic_certificate``) is the stronger property of
    the *Gröbner basis* (the tips), which may acquire longer words on completion.

    A semisimple algebra (no arrows) and a hereditary path algebra kQ (no relations)
    are quadratic *trivially*: I = 0 is generated in degree 2 vacuously (kQ is the
    quadratic algebra with an empty relation space).
    """
    _require_presentation(A, "is_quadratic")
    # No relations ⇒ hereditary kQ (or semisimple) ⇒ quadratic vacuously.
    # Otherwise every relation must be homogeneous of degree exactly 2: admissibility
    # already forces min_length >= 2, so requiring max_length == 2 pins every term to
    # length 2. An inhomogeneous or cubic minimal relation makes A non-quadratic.
    return all(rel.min_length == 2 and rel.max_length == 2 for rel in (A.relations or []))


def g_quadratic_certificate(A) -> bool:
    """True iff the **confluent** reduction system of A has ALL tips of length 2.

    This is Priddy's PBW / G-quadratic condition: a quadratic Gröbner basis proves
    A is **Koszul** (``priddy``). It CERTIFIES Koszulity; a ``False`` return is
    *inconclusive*, never a disproof -- A may still be Koszul via a different order
    or by non-PBW means. Covers the whole oracle battery (quadratic monomial,
    rad^2 = 0, hereditary, quantum CI, commutative square).

    The reduction system is re-derived by the public ``reduction_system_of`` (the
    same completed system CS consumes). Hereditary/semisimple A have an empty rule
    set, so "all tips length 2" holds vacuously -> True. If the system cannot be
    certified (completion did not stabilize under the degree bound, or the algebra
    is not certified finite-dimensional), we cannot exhibit the PBW basis, so the
    certificate is inconclusive -> ``False``.
    """
    _require_presentation(A, "g_quadratic_certificate")
    from quiverlab.resolutions_cs.build import reduction_system_of
    try:
        rs = reduction_system_of(A)
    except (AdmissibilityError, NotFiniteDimensionalError):
        # No certified confluent PBW basis in hand -> cannot certify Koszul.
        return False
    # Vacuously True on an empty rule set (hereditary kQ / semisimple).
    return all(len(rule.lead) == 2 for rule in rs.rules)


# ---------------------------------------------------------------------------
# The quadratic dual A^!
# ---------------------------------------------------------------------------
def _length2_paths(quiver):
    """All length-2 paths of ``quiver`` as words (a, b) with target(a) == source(b),
    read LEFT TO RIGHT (Assem-Simson-Skowroński)."""
    arrows = quiver.arrows
    return [(a, b)
            for a, (_sa, ta) in arrows.items()
            for b, (sb, _tb) in arrows.items()
            if ta == sb]


def _opposite_quiver(quiver):
    """Q^op: SAME arrow names, reversed endpoints -- matching
    ``modules/opposite.py`` so the dual's labels agree with the repo's A^op."""
    return Quiver(list(quiver.vertices),
                  {name: (t, s) for name, (s, t) in quiver.arrows.items()})


def _to_domain(dom, c):
    """A relation coefficient as an element of ``dom``.

    Relation coefficients are normally the raw parse values (``Fraction`` / sympy
    ``Expr`` / int) that ``dom.coerce`` accepts. But a dual algebra built HERE
    stores coefficients that are ALREADY domain elements (e.g. CC's internal
    ``mpq``, which ``ComplexField.coerce`` rejects); those are valid domain
    elements already, so pass them through unchanged. ``coerce`` is idempotent on
    the values it accepts, so this never silently changes a coercible coefficient.
    """
    try:
        return dom.coerce(c)
    except (FieldError, TypeError, ValueError):
        return c


def _dual_data(A):
    """The quadratic-dual presentation as raw data ``(Qop, combs, dom)``.

    ``combs`` is a list of relation elements (each a dict ``word -> dom-element``
    over Q^op) spanning R^perp corner by corner:

    * Length-2 Q-paths are grouped by corner (i, k) = (source, target).
    * Inside a corner, R_{ik} is spanned by the coordinate vectors of A's quadratic
      relations living there; R^perp_{ik} is its ordinary orthogonal complement
      (Gram = I in the corresponding bases), i.e. ``nullspace`` of the relation
      rows. A corner with NO relations contributes its whole length-2 space to
      R^perp (every path becomes a dual relation -- the hereditary case ⇒ kQ^op/J^2);
      a corner whose relations span everything contributes nothing (rad^2 = 0 ⇒ free).
    * Each R^perp basis vector is relabelled to Q^op by REVERSING every path word
      (``(alpha, beta) -> (beta, alpha)``): the dual relation ``sum v[s]*reverse(B[s])``
      is a parallel Q^op path k -> i.
    """
    _require_presentation(A, "quadratic_dual")
    if not is_quadratic(A):
        raise QuiverlabError(
            "quadratic_dual is defined only for a quadratic algebra (ideal generated "
            "in degree 2)",
            hint="is_quadratic(A) is False: a minimal relation has length != 2, so A^! "
                 "is not the Koszul dual; compute the Yoneda algebra directly instead")
    dom = A.domain
    arrows = A.quiver.arrows

    # Group length-2 Q-paths by corner (i, k).
    corners = {}
    for w in _length2_paths(A.quiver):
        corners.setdefault((arrows[w[0]][0], arrows[w[1]][1]), []).append(w)

    # Relation rows per corner (coordinate vectors over that corner's path basis).
    rel_rows = {}
    for rel in (A.relations or []):
        key = (rel.source, rel.target)
        basis = corners.get(key)
        if basis is None:            # a relation with no length-2 paths cannot arise
            continue                 # for a quadratic ideal, but stay defensive
        index = {w: s for s, w in enumerate(basis)}
        vec = [dom.zero()] * len(basis)
        for c, w in rel.terms:
            vec[index[w]] = _to_domain(dom, c)
        rel_rows.setdefault(key, []).append(vec)

    combs = []
    for (i, k), basis in corners.items():
        rows = rel_rows.get((i, k), [])
        if rows:
            perp = nullspace(rows, dom)
        else:
            # No relations in this corner: R_{ik} = 0, so R^perp_{ik} is the whole
            # space -- every length-2 path is a standard-basis vector of R^perp.
            # (nullspace([]) returns [] -- "zero space" -- which is the WRONG
            # complement here, so we build the full basis explicitly.)
            perp = [[dom.one() if s == t else dom.zero() for t in range(len(basis))]
                    for s in range(len(basis))]
        for v in perp:
            comb = {(basis[s][1], basis[s][0]): v[s]     # reverse the path word
                    for s in range(len(basis)) if not dom.is_zero(v[s])}
            if comb:
                combs.append(comb)
    return _opposite_quiver(A.quiver), combs, dom


def _relation_of_comb(quiver, comb):
    """A ``Relation`` from a dual element ``comb`` (word -> dom-element). Terms carry
    domain-element coefficients (source/target read off Q^op) so the resulting
    algebra is itself re-dualizable."""
    terms = tuple((coeff, word) for word, coeff in comb.items())
    w0 = terms[0][1]
    return Relation(terms, source=quiver.word_source(w0), target=quiver.word_target(w0))


def _completed_rules(quiver, combs, dom, degree_bound):
    """Completed reduction rules of the dual presentation, from domain-element combs.

    We build rules DIRECTLY (``rule_from_comb`` + ``complete``) rather than through
    ``build_reduction_system`` because the coefficients are already domain elements
    (see ``_to_domain``): re-parsing them would be a needless -- and, for CC's
    ``mpq``, impossible -- round-trip."""
    order = path_order(quiver)
    init = [rule_from_comb(comb, order, quiver, dom) for comb in combs]
    return order, complete(init, order, quiver, dom, degree_bound)


def _graded_count_matrices(quiver, tips, through):
    """``mats[n][u][v]`` = number of length-n paths u -> v (vertices in
    ``quiver.vertices`` order) avoiding every tip as a contiguous factor -- the
    graded k-dimensions ``dim e_u A_n e_v`` of the monomial (associated-graded)
    algebra with these tips. Field-free integer combinatorics.

    Built by extending tip-free words one arrow at a time: an irreducible word stays
    irreducible under an extension exactly when the extension does not create a tip
    as a SUFFIX (any earlier occurrence was already forbidden)."""
    verts = quiver.vertices
    idx = {v: i for i, v in enumerate(verts)}
    nv = len(verts)
    tips = [tuple(t) for t in tips]

    m0 = [[0] * nv for _ in range(nv)]
    for v in verts:
        m0[idx[v]][idx[v]] = 1                      # the length-0 path e_v
    mats = [m0]
    frontier = [((), v, v) for v in verts]          # (word, source vertex, current end)
    for _n in range(1, through + 1):
        mn = [[0] * nv for _ in range(nv)]
        new_frontier = []
        for word, src, cur in frontier:
            for a, (a_s, a_t) in quiver.arrows.items():
                if a_s != cur:
                    continue
                w2 = word + (a,)
                if any(len(t) <= len(w2) and w2[-len(t):] == t for t in tips):
                    continue
                new_frontier.append((w2, src, a_t))
                mn[idx[src]][idx[a_t]] += 1
        mats.append(mn)
        frontier = new_frontier
    return mats


@dataclass(frozen=True)
class QuadraticDual:
    """A lightweight record for A^! when it is (or may be) infinite-dimensional as an
    algebra -- e.g. ``A^!(k[x]/x^2) = k[y]`` or the quantum plane. It knows its
    quiver, its quadratic relations, and (combinatorially, via the completed
    reduction system of the dual presentation) its graded dimensions, which is all
    the Fröberg / cross-check consumers need. When A^! is finite-dimensional,
    ``quadratic_dual`` returns a genuine ``Algebra`` instead."""
    quiver: object
    domain: object
    _combs: tuple      # each element a dict Q^op-word -> domain element

    @property
    def relations(self):
        """The quadratic relations R^perp of A^! (parallel Q^op paths)."""
        return [_relation_of_comb(self.quiver, dict(c)) for c in self._combs]

    def graded_dims_matrix(self, through):
        """``mats[n][u][v]`` = dim of the corner (u, v) piece of ``(A^!)_n`` (vertices
        in ``quiver.vertices`` order)."""
        # Complete generously so every tip of length <= through is found (an overlap
        # of length-<=through leads reaches length 2*through-1).
        _order, rules = _completed_rules(self.quiver, list(self._combs), self.domain,
                                         max(2 * through, 8))
        return _graded_count_matrices(self.quiver, [r.lead for r in rules], through)

    def graded_dims(self, through):
        """Total graded dimensions ``[dim (A^!)_0, ..., dim (A^!)_through]``."""
        mats = self.graded_dims_matrix(through)
        return [sum(sum(row) for row in mn) for mn in mats]

    def __repr__(self):
        rels = self.relations
        head = (f"QuadraticDual A^! over {self.domain.name} "
                f"(likely infinite-dimensional as an algebra)")
        body = ["  vertices: " + ", ".join(str(v) for v in self.quiver.vertices),
                "  arrows: " + "; ".join(f"{n}: {s} -> {t}"
                                         for n, (s, t) in self.quiver.arrows.items())]
        if rels:
            body.append("  relations R^perp: " + "; ".join(repr(r) for r in rels))
        else:
            body.append("  relations R^perp: (none -- free path algebra kQ^op)")
        return "\n".join([head] + body)


def _finite_dual_algebra(quiver, combs, dom):
    """Build A^! as a genuine ``Algebra`` from the dual presentation, or raise
    ``NotFiniteDimensionalError`` (infinite) / ``AdmissibilityError`` (completion did
    not stabilize under the certificate's degree bound). Mirrors
    ``groebner.lower.groebner_algebra`` but consumes domain-element combs directly."""
    # Quadratic init (leads of length 2): the certificate's default degree bound is
    # max(8, 2*2+4) = 8. If completion outgrows it, check_degree_bound raises loudly.
    degree_bound = 8
    order, rules = _completed_rules(quiver, combs, dom, degree_bound)
    check_degree_bound(rules, degree_bound)
    irreducibles = certified_irreducibles(quiver, rules)   # loud if infinite

    zero, one = dom.zero(), dom.one()
    basis = [("e", v) for v in quiver.vertices] + [("p", w) for w in irreducibles]
    index = {b: i for i, b in enumerate(basis)}
    m = len(basis)

    def src(b):
        return b[1] if b[0] == "e" else quiver.word_source(b[1])

    def tgt(b):
        return b[1] if b[0] == "e" else quiver.word_target(b[1])

    def product_vector(bi, bj):
        vec = [zero] * m
        if tgt(bi) != src(bj):
            return vec
        if bi[0] == "e":
            vec[index[bj]] = one
            return vec
        if bj[0] == "e":
            vec[index[bi]] = one
            return vec
        nf = reduce_comb({bi[1] + bj[1]: one}, rules, order, dom)
        for word, coeff in nf.items():
            vec[index[("p", word)]] = coeff
        return vec

    T = [[product_vector(bi, bj) for bj in basis] for bi in basis]
    unit = [zero] * m
    for v in quiver.vertices:
        unit[index[("e", v)]] = one
    labels = [f"e_{b[1]}" if b[0] == "e" else "*".join(b[1]) for b in basis]
    relations = [_relation_of_comb(quiver, comb) for comb in combs]
    return Algebra(dom, T, unit, basis_labels=labels, _quiver=quiver,
                   _relations=relations)


def quadratic_dual(A):
    """The quadratic (Koszul) dual A^! = kQ^op/(R^perp) of a quadratic ``A = kQ/I``.

    Returns a genuine ``Algebra`` when A^! is finite-dimensional (admissible), and a
    :class:`QuadraticDual` record otherwise -- A^! is frequently infinite-dimensional
    (``A^!(k[x]/x^2) = k[y]``; the QuantumCI dual is the quantum plane). Both support
    graded dimensions, which is what the Fröberg / cross-check consumers use.

    Raises ``QuiverlabError`` for non-quadratic A (then A^! is not the Koszul dual).

    Convention (Polishchuk-Positselski, ``polishchuk_positselski``): Q^op with the
    orthogonal complement R^perp under the pairing that sends a length-2 Q-path to
    the reversed length-2 Q^op-path (Gram = I in the reversed bases). MIND the
    reversal: ``dual(k[x]/x^2)`` is the FREE loop k[y] (R = span{x^2}, R^perp = 0),
    while ``dual(kQ)`` is ``kQ^op/J^2`` (R = 0, R^perp = all length-2 paths).
    """
    Qop, combs, dom = _dual_data(A)
    try:
        return _finite_dual_algebra(Qop, combs, dom)
    except (NotFiniteDimensionalError, AdmissibilityError):
        # A^! is infinite-dimensional (or its finiteness is not certifiable under the
        # v1 degree bound): return the record, which still delivers graded_dims.
        return QuadraticDual(Qop, dom, tuple(combs))


# ---------------------------------------------------------------------------
# The Fröberg numeric Koszulity test  P(t) * C_A(-t) = I
# ---------------------------------------------------------------------------
def _is_length_graded(A):
    """A is length-graded (I is a homogeneous ideal) iff every defining relation is
    length-homogeneous. Then #irreducible paths of length n = dim A_n exactly (the
    length-lex Gröbner tips are homogeneous)."""
    return all(rel.min_length == rel.max_length for rel in (A.relations or []))


def _algebra_graded_matrices(A, through):
    """``C_b[i][j] = dim e_i A_b e_j`` for b = 0..through (vertices in
    ``A.quiver.vertices`` order), from A's own reduction-system tips."""
    from quiverlab.resolutions_cs.build import reduction_system_of
    tips = reduction_system_of(A).leading_words()
    return _graded_count_matrices(A.quiver, tips, through)


def _matmul(P, C, nv):
    """Integer matrix product P @ C (nv x nv)."""
    out = [[0] * nv for _ in range(nv)]
    for i in range(nv):
        Pi = P[i]
        for t in range(nv):
            pit = Pi[t]
            if pit == 0:
                continue
            Ct = C[t]
            for j in range(nv):
                out[i][j] += pit * Ct[j]
    return out


def _check_hilbert_shape(hilbert_matrix, nv, d, what):
    """hilbert_matrix must supply an nv x nv integer matrix for every degree 0..d."""
    if len(hilbert_matrix) < d + 1:
        raise QuiverlabError(
            f"{what} needs the Ext dimension matrix through degree {d}: got "
            f"{len(hilbert_matrix)} degrees, need {d + 1}",
            hint="pass hilbert_matrix[n] for n = 0..d, each a |Q0|x|Q0| int matrix")
    for n, mat in enumerate(hilbert_matrix[:d + 1]):
        if len(mat) != nv or any(len(row) != nv for row in mat):
            raise QuiverlabError(
                f"{what}: hilbert_matrix[{n}] is not {nv}x{nv} (|Q0| = {nv})",
                hint="corners are indexed by A.quiver.vertices order")


def froberg_obstruction(A, hilbert_matrix, d):
    """First degree n <= d where the Fröberg identity ``P(t) * C_A(-t) = I`` FAILS
    coefficientwise, else ``None`` (Fröberg, ``froberg_koszul``).

    ``hilbert_matrix[n][i][j] = dim E^n_{ij}`` (the graded Betti / Ext^n(S_i, S_j)
    dimensions -- the Yoneda algebra's Hilbert matrices), a list of |Q0|x|Q0| int
    matrices for n = 0..d. ``C_A(t)_{ij} = dim e_i A_n e_j`` is A's graded
    path-count matrix. Writing ``P = sum P_n t^n`` and ``C_A(-t) = sum (-1)^n C_n t^n``,
    the product's degree-N matrix is ``M_N = sum_{a+b=N} (-1)^b P_a @ C_b``, which the
    identity forces to be I at N = 0 and 0 for N >= 1. We return the first N where
    ``M_N`` differs, else ``None`` (no obstruction through d -- necessary for Koszul,
    not sufficient; a numeric pass is per-field only).

    A must be length-graded (homogeneous I) for ``C_A`` to mean dim A_n; if it is
    not, we return ``None`` (the length-associated-graded is available but the
    identity is not meaningful off the graded case -- documented honest limitation).

    Hand-check on kA_2: C = [[1, t], [0, 1]], P = [[1, t], [0, 1]], and
    P*C(-t) = [[1, 0], [0, 1]] = I -> no obstruction.
    """
    _require_presentation(A, "froberg_obstruction")
    if not _is_length_graded(A):
        return None
    nv = len(A.quiver.vertices)
    _check_hilbert_shape(hilbert_matrix, nv, d, "froberg_obstruction")
    C = _algebra_graded_matrices(A, d)     # C[b], b = 0..d
    for N in range(d + 1):
        MN = [[0] * nv for _ in range(nv)]
        for a in range(N + 1):
            b = N - a
            prod = _matmul(hilbert_matrix[a], C[b], nv)
            sign = -1 if b % 2 else 1
            for i in range(nv):
                for j in range(nv):
                    MN[i][j] += sign * prod[i][j]
        expected = [[1 if (N == 0 and i == j) else 0 for j in range(nv)]
                    for i in range(nv)]
        if MN != expected:
            return N
    return None


def dual_dims_crosscheck(A, hilbert_matrix, d):
    """For quadratic A, compare the given Ext dimensions against the graded dimensions
    of ``(A^!)^op`` through degree d: returns ``True``/``False``, or ``None`` when not
    computable (A not quadratic, or too few degrees supplied).

    ``E(A) ~= (A^!)^op`` (the locked Plan-27 convention), so ``dim E^n_{ij}`` must
    equal ``dim (A^!)^op`` in corner (i, j) = ``dim A^!`` in corner (j, i) -- the op
    TRANSPOSES corners. Concretely: ``hilbert_matrix[n][i][j] == D_n[j][i]`` where
    ``D_n[u][v]`` counts length-n Q^op paths u -> v. Getting the transpose right is
    the whole point: on the rad^2 = 0 line 1 -> 2 -> 3, ``E = kQ`` (NOT kQ^op), and
    the identification only closes after the corner flip.
    """
    _require_presentation(A, "dual_dims_crosscheck")
    if not is_quadratic(A):
        return None
    nv = len(A.quiver.vertices)
    if len(hilbert_matrix) < d + 1:
        return None
    Qop, combs, dom = _dual_data(A)
    _order, rules = _completed_rules(Qop, combs, dom, max(2 * d, 8))
    D = _graded_count_matrices(Qop, [r.lead for r in rules], d)   # D_n[u][v]
    for n in range(d + 1):
        for i in range(nv):
            for j in range(nv):
                if hilbert_matrix[n][i][j] != D[n][j][i]:          # op: transpose
                    return False
    return True
