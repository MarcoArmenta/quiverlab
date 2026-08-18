"""The combinatorial Tits quadratic form and the tame/wild certificate
(Plan 62 / R19).

The combinatorial Tits form of a triangular algebra A = kQ/I is

    q_A(x) = sum_{i in Q_0} x_i^2 - sum_{(i->j) in Q_1} x_i x_j
                                  + sum_{i,j} r_{ij} x_i x_j,

with r_{ij} = the number of MINIMAL relations from i to j
(= dim_k e_j (I/(rad.I + I.rad)) e_i, the P56 count; independent of the chosen
minimal generating set), which equals dim_k Ext^2_A(S_i, S_j). It truncates the
homological Euler form at Ext^2 and is defined for EVERY admissible presentation.
The computation of q_A and of its weak positivity / weak nonnegativity is
FIELD-FREE integer arithmetic.

The representation-type verdict on top of the form requires an algebraically
closed base field (`_is_alg_closed`: the CC working domain, or -- absent P61's
`is_algebraically_closed` field flag -- any characteristic-0 field read by base
change to the algebraic closure, the form and the representation type being
field-independent in characteristic 0). It is gated on the P56 certificate:

    * A representation-finite  <=>  q_A weakly positive     (Bongartz 1984,
      simply connected);
    * A tame                   <=>  q_A weakly nonnegative  (Brustle-de la
      Pena-Skowronski 2011, STRONGLY simply connected).

giving rep-finite subset tame subset (weakly nonnegative); wild = not weakly
nonnegative. Every "no" carries an EXACT witness d >= 0 with q_A(d) <= 0 (isotropic,
tame boundary) or q_A(d) < 0 (wild). P56's three-valued None propagates to a None
verdict -- never a fabricated tame/wild.

Distinct from Plan 38's `forms.tits_form` (the homological Euler form via C^{-1},
definiteness over all of R^n). This module's public method is
`Algebra.tits_form_combinatorial` (P38's `Algebra.tits_form` is left untouched);
the two coincide iff gl.dim A <= 2.

Refs: Bongartz Math. Ann. 269 (1984) 1-12; Brustle-de la Pena-Skowronski Adv.
Math. 226 (2011) 887-951; Ovsienko (1978); von Hohne Proc. LMS 73 (1996) 47-67;
BJP (Springer AA 25, 2019); Encyclopedia of Mathematics "Tits quadratic form".
"""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError

_BOX = 6                    # Ovsienko's bound: positive roots have entries <= 6


# ---------------------------------------------------------------------------
# the unit form value object
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class UnitForm:
    """An integer symmetric unit form: n variables, symmetric integer Gram matrix
    `gram` with diagonal 2 (q(d) = (1/2) d^T gram d), and `labels` = the vertex
    labels in index order (for witness readout)."""
    n: int
    gram: tuple
    labels: tuple

    def evaluate(self, d) -> int:
        """q(d) = (1/2) d^T gram d, exact int. `d` is a sequence in index order."""
        g = self.gram
        n = self.n
        dd = list(d)
        return sum(g[i][j] * dd[i] * dd[j] for i in range(n) for j in range(n)) // 2

    def restrict(self, support) -> "UnitForm":
        """Principal submatrix on the index list `support` (a subset of range(n))."""
        idx = list(support)
        sub = tuple(tuple(self.gram[i][j] for j in idx) for i in idx)
        return UnitForm(len(idx), sub, tuple(self.labels[i] for i in idx))

    def underlying_graph(self) -> dict:
        """Adjacency: i ~ j iff ``gram[i][j] != 0`` (i != j)."""
        g = self.gram
        n = self.n
        return {i: {j for j in range(n) if j != i and g[i][j] != 0} for i in range(n)}


# ---------------------------------------------------------------------------
# building the form from an algebra
# ---------------------------------------------------------------------------
def _require_quiver(A, what):
    if A.quiver is None:
        raise QuiverlabError(
            f"{what} needs the quiver presentation",
            hint="build the algebra via Quiver.algebra(...); structure-constant "
                 "algebras carry no quiver, so the combinatorial Tits form is undefined")


def _arrow_counts(quiver):
    counts: dict = {}
    for (s, t) in quiver.arrows.values():
        counts[(s, t)] = counts.get((s, t), 0) + 1
    return counts


def tits_form_combinatorial(A, d) -> int:
    """q_A(d) = sum d_i^2 - sum_{arrows i->j} d_i d_j + sum_{(i,j)} r_ij d_i d_j,
    with r_ij = minimal_relation_counts(A). `d` is a vertex-order sequence or a
    {vertex: value} dict. Field-free exact int. Loud on presentation-less input."""
    _require_quiver(A, "tits_form_combinatorial")
    from quiverlab.invariants.coverings import minimal_relation_counts
    verts = list(A.quiver.vertices)
    if isinstance(d, dict):
        val = {v: int(d.get(v, 0)) for v in verts}
    else:
        dd = list(d)
        if len(dd) != len(verts):
            raise QuiverlabError(
                f"tits_form_combinatorial: dimension vector has {len(dd)} entries "
                f"but the quiver has {len(verts)} vertices",
                hint="pass a value per vertex (vertex order) or a {vertex: value} dict")
        val = {v: int(dd[i]) for i, v in enumerate(verts)}
    total = sum(val[v] * val[v] for v in verts)
    for (s, t) in A.quiver.arrows.values():
        total -= val[s] * val[t]
    for (s, t), r in minimal_relation_counts(A).items():
        total += r * val[s] * val[t]
    return total


def tits_matrix_combinatorial(A):
    """The integer symmetric Gram matrix G with q_A(d) = (1/2) d^T G d, as a sympy
    Matrix (for display / the payload):
        G_ii = 2 - 2 a_ii + 2 r_ii ;  G_ij = r_ij + r_ji - a_ij - a_ji  (i != j),
    with a = arrow counts, r = minimal-relation counts. For a triangular (acyclic)
    quiver G_ii = 2 automatically. Loud on presentation-less input."""
    import sympy as sp
    _require_quiver(A, "tits_matrix_combinatorial")
    from quiverlab.invariants.coverings import minimal_relation_counts
    verts = list(A.quiver.vertices)
    pos = {v: i for i, v in enumerate(verts)}
    n = len(verts)
    a = _arrow_counts(A.quiver)
    r = minimal_relation_counts(A)
    G = [[0] * n for _ in range(n)]
    for i, v in enumerate(verts):
        G[i][i] = 2 - 2 * a.get((v, v), 0) + 2 * r.get((v, v), 0)
    for (s, t), c in a.items():
        if s != t:
            G[pos[s]][pos[t]] -= c
            G[pos[t]][pos[s]] -= c
    for (s, t), c in r.items():
        if s != t:
            G[pos[s]][pos[t]] += c
            G[pos[t]][pos[s]] += c
    return sp.Matrix(G)


def as_unit_form(A) -> UnitForm:
    """The combinatorial Tits form of A as a UnitForm. RAISES loudly unless
    A.quiver.is_acyclic() -- the Tits-form tame/wild theory is defined for
    TRIANGULAR algebras only (M2: the G_ii == 2 heuristic is unsound, it wrongly
    passes k[x]/(x^3), so is_acyclic() is the honest gate)."""
    _require_quiver(A, "as_unit_form")
    if not A.quiver.is_acyclic():
        raise QuiverlabError(
            "the combinatorial Tits form is a unit form only for a triangular "
            "(acyclic) quiver; this quiver has an oriented cycle/loop",
            hint="the Bongartz / Brustle-de la Pena-Skowronski tame/wild theory is "
                 "stated for triangular algebras (Skowronski 1993)")
    from quiverlab.invariants.coverings import minimal_relation_counts
    verts = list(A.quiver.vertices)
    pos = {v: i for i, v in enumerate(verts)}
    n = len(verts)
    a = _arrow_counts(A.quiver)
    r = minimal_relation_counts(A)
    G = [[0] * n for _ in range(n)]
    for i in range(n):
        G[i][i] = 2                          # acyclic: no loops, no i->i relations
    for (s, t), c in a.items():
        G[pos[s]][pos[t]] -= c
        G[pos[t]][pos[s]] -= c
    for (s, t), c in r.items():
        G[pos[s]][pos[t]] += c
        G[pos[t]][pos[s]] += c
    gram = tuple(tuple(row) for row in G)
    return UnitForm(n, gram, tuple(verts))


# ---------------------------------------------------------------------------
# form verdicts
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class FormVerdict:
    holds: object            # True / False (witness found) / None (undecided)
    witness: object          # on False: the exact d >= 0 (index-order tuple)
    witness_value: object    # q(witness)
    reason: str
    checked: int             # lattice points evaluated (0 for a pure list/PSD decision)


def _sympy_gram(form):
    import sympy as sp
    return sp.Matrix([[int(x) for x in row] for row in form.gram])


def _is_pd(form) -> bool:
    return bool(_sympy_gram(form).is_positive_definite)


def _is_psd(form) -> bool:
    return bool(_sympy_gram(form).is_positive_semidefinite)


def _radical_witness(form):
    """A sign-definite integer generator of the radical (kernel of the Gram
    matrix), normalised to be nonnegative and primitive, or None. Such a vector z
    has q(z) = 0 -- an exact isotropic witness on the nonnegative cone."""
    import sympy as sp
    M = _sympy_gram(form)
    ns = M.nullspace()
    for v in ns:
        denom = sp.ilcm(*[c.q for c in v]) if any(c.q != 1 for c in v) else 1
        z = [int(c * denom) for c in v]
        if all(x <= 0 for x in z):
            z = [-x for x in z]
        if all(x >= 0 for x in z) and any(x > 0 for x in z):
            g = 0
            for x in z:
                g = sp.igcd(g, x)
            if g and g > 1:
                z = [x // g for x in z]
            return tuple(z)
    return None


class _Budget(Exception):
    pass


def _connected_supports(graph, n):
    """Yield connected index-subsets (frozensets) of the underlying graph, smallest
    first, deduplicated. Includes the singletons."""
    seen = set()
    frontier = []
    for v in range(n):
        s = frozenset({v})
        if s not in seen:
            seen.add(s)
            frontier.append(s)
            yield s
    while frontier:
        nxt = []
        for s in frontier:
            border = set()
            for v in s:
                border |= graph.get(v, set())
            for w in border - s:
                t = s | {w}
                if t not in seen:
                    seen.add(t)
                    nxt.append(t)
                    yield t
        frontier = nxt


def _search_support(gram, support, strict, state):
    """Search sincere d in {1..6}^support for q(d) <= 0 (strict False) / < 0
    (strict True). Coordinate-by-coordinate with a completion lower-bound prune.
    Increments state['checked'] per full evaluation; raises _Budget when it would
    exceed state['budget']. Returns a witness dict {index: value} or None."""
    idx = list(support)
    L = len(idx)
    H = [[gram[idx[p]][idx[r]] for r in range(L)] for p in range(L)]
    hi = _BOX
    d = [0] * L

    def lb(k):
        # lower bound on q over completions with positions k..L-1 free in [1..hi]
        val = 0
        for p in range(k):
            val += d[p] * d[p]
            for r in range(p + 1, k):
                val += H[p][r] * d[p] * d[r]
            for r in range(k, L):
                term = H[p][r] * d[p]
                val += min(term * 1, term * hi)
        for p in range(k, L):
            val += 1
            for r in range(p + 1, L):
                h = H[p][r]
                val += min(h * 1 * 1, h * hi * hi)
        return val

    def rec(k):
        if k == L:
            if state["checked"] >= state["budget"]:
                raise _Budget()
            state["checked"] += 1
            val = sum(H[p][r] * d[p] * d[r] for p in range(L) for r in range(L)) // 2
            if (val < 0) if strict else (val <= 0):
                return {idx[p]: d[p] for p in range(L)}
            return None
        bound = lb(k)
        if (bound >= 0) if strict else (bound > 0):
            return None                       # prune: no completion can qualify
        for value in range(1, hi + 1):
            d[k] = value
            got = rec(k + 1)
            if got is not None:
                return got
        d[k] = 0
        return None

    return rec(0)


def _box_search(form, *, strict, budget):
    """Connected-support branch-and-bound over {0..6}^n. Returns
    (witness_full_tuple_or_None, value_or_None, completed_bool, checked_int).
    `completed` iff the whole box was swept within budget (a complete Ovsienko
    decision for the non-strict / weak-positivity use)."""
    graph = form.underlying_graph()
    state = {"checked": 0, "budget": budget}
    try:
        for support in _connected_supports(graph, form.n):
            got = _search_support(form.gram, support, strict, state)
            if got is not None:
                w = tuple(got.get(i, 0) for i in range(form.n))
                return w, form.evaluate(w), True, state["checked"]
    except _Budget:
        return None, None, False, state["checked"]
    return None, None, True, state["checked"]


def _restriction_isotropic_witness(form, budget):
    """Scan connected restrictions for a positive-SEMI-definite-not-definite one
    (an extended-Dynkin / Euclidean restriction). Its sign-definite radical z gives
    q(z) = 0 on the nonnegative cone -- an EXACT isotropic weak-positivity witness,
    found by linear algebra (no box sweep). Returns (witness_full, checked) or
    (None, checked). This is the fast route for large-support isotropic witnesses
    (e.g. the ~E8 subform inside T_{2,3,7})."""
    checked = 0
    for support in _connected_supports(form.underlying_graph(), form.n):
        if len(support) < 2:
            continue
        checked += 1
        if checked > budget:
            return None, checked
        idx = sorted(support)
        sub = form.restrict(idx)
        if _is_pd(sub):
            continue                          # positive definite: no isotropic vector
        z = _radical_witness(sub)             # PSD-not-PD => sign-definite radical
        if z is not None:
            w = tuple(z[idx.index(i)] if i in support else 0 for i in range(form.n))
            return w, checked
    return None, checked


def is_weakly_positive(form, *, budget=3_000_000) -> FormVerdict:
    """Exact weak-positivity decision (Ovsienko's box-6, cited complete):
    q weakly positive iff q(d) > 0 for every 0 != d in {0..6}^n. Fast certificates
    first -- positive-definite => weakly positive; an isotropic radical vector of
    the form OR of any connected restriction => a q = 0 witness => not weakly
    positive -- then a connected-support branch-and-bound sweep for the residual
    indefinite (strictly-negative) witnesses. None only on budget (never a guessed
    True)."""
    if _is_pd(form):
        return FormVerdict(True, None, None, "positive definite (=> weakly positive)", 0)
    z = _radical_witness(form)
    if z is not None:
        return FormVerdict(False, z, form.evaluate(z),
                           f"isotropic radical witness d={list(z)} q(d)=0", 0)
    w, checked = _restriction_isotropic_witness(form, budget)
    if w is not None:
        return FormVerdict(False, w, form.evaluate(w),
                           f"isotropic radical of a Euclidean restriction "
                           f"d={list(w)} q(d)=0", checked)
    w, val, completed, box_checked = _box_search(form, strict=False, budget=budget)
    checked += box_checked
    if w is not None:
        return FormVerdict(False, w, val,
                           f"witness d={list(w)} q(d)={val} (<= 0)", checked)
    if not completed:
        return FormVerdict(None, None, None,
                           f"budget_exceeded (n={form.n})", checked)
    return FormVerdict(True, None, None,
                       "weakly positive (Ovsienko box-6 exhausted, no root <= 0)", checked)


def is_weakly_nonnegative(form, *, budget=3_000_000) -> FormVerdict:
    """Weak-nonnegativity decision. PRIMARY route = the classified hypercritical
    list (_tits_lists); a False is always a FOUND exact witness (q < 0):
      (1) positive-semidefinite => True (nonnegative on the whole space, hence on
          the cone; a PSD form has NO hypercritical restriction, so this is a
          genuine completeness certificate);
      (2) a connected restriction matching a listed hypercritical form => False
          with that entry's tabulated defect (catches out-of-box witnesses, e.g.
          T_{2,3,7});
      (3) a bounded branch-and-bound witness-finder finds d >= 0 with q(d) < 0
          => False (small witnesses, e.g. the 3-Kronecker (1,1));
      (4) else, if HYPERCRITICAL_COVERAGE covers n and the box was fully swept
          => True (list complete for n);
      (5) else None -- budget_exceeded, or "hypercritical list partially
          transcribed" (honest; NEVER a guessed True)."""
    from quiverlab.invariants._tits_lists import (HYPERCRITICAL_COVERAGE,
                                                  matches_hypercritical)
    if _is_psd(form):
        return FormVerdict(True, None, None,
                           "positive semidefinite (=> weakly nonnegative; no "
                           "hypercritical restriction can exist)", 0)
    # (2) list match: the whole form, then connected restrictions
    m = matches_hypercritical(form)
    if m is not None:
        w = m["defect"]
        return FormVerdict(False, w, form.evaluate(w),
                           f"matches hypercritical {m['name']} "
                           f"(defect d={list(w)} q(d)={form.evaluate(w)})", 0)
    graph = form.underlying_graph()
    checked_restr = 0
    for support in _connected_supports(graph, form.n):
        if len(support) == form.n:
            continue
        checked_restr += 1
        if checked_restr > budget:
            break
        sub = form.restrict(sorted(support))
        m = matches_hypercritical(sub)
        if m is not None:
            idx = sorted(support)
            w = tuple(m["defect"][idx.index(i)] if i in support else 0
                      for i in range(form.n))
            return FormVerdict(False, w, form.evaluate(w),
                               f"restriction on {[form.labels[i] for i in idx]} "
                               f"matches hypercritical {m['name']}", checked_restr)
    # (3) sound witness-finder
    w, val, completed, checked = _box_search(form, strict=True, budget=budget)
    if w is not None:
        return FormVerdict(False, w, val,
                           f"witness d={list(w)} q(d)={val} (< 0)", checked)
    # (4) certified True only under recorded completeness
    if completed and form.n in HYPERCRITICAL_COVERAGE:
        return FormVerdict(True, None, None,
                           "no hypercritical restriction; list complete for "
                           f"n={form.n}", checked)
    if not completed:
        return FormVerdict(None, None, None, f"budget_exceeded (n={form.n})", checked)
    return FormVerdict(None, None, None,
                       "hypercritical list partially transcribed -- no certified "
                       f"weakly-nonnegative verdict for n={form.n} vertices", checked)


# ---------------------------------------------------------------------------
# the tame/wild certificate (P56-gated; verdict on an algebraically closed base)
# ---------------------------------------------------------------------------
def _is_alg_closed(A) -> bool:
    """The representation-type verdict gate: True exactly when the base field
    admits the tame/wild reading. Post-P61 this reads the domain's
    ``is_algebraically_closed`` flag (True only on the CC working domain; QQ,
    QQ(i) and every GF(p)/GF(p^n) are False). On a branch predating that flag it
    falls back to ``characteristic == 0`` -- a char-0 verdict is then justified by
    base change to the algebraic closure (the combinatorial Tits form and the
    representation type are field-independent in characteristic 0). A
    positive-characteristic field is refused in either regime (finite, not
    algebraically closed, and quiverlab does not model its algebraic closure)."""
    return getattr(A.domain, "is_algebraically_closed", A.domain.characteristic == 0)


@dataclass(frozen=True)
class TameWildCertificate:
    # -- form layer (field-free; present for every triangular presented algebra) --
    gram: tuple
    minimal_relation_counts: dict
    is_unit_form: bool
    weakly_positive: object            # True / False / None
    weakly_nonnegative: object
    witness: object                    # exact d >= 0 with q(d) <= 0 / < 0
    witness_value: object
    # -- P56 certificate trail --
    simply_connected: object
    strongly_simply_connected: object
    strong_certificate: object
    # -- verdict layer (algebraically-closed-base gate + P56 gate) --
    field_alg_closed: bool             # the verdict field gate (flag if present, else char 0)
    rep_type: object                   # "rep-finite" | "tame" | "wild" | None
    reason: str
    scope_note: str
    certified: object = None           # "rep_infinite" when rep_type is None but
    #                                    Bongartz certifies representation-infinite
    #                                    (not weakly positive + simply connected)


_SCOPE_NOTE = (
    "The combinatorial Tits form q_A and its weak positivity / weak nonnegativity "
    "are field-free integer arithmetic. The representation-type verdict requires an "
    "algebraically closed base field (quiverlab's CC domain; over a characteristic-0 "
    "field it is read by base change to the algebraic closure, the Tits form being "
    "field-independent in characteristic 0 -- GF(p)/GF(p^n) is refused) and is gated "
    "on the P56 certificate: rep-finite <=> weakly positive for a SIMPLY connected "
    "algebra (Bongartz 1984); tame <=> weakly nonnegative for a STRONGLY simply "
    "connected algebra (Brustle-de la Pena-Skowronski 2011), decided here by the "
    "exact positive-semidefinite (Euclidean) certificate. wild = a found d >= 0 with "
    "q_A(d) < 0. Off scope the form is still computed and the verdict is None."
)


def tame_wild_certificate(A, *, convex_budget=20000,
                          search_budget=3_000_000) -> TameWildCertificate:
    """The Tits-form representation-type certificate (Plan 62 / R19): the
    rep-finite / tame / wild trichotomy gated on the P56 strong-simple-connectivity
    certificate over an algebraically closed base field (a characteristic-0 field is
    read by base change to the algebraic closure). Loud on presentation-less input;
    raises (via as_unit_form) on a non-triangular quiver. P56's None propagates to
    a None verdict -- never a fabricated tame/wild."""
    _require_quiver(A, "tame_wild_certificate")
    form = as_unit_form(A)                    # raises unless acyclic (M2 gate)
    from quiverlab.invariants.coverings import minimal_relation_counts
    mrc = minimal_relation_counts(A)

    wp = is_weakly_positive(form, budget=search_budget)
    wnn = is_weakly_nonnegative(form, budget=search_budget)
    if wnn.witness is not None:
        witness, witness_value = wnn.witness, wnn.witness_value
    else:
        witness, witness_value = wp.witness, wp.witness_value

    sc = A.is_simply_connected(strong=True, convex_budget=convex_budget)
    scv = sc.verdict
    ssc_obj = sc.strongly
    sscv = ssc_obj.verdict if ssc_obj is not None else None

    field_ok = _is_alg_closed(A)

    rep_type = None
    certified = None
    if not field_ok:
        p = getattr(A.domain, "characteristic", 0)
        why = (f"this field has characteristic {p}" if p else
               "this characteristic-0 field is not flagged algebraically closed")
        reason = ("the representation-type verdict requires an algebraically closed "
                  f"base field (quiverlab's CC working domain); {why}, so the tame/wild "
                  "reading is out of scope. The Tits form above is field-free and "
                  "reported.")
    elif sscv is True:
        # tame/wild axis (BdlPS, strong simple connectivity)
        if wp.holds is True:
            rep_type = "rep-finite"
            reason = ("weakly positive => representation-finite (Bongartz 1984); "
                      "strongly simply connected, verdict read by base change to the "
                      "algebraic closure (the Tits form is field-independent in "
                      "characteristic 0).")
        elif wnn.holds is False:
            rep_type = "wild"
            reason = (f"a witness d={list(witness)} with q_A(d)={witness_value} < 0 was "
                      f"found => not weakly nonnegative => wild "
                      f"(Brustle-de la Pena-Skowronski 2011); strongly simply connected.")
        elif wnn.holds is True and wp.holds is False:
            rep_type = "tame"
            reason = ("weakly nonnegative but not weakly positive => tame "
                      "(Brustle-de la Pena-Skowronski 2011); strongly simply "
                      "connected. The isotropic witness records the tame direction.")
        elif wp.holds is False:
            reason = ("strongly simply connected and NOT weakly positive => "
                      "representation-infinite (Bongartz 1984); the tame/wild split is "
                      f"undetermined (weak nonnegativity = {wnn.reason}). Verdict "
                      "withheld between tame and wild (no guess).")
        else:
            reason = ("strongly simply connected, but the form verdict is "
                      f"undetermined: weak positivity = {wp.reason}; weak "
                      f"nonnegativity = {wnn.reason}. Verdict withheld (no guessed "
                      "tame/wild).")
    elif scv is True:
        # Bongartz rep-finite axis only (simple, not strong)
        if wp.holds is True:
            rep_type = "rep-finite"
            reason = ("weakly positive => representation-finite (Bongartz 1984, "
                      "simply connected). The tame/wild split is withheld: strong "
                      "simple connectivity is not certified "
                      f"({ssc_obj.reason if ssc_obj else 'not computed'}).")
        elif wp.holds is False:
            reason = ("not weakly positive => representation-INFINITE (Bongartz "
                      "1984, simply connected), but the tame/wild split needs "
                      "STRONG simple connectivity, which is not certified here "
                      f"({ssc_obj.reason if ssc_obj else 'not computed'}). Verdict "
                      "withheld (rep-infinite, tame-vs-wild undetermined).")
        else:
            reason = ("simply connected, but weak positivity is undetermined "
                      f"({wp.reason}). Verdict withheld.")
    else:
        # neither strongly nor simply connected is certified True
        if sscv is False:
            reason = (f"not strongly simply connected (P56: {ssc_obj.reason}) -- the "
                      f"Tits-form tame/wild verdict is out of scope.")
        elif scv is False:
            reason = (f"not simply connected (P56: {sc.reason}) -- the Tits-form "
                      f"representation-type verdict is out of scope.")
        else:
            reason = (f"simple connectivity undecided (P56: {sc.reason}) -- verdict "
                      f"None propagated; never upgraded to a tame/wild claim.")

    # Bongartz-certain representation-INFINITE (Plan 62 ruling 4): not weakly positive
    # on a (strongly) simply connected algebra over the verdict field certifies
    # representation-infinite, even when the tame/wild split stays undetermined. Surfaced
    # explicitly so the GUI/report can STATE the certainty; rep_type stays None because
    # tame-vs-wild is genuinely unknown here (never conflated with the undecided case).
    if (rep_type is None and field_ok and wp.holds is False
            and (scv is True or sscv is True)):
        certified = "rep_infinite"

    return TameWildCertificate(
        gram=form.gram,
        minimal_relation_counts=dict(mrc),
        is_unit_form=True,
        weakly_positive=wp.holds,
        weakly_nonnegative=wnn.holds,
        witness=witness,
        witness_value=witness_value,
        simply_connected=scv,
        strongly_simply_connected=sscv,
        strong_certificate=ssc_obj,
        field_alg_closed=field_ok,
        rep_type=rep_type,
        reason=reason,
        scope_note=_SCOPE_NOTE,
        certified=certified,
    )
