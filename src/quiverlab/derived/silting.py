"""Silting objects in K^b(proj A): verifier, single mutation, bounded exploration,
and the co-t-structure dictionary (Plan 67 / Aihara-Iyama arXiv:1009.3370, J. LMS 85
(2012) 633-668).

A silting object T satisfies Hom_{D^b}(T, T[n]) = 0 for all n > 0 (PRESILTING) and
thick(T) = K^b(proj A) (GENERATION). This is WEAKER than tilting (which also needs the
n < 0 vanishing); tilting => silting (AI Def 2.1). Presilting is DECIDED on the exact
positive window (perfect => bounded => finite window; outside it hyper-Hom is provably
the zero cochain group). Generation is three-valued: certified True on the CERTIFIED
classes -- a 2-term K0-basis object (IJY/AIR completion) incl. the regular object A (the
"tilting (AI Ex 2.2)" rung fires only on ``is_tilting_complex(...).is_tilting is True``,
itself certified-only), and a local indecomposable (AI Thm 2.26) -- else 'unknown'.
det(g_proj) = +-1 is NECESSARY (AI Thm 2.27) but NOT sufficient in general: for a WIDE
(non-2-term) rigid K0-basis object, "rigid + (#summands = rk K0) => generation" is exactly
Rickard's rank QUESTION, still OPEN (thick subcategories are not K0-classified; Krah
phantom, arXiv:2302.12502). Such an object falls through to the K0-basis-only rung and
returns 'unknown' -- the SAME verdict ``is_tilting_complex`` gives it (Plan 67 fix round:
the tilting rung and the K0-only rung are reconciled; neither certifies from bare K0). The
g-matrix is computed in the
PROJECTIVE basis K0(K^b proj A) = (+)_v Z[P_v] via ``derived.tilting.g_proj`` (Plan 67
Task 0), NOT the composition-factor basis -- ``_chi`` gives det(Cartan . g) and is wrong
on non-unimodular Cartan (self-injective/symmetric; Example 2.47).

Single mutation (AI Def 2.30/2.34) is one approximation triangle: at an indecomposable
summand X, the cone of the minimal left add(T/X)-approximation (left mutation mu^+) or the
cocone of the right one (right mutation mu^-). Every mutant is RE-VERIFIED PRESILTING (the
decidable positive-window half, NOT full generation), shares exactly n-1 summand PROFILES
with its input, and mu^- o mu^+ = id (AI Thm 2.31 / Prop 2.33). The
bounded-radius exploration truncates loudly (the silting quiver can be infinite -- kA2
already is; transitivity is proven only for local/hereditary/canonical, AI Thm 1.2, so
there is NO general BFS/enumeration claim). Float-free; exact linear algebra only."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.modules.complexes import hyper_hom_dims
from quiverlab.derived.tilting import (is_tilting_complex, _direct_sum_complex, _span,
                                       _is_two_term, g_proj)


# --------------------------------------------------------------------------- #
# the silting-object verifier
# --------------------------------------------------------------------------- #
@dataclass
class SiltingReport:
    is_silting: object          # True | False | "unknown"  (three-valued; see ruling)
    is_presilting: bool         # Hom_{D^b}(T, T[n]) == 0 for all n > 0 in the window
    k0_basis: bool              # g_proj square (#summands == #simples) AND det(g_proj) +-1
    window: tuple               # (1, n_max): the EXACT positive rigidity-check range
    g_matrix: list              # rows = summand K0 classes in the PROJECTIVE basis
                                # (+)_v Z[P_v] via g_proj -- NOT the composition-factor
                                # basis (see Task 0 / the sufficiency ruling)
    det: int                    # det(g_proj) -- unimodularity in the projective basis
    generation_certified_by: str


def _positive_window(summands):
    """The EXACT positive rigidity window (1, n_max). For perfect summands spanning
    degrees [lo_i, hi_i], Hom^n(T_i, T_j) is the zero cochain group for
    n > max_i hi_i - min_j lo_j, so scanning [1, n_max] fully DECIDES presilting."""
    spans = [_span(T) for T in summands]
    return 1, max(hi for _, hi in spans) - min(lo for lo, _ in spans)


def is_silting_object(summands):
    """Decide whether ``summands`` (a list of INDECOMPOSABLE perfect complexes) assemble a
    silting object of K^b(proj A). Presilting is DECIDED on the exact positive window
    (perfect => bounded => finite; outside it hyper-Hom is provably 0). Generation is
    three-valued per the sufficiency ladder (True on tilting/2-term/local, 'unknown' on
    K0-basis-only, False otherwise). Returns a :class:`SiltingReport`; loud if any summand
    is not a certified perfect complex."""
    import sympy as sp
    if not summands:
        raise QuiverlabError("is_silting_object: need at least one summand")
    for T in summands:
        if not T.is_perfect():
            raise QuiverlabError("is_silting_object: every summand must be a "
                                 "certified perfect complex")
    A = summands[0].algebra
    verts = list(A.quiver.vertices)
    n_lo, n_max = _positive_window(summands)
    Tsum = _direct_sum_complex(summands)
    hh = hyper_hom_dims(Tsum, Tsum, n_lo, n_max) if n_max >= n_lo else {}
    presilting = all(hh.get(n, 0) == 0 for n in range(n_lo, n_max + 1))
    g = [g_proj(T, verts) for T in summands]       # K0 in the PROJECTIVE basis (Task 0)
    det = int(sp.Matrix(g).det()) if len(g) == len(verts) else 0
    k0_basis = (len(g) == len(verts)) and det in (1, -1)
    # sufficiency ladder (see the module docstring / the plan sufficiency ruling)
    is_silting, why = False, ""
    if is_tilting_complex(summands).is_tilting is True:   # `is True`: the tilting verifier
        # is now itself three-valued -- is_tilting is True ONLY when generation is CERTIFIED
        # (2-term/regular, IJY completion), "unknown" for a non-2-term rigid K0-basis object
        # (Rickard's rank question, OPEN). We fire this rung on hard True ONLY; an "unknown"
        # tilting verdict deliberately falls through to the K0-basis-only rung below, where
        # it lands on the SAME honest "unknown" -- the two verifiers agree exactly in the
        # Rickard-open regime (Plan 67 fix round: rung 1 and the K0-only rung reconciled).
        is_silting, why = True, "tilting (AI Ex 2.2)"
    elif presilting and _is_two_term(summands) and k0_basis:
        # IJY: a 2-term PRESILTING object whose #(distinct indec) summands == #simples is
        # silting. k0_basis = (#summands == #simples) AND det(g_proj) in {+-1}: the det
        # condition is automatic for a genuine (basic) 2-term silting object AND it
        # correctly REJECTS degenerate non-basic inputs (a repeated summand [P, P] has
        # #summands == #simples but det(g_proj) = 0 -> falls through to False, never a
        # false True or a crash). Do NOT drop det from the condition.
        is_silting, why = True, "2-term (AIR/IJY completion)"
    elif presilting and len(verts) == 1 and len(summands) == 1:
        is_silting, why = True, "local (AI Thm 2.26)"
    elif presilting and k0_basis:
        is_silting = "unknown"
        why = "K0-basis only (necessary; not sufficient -- Krah phantom)"
    return SiltingReport(is_silting=is_silting, is_presilting=presilting,
                         k0_basis=k0_basis, window=(n_lo, n_max), g_matrix=g,
                         det=det, generation_certified_by=why)


# --------------------------------------------------------------------------- #
# the co-t-structure dictionary (AI Prop 2.23(b) / Jorgensen)
# --------------------------------------------------------------------------- #
_CT_REFS = ["aihara_iyama_silting", "jorgensen_cotstructures"]


def co_t_structure_of(summands):
    """The bounded co-t-structure with coheart add(T) (AI Prop 2.23(b) / Jorgensen). A
    DOCUMENTATION record -- NOT a computed subcategory (the aisles are infinite): it ships
    the coheart's finite per-summand dim-vectors plus descriptor strings + references.
    Refuses loudly if ``summands`` is not a (certified or candidate) silting object
    (``is_silting`` in {True, 'unknown'})."""
    rep = is_silting_object(summands)
    if rep.is_silting not in (True, "unknown"):
        raise QuiverlabError(
            "co_t_structure_of: input is not a (certified or candidate) silting object; "
            "no co-t-structure to report", hint=str(rep.generation_certified_by))
    coheart = [{"dimvecs": {n: dict(T.term(n).dimension_vector()) for n in T.degrees()}}
               for T in summands]
    return {"coheart": coheart, "bounded": True,
            "aisle": "T_M^{<=0} = { X : Hom_{D^b}(T, X[>0]) = 0 }  (AI Def 2.12)",
            "coaisle": "^perp(T_M^{<=0})  (AI Prop 2.23(b) torsion pair)",
            "coheart_is_addT": True, "references": list(_CT_REFS)}


# --------------------------------------------------------------------------- #
# single silting mutation via one approximation triangle (AI Def 2.30/2.34)
#
# The mutation MATH (minimal add(T/X)-approximation in K^b(proj), its cone / the dual
# cocone, and the minimal-complex reduction that the involution + exploration-dedup
# oracles REQUIRE) is exactly P45's retract-prune already lifted to K^b(proj) over the
# projective-vertex complex representation (``tautilting._twoterm``): ``min_left_approx`` /
# ``min_right_approx`` (P44 retract-pruning, lifted), ``left_mutation_summand`` /
# ``right_mutation_summand`` = cone / cocone[-1] (AI Def 2.34), ``reduce_complex`` =
# delooping/Gaussian elimination of iso blocks (the minimal representative). We reuse that
# TESTED engine rather than re-deriving it on ``ChainComplex`` (the plan's Task-2 sketch
# omitted the minimizer, which the involution oracle needs); the public surface stays
# ``ChainComplex`` and the mutant is re-verified PRESILTING (the decidable positive-window
# half) + shares-n-1 by the verifier above -- full generation is not re-derived here.
# The two representations differ only by the degree convention (PComplex is cohomological
# d^i: C^i -> C^{i+1}; ChainComplex is homological d_n: C_n -> C_{n-1}), so n = -i is the
# exact bijection and the differential matrix is byte-identical (same shape, same
# concatenated builders.projective basis).
# --------------------------------------------------------------------------- #
def _vertex_lists(cx):
    """Per (homological) degree, the ordered projective-summand vertex list of ``C_n``'s
    basis blocks. Uses the ``_proj_vertices`` provenance when present (every engine-built
    complex -- from_projective_resolution / cone-mutants / two_term_silting_from_presentation
    -- carries it); else recovers it for a term that is a SINGLE indecomposable projective
    (``top`` is 1-dim => one vertex), which covers the ``ChainComplex.stalk(A.projective(v))``
    seeds. Raises loudly on a multi-block term with no provenance."""
    prov = getattr(cx, "_proj_vertices", None)
    out = {}
    for n in cx.degrees():
        if prov is not None and n in prov:
            out[n] = list(prov[n])
            continue
        tv = cx.term(n).top().dimension_vector()
        blocks = [(v, m) for v, m in tv.items() if m]
        if sum(m for _, m in blocks) == 1:
            out[n] = [blocks[0][0]]
        else:
            raise QuiverlabError(
                "silting_mutate: a summand term carries no projective-block provenance "
                "and is not a single indecomposable projective; supply a complex built by "
                "the silting engine (stalks of A.projective(v) or prior mutants)")
    return out


def _cx_to_pc(cx):
    """A perfect ``ChainComplex`` (homological) -> the ``tautilting._twoterm.PComplex``
    (cohomological) it equals under ``i = -n``: same vertex-lists, same differential
    matrices."""
    from quiverlab.tautilting import _twoterm as tt
    vl = _vertex_lists(cx)
    terms = {-n: list(vl[n]) for n in vl}
    diffs = {-n: mat for n, mat in cx._dmats.items() if mat and mat[0]}
    return tt.PComplex(cx.algebra, terms, diffs)


def _pc_to_cx(pc):
    """A ``PComplex`` (cohomological) -> the perfect ``ChainComplex`` (homological) it
    equals under ``n = -i``. Terms are ``(+)_v builders.projective(A, v)`` in the vertex
    order (the same basis PComplex uses), ``check=True`` re-certifies d.d = 0, and
    ``_proj_vertices`` provenance is carried so a subsequent mutation round-trips."""
    from quiverlab.tautilting import _twoterm as tt
    from quiverlab.modules.complexes import ChainComplex
    A = pc.algebra
    terms, dmats, prov = {}, {}, {}
    for i in pc.degrees():
        M = tt._proj_sum(A, pc.terms[i])
        if M.dim:
            terms[-i] = M
            prov[-i] = list(pc.terms[i])
    for i, d in pc.diffs.items():
        n = -i
        if d and d[0] and (n in terms) and ((n - 1) in terms):
            dmats[n] = d
    if not terms:
        raise QuiverlabError("silting mutation produced the zero complex (degenerate) -- "
                             "no silting summand to report")
    cx = ChainComplex(terms, dmats, check=True)
    cx._perfect = True
    cx._proj_vertices = prov
    return cx


def _summand_key(cx):
    """A degree-sensitive per-summand fingerprint (per-degree dim-vectors). This is a
    NECESSARY, not canonical, identifier of a summand: two non-isomorphic perfect complexes
    CAN share a per-degree dim-vector profile (a canonical iso-class key for complexes is
    not cheap -- cf. ``modules/hom.py::identify_standard`` for modules). ``_assert_neighbour``
    handles that ambiguity explicitly (loudly) rather than silently mislabelling it."""
    return tuple(sorted((n, tuple(sorted(cx.term(n).dimension_vector().items())))
                        for n in cx.degrees()))


def _assert_neighbour(T, mutant, n):
    """Self-certificate for a single mutation (never trusts the construction). Two DECIDABLE
    checks: the mutant has ``n`` summands and RE-VERIFIES **PRESILTING** (AI Thm 2.31 -- a
    mutation of a silting object is silting; presilting is the positive-window half we
    actually decide here, NOT full generation), plus a fingerprint check that it shares
    EXACTLY ``n-1`` summand PROFILES with the input (a Hasse neighbour differing from it).

    The fingerprint (``_summand_key``, per-degree dim-vectors) is not canonical, so the
    shares-count is split into two honest branches:
    - ``shared < n-1`` -- the mutation DEGENERATED: too few surviving profiles, so the
      approximation / cone / shift is wrong (a genuine failure).
    - ``shared > n-1`` (== n) -- FINGERPRINT AMBIGUITY: the new cone ``N_X`` carries the
      SAME per-degree dim-vector profile as the mutated summand ``X`` (or the mutant equals
      the input up to reordering), so the dim-vector fingerprint cannot certify the mutant
      DIFFERS from the input. This is NOT proof the mutation is wrong (the mutant is still
      presilting with ``n`` summands) -- it is a limitation of the coarse fingerprint. We
      refuse LOUDLY and honestly rather than mislabel a possibly-valid mutation as a
      non-neighbour; a canonical iso-class key for complexes would resolve it but is not
      cheap. (Never observed on the shipped batteries -- a defensive guard.)"""
    if len(mutant) != n:
        raise QuiverlabError(
            f"silting_mutate: mutant has {len(mutant)} summands, expected {n}")
    if not is_silting_object(mutant).is_presilting:
        raise QuiverlabError(
            "silting_mutate: the mutant is not presilting (AI Thm 2.31 says a mutation of "
            "a silting object is silting) -- the approximation / cone / shift is wrong")
    tkeys = [_summand_key(c) for c in T]
    used = [False] * len(tkeys)
    shared = 0
    for c in mutant:
        ck = _summand_key(c)
        for k, tk in enumerate(tkeys):
            if not used[k] and tk == ck:
                used[k] = True
                shared += 1
                break
    if shared < n - 1:
        raise QuiverlabError(
            f"silting_mutate: mutant shares only {shared} summand profiles with the input "
            f"(expected n-1 = {n - 1}) -- the mutation DEGENERATED (the approximation / "
            "cone / shift is wrong)")
    if shared > n - 1:
        raise QuiverlabError(
            f"silting_mutate: the mutant's per-degree dim-vector profiles coincide with the "
            f"input's ({shared} == n = {n}), so the fingerprint cannot certify the mutant "
            "DIFFERS from the input -- a FINGERPRINT AMBIGUITY (the new cone shares the "
            "mutated summand's profile), NOT necessarily a degeneration; report the algebra "
            "+ summand index")


def silting_mutate(summands, i, direction="left"):
    """Irreducible mutation at summand index ``i`` (AI Def 2.34), ``direction`` in
    ``{"left", "right"}``. SELF-CERTIFIED: the input must be silting or presilting +
    K0-basis (loud otherwise); the result RE-VERIFIES **PRESILTING** (AI Thm 2.31 -- the
    decidable positive-window half, NOT full generation) and shares exactly ``n-1`` summand
    PROFILES with the input, differing from it. Returns the new list ``[rest..., N_X]`` of
    minimal perfect complexes. Left mutation ``mu^+`` is the cone of the minimal left
    ``add(T/X)``-approximation; right mutation ``mu^-`` is the cocone of the right one;
    ``mu^- o mu^+ = id`` (AI Prop 2.33)."""
    from quiverlab.tautilting import _twoterm as tt
    rep = is_silting_object(summands)
    if rep.is_silting not in (True, "unknown"):
        raise QuiverlabError("silting_mutate: input is not a (candidate) silting object",
                             hint=rep.generation_certified_by)
    n = len(summands)
    if not (0 <= i < n):
        raise QuiverlabError(f"silting_mutate: summand index {i} out of range 0..{n - 1}")
    pcs = [_cx_to_pc(T) for T in summands]
    Xk = pcs[i]
    U = pcs[:i] + pcs[i + 1:]
    if direction == "left":
        new_pc = tt.left_mutation_summand(Xk, U)              # cone(min left approx)
    elif direction == "right":
        new_pc = tt.right_mutation_summand(Xk, U)             # cocone(min right approx)
    else:
        raise QuiverlabError("silting_mutate: direction must be 'left' or 'right'")
    N = _pc_to_cx(tt.reduce_complex(new_pc))                  # minimal representative
    mutant = summands[:i] + summands[i + 1:] + [N]
    _assert_neighbour(summands, mutant, n)                    # self-cert (AI Thm 2.31)
    return mutant


def silting_neighbors(summands, direction="left"):
    """The (up to) ``n`` irreducible mutations -- the single-step neighbourhood, as a list
    of ``(i, mutant)``. A summand whose approximation degenerates (mutation not defined /
    refuses there) is reported with ``mutant=None`` -- never a silent skip."""
    out = []
    for i in range(len(summands)):
        try:
            out.append((i, silting_mutate(summands, i, direction=direction)))
        except QuiverlabError:
            out.append((i, None))
    return out


# --------------------------------------------------------------------------- #
# bounded-radius exploration with LOUD truncation (AI Thm 1.2 -- no general BFS)
#
# The silting quiver can be INFINITE (kA2 already is, AI Example 2.45) and
# mutation-transitivity is proven only for local / hereditary / canonical (AI Thm 1.2,
# and it FAILS for a symmetric algebra [AGI]). So there is no general BFS/enumeration
# claim: this is a bounded-radius walk from a silting object (the regular A, or a supplied
# one) through irreducible left + right mutations, deduped by a shift-insensitive silting
# key, stopping with a loud ``status`` -- ``"complete"`` ONLY for the local class (Thm
# 2.26: silt = {A[i]}, radius 0 already closes; ``finite_class in {"local", None}``, there
# is deliberately NO ``"two_term"`` -- a general mutation walk cannot be restricted to the
# 2-term slice, so 2-term finiteness is P45's, cross-checked directly, not claimed here).
# --------------------------------------------------------------------------- #
@dataclass
class SiltingExploration:
    vertices: list        # {"summands", "key", "is_initial", "silting"}
    arrows: dict          # {(i, j): {"direction": "left"|"right", "summand": <index>}}
    status: str           # "complete" | "radius" | "budget"
    radius: int           # the radius actually reached
    finite_class: object  # "local" | None (local is the ONLY certified-complete class)


def _silting_key(summands):
    """A SHIFT-INSENSITIVE fingerprint of a silting object: normalise the common shift
    (anchor the whole object at minimal degree 0 -- AI identifies T with T[i]) and take the
    frozenset of the per-summand (shift-normalised) per-degree dim-vector tuples. Dedups
    ``A`` and ``A[3]``."""
    degs = [d for T in summands for d in T.degrees()]
    base = min(degs) if degs else 0
    parts = []
    for T in summands:
        parts.append(tuple(sorted(
            (n - base, tuple(sorted(T.term(n).dimension_vector().items())))
            for n in T.degrees())))
    return frozenset(parts)


def _explore_record(summands, is_initial):
    return {"summands": summands, "key": _silting_key(summands),
            "is_initial": is_initial, "silting": is_silting_object(summands).is_silting}


def bounded_silting_exploration(start, radius=3, budget=256):
    """Bounded-radius BFS of the silting quiver from ``start`` (an :class:`Algebra` ->
    begin at the regular silting ``A``, or a silting summand list) through irreducible
    left+right mutations, deduped by :func:`_silting_key`. LOUD ``status``: ``"complete"``
    only for a proven-finite class (LOCAL only -- Thm 2.26), ``"radius"`` when the radius
    bound is reached (or the ball closes with no finiteness theorem -- a closed ball is not
    a completeness proof), ``"budget"`` when the vertex cap trips first. NEVER claims
    completeness outside the local class."""
    from quiverlab.modules.complexes import ChainComplex
    if hasattr(start, "quiver"):                      # an Algebra: the regular object A
        A = start
        init = [ChainComplex.stalk(A.projective(v), 0) for v in A.quiver.vertices]
    else:
        init = list(start)
        A = init[0].algebra
    verts = list(A.quiver.vertices)
    # LOCAL (one simple): silt = {A[i]} (Thm 2.26) -- radius-0 already closes -> complete.
    if len(verts) == 1:
        return SiltingExploration(vertices=[_explore_record(init, True)], arrows={},
                                  status="complete", radius=0, finite_class="local")
    records = [_explore_record(init, True)]
    index = {records[0]["key"]: 0}
    depth = {0: 0}
    arrows = {}
    frontier = [0]
    while frontier:
        i = frontier.pop(0)
        if depth[i] >= radius:                        # do not expand beyond the radius
            continue
        for direction in ("left", "right"):
            for (k, mut) in silting_neighbors(records[i]["summands"], direction):
                if mut is None:
                    continue
                key = _silting_key(mut)
                j = index.get(key)
                if j is None:
                    if len(records) >= budget:        # LOUD cap before the (budget+1)-th
                        return SiltingExploration(vertices=records, arrows=arrows,
                                                  status="budget", radius=radius,
                                                  finite_class=None)
                    j = len(records)
                    index[key] = j
                    depth[j] = depth[i] + 1
                    records.append(_explore_record(mut, False))
                    frontier.append(j)
                if i != j:
                    e = (min(i, j), max(i, j))
                    if e not in arrows:
                        arrows[e] = {"direction": direction, "summand": k}
    # frontier emptied. A closed ball is NOT a proof of completeness without a finiteness
    # theorem (only the local class, handled above, upgrades to "complete"): status stays
    # "radius" for every non-local algebra (honest -- AI Thm 1.2).
    return SiltingExploration(vertices=records, arrows=arrows, status="radius",
                              radius=radius, finite_class=None)
