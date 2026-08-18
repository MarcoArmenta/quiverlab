"""Plan 78 (R13) -- the char-0 formal-deformation surface for A = kQ/I, built on
Hochschild-cohomology deformation theory and the L-infinity machinery of
Redondo-Rossi Bertone / Mueller-Redondo-Rossi Bertone-Suarez.

The computational backbone is Chouhy-Solotar over the algebra's exact Domain (QQ
included): the shipped native CS differential (l_1) and the native CS Gerstenhaber
bracket (l_2, Plan 51). By Gerstenhaber's theorem C(A) is a genuine differential
graded Lie algebra, so the deformation functor -- infinitesimal HH^2, the primary
obstruction [alpha,alpha] in HH^3, the order-by-order Maurer-Cartan equation, the
presented deformed algebra A_alpha and its Ext-algebra -- needs only delta + l_2 and
never a higher bracket (section 4 of the plan). The B(A)[1] L-infinity structure is
a small-model companion whose exact v1 instances are the rad^2=0 dg-Lie collapse
certificate and a feasibility-gated char-0 Bardzell l_3.

SCOPE (loud typed refusals):
  * the raw HH^2 / [alpha,alpha] block computes over ANY exact field;
  * the deformation / L-infinity INTERPRETATION (maurer_cartan, is_l_infinity_nilpotent,
    dg_lie_certificate, l3_bracket, the A_alpha narrative) is CHARACTERISTIC 0 (RRB/MRRS
    are char-0 theorems; the L-infinity 1/n! and the DGLA 1/2 need char 0 / char != 2).

Bracket / obstruction CONSTANTS are basis-dependent and carry basis="cs/<domain>"
provenance; cross-engine comparisons use only basis-independent data (dims, verdicts,
obstruction rank / class nonvanishing) -- the Plan-35 rule.

Exact only; no floats in src/. No new numeric kernels -- a thin exact layer over
resolutions_cs (bracket/differential/HH), Quiver.algebra, and ext_algebra.
"""
from dataclasses import dataclass, field

from quiverlab.errors import QuiverlabError

# The plan's OWN dim backstop (NOT P70's DEFAULT_MAXDIM=48): the obstruction cost law
# is HH^2/HH^3-richness x resolution size, not A.dim, so a dim cap is only a coarse DoS
# guard bounding the cheap HH pre-probe (kZ16/J^2 HH ~ 2.4s at dim 32). H1 re-benchmark.
DEFORM_MAXDIM = 32

# every deformation block cites the five deformation papers (RRB/MRRS/RRRV/RRR/Chouhy);
# the obstruction bracket additionally rests on Gerstenhaber + the CS/NW bracket keys.
_DEFORM_REFERENCES = [
    "rrb_linfty_bardzell", "mrrs_mc_gentle", "rrrv_morita_deform",
    "rrr_ext_deform", "chouhy_degeneration",
]
_BRACKET_REFERENCES = ["gerstenhaber1963", "chouhy_solotar", "oke_koszul"]


# ---------------------------------------------------------------------------
# frozen data reports (like Plan-35's HHProducts; __bool__ is NOT defined)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Obstruction:
    hh2_dim: int
    hh3_dim: int
    basis: str                       # "cs/<domain>"
    unobstructed: bool               # the map alpha |-> [alpha,alpha] vanishes on all HH^2
    witness: str | None              # a basis-DEPENDENT display of an obstructed direction
    obstruction_constants: tuple     # the symmetric (2,2) self-bracket HH^3-classes, EXACT
                                     #   strings (basis-dependent; never compared across bases)
    characteristic: int
    char0_note: str | None
    references: tuple = ()


@dataclass(frozen=True)
class MCReport:
    characteristic: int
    hh2_dim: int
    hh3_dim: int
    basis: str
    nilpotent_regime: bool | None    # MRRS Thm 5.4 verdict (None if hypotheses not met)
    unobstructed: bool | None        # None over char p (field-general HH only)
    obstruction_witness: str | None
    mc_description: str
    mc_order_certified: int          # highest order the DGLA MC was solved + certified (>=1)
    z2_dim: int | None               # dim Z^2 (the MC set in the nilpotent regime)
    order_requested: int
    window_note: str | None
    char0_note: str | None
    note: str
    references: tuple = ()


@dataclass(frozen=True)
class L3Report:
    status: str                      # "dg_lie" | "computed" | "deferred" | "n/a"
    monomial: bool
    rad_square_zero: bool | None
    dg_lie: bool | None
    note: str
    references: tuple = ()


@dataclass(frozen=True)
class Deformations:
    characteristic: int
    hh2_dim: int
    hh3_dim: int
    basis: str
    unobstructed: bool | None
    obstruction_witness: str | None
    nilpotent_regime: bool | None
    mc_description: str
    mc_order_certified: int
    dg_lie: bool | None
    l3_status: str
    a_alpha: dict | None
    ext_algebra_summary: str | None
    base_change_note: str | None
    char0_note: str | None
    window_note: str | None
    status: str                      # "complete" | "budget" | "unsupported"
    note: str
    references: tuple = ()


# ---------------------------------------------------------------------------
# small exact helpers
# ---------------------------------------------------------------------------
def _is_zero_vec(v, dom):
    return all(dom.is_zero(x) for x in v)


def _rad_square_zero(A):
    """rad^2 = 0 iff A has no basis path of length >= 2, i.e. dim A = |Q_0| + |Q_1|
    (A = k^{Q_0} (+) span(arrows)). Cheap, char-independent, structural. Requires a
    quiver presentation; None when there is none."""
    q = A.quiver
    if q is None:
        return None
    return A.dim == len(list(q.vertices)) + len(q.arrows)


def _char0_note(A):
    if A.domain.characteristic == 0:
        return None
    return ("characteristic %d: the raw HH^2 / [alpha,alpha] data is field-general, but the "
            "formal-deformation / L-infinity interpretation (Maurer-Cartan, nilpotent regime, "
            "dg-Lie, presented A_alpha) is a characteristic-0 theory (RRB / MRRS; the 1/n! and "
            "1/2 need char 0 / char != 2)." % A.domain.characteristic)


def _require_presented(A, what):
    if A.quiver is None or A.relations is None:
        raise QuiverlabError(
            "%s needs a quiver-and-relations presentation (the Chouhy-Solotar engine "
            "cannot run on a presentation-less structure-constant algebra)" % what,
            hint="build A via Quiver.algebra(relations=[...])")


def _require_budget(A, budget, what):
    if A.dim > budget:
        raise QuiverlabError(
            "%s: A.dim = %d exceeds the deformation budget %d -- the CS obstruction bracket "
            "cost tracks HH^2/HH^3 richness x resolution size (H1 re-benchmark); this dim cap "
            "is a coarse DoS backstop" % (what, A.dim, budget),
            hint="raise the budget explicitly if you know the algebra is HH-poor")


# ---------------------------------------------------------------------------
# the shared (2,2) obstruction computation (targeted, cached liftings)
# ---------------------------------------------------------------------------
def _obstruction_data(A, *, max_cells=4_000_000):
    """The targeted (2,2) Gerstenhaber self-bracket HH^2 x HH^2 -> HH^3 on the CS HH
    basis over A.domain. Returns the raw data (Domain elements, not strings) so both
    obstruction_map and the self-cert tests read the SAME numbers:

        {"res", "dom", "coh2", "coh3", "img3", "B", "hh2_dim", "hh3_dim"}

    B[i][j] (== B[j][i]) is the HH^3-class coordinate vector of [alpha_i, alpha_j]
    (symmetric: degree-2 antisymmetry gives sign +1). The homotopy lifting of each HH^2
    class is built once and reused across its whole row/column (the Plan-51 cache; the
    canonical solve keeps it byte-identical)."""
    from quiverlab.resolutions_cs.build import reduction_system_of
    from quiverlab.resolutions_cs.homology import (
        cs_hh_basis, _require_admissible, _columns)
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    from quiverlab.resolutions_cs.bracket import native_bracket
    from quiverlab.resolutions_cs.homotopy_lifting import homotopy_lifting
    from quiverlab.hochschild.products import _class_coords

    rs = reduction_system_of(A)
    _require_admissible(rs)
    dom = A.domain
    res = ChouhySolotarResolution(A, rs, max_degree=4, max_cells=max_cells)
    coh2 = cs_hh_basis(A, 2, "coh", max_cells=max_cells)
    coh3 = cs_hh_basis(A, 3, "coh", max_cells=max_cells)
    img3 = _columns(res.matrix(2, "coh"))
    n = len(coh2)
    lifts = [homotopy_lifting(res, coh2[i], 2) for i in range(n)]
    B = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            vec = native_bracket(res, coh2[i], 2, coh2[j], 2,
                                 psi_f=lifts[i], psi_g=lifts[j])
            coords = _class_coords(vec, coh3, img3, dom, stringify=False)
            B[i][j] = coords
            B[j][i] = coords
    return {"res": res, "dom": dom, "coh2": coh2, "coh3": coh3, "img3": img3,
            "B": B, "hh2_dim": n, "hh3_dim": len(coh3)}


def self_bracket_class(data, coords):
    """The HH^3 class of [alpha, alpha] for alpha = sum coords_i * alpha_i, as a
    coordinate vector over the recorded HH^3 basis: sum_{i,j} coords_i coords_j ``B[i][j]``.
    Basis-dependent; the LOAD-BEARING fact is only whether it is zero (a lift to second
    order exists iff this class vanishes -- the DGLA recursion delta mu_2 = -1/2[a,a])."""
    dom, B, n = data["dom"], data["B"], data["hh2_dim"]
    out = [dom.zero()] * data["hh3_dim"]
    for i in range(n):
        ci = dom.coerce(coords[i])
        if dom.is_zero(ci):
            continue
        for j in range(n):
            cj = dom.coerce(coords[j])
            if dom.is_zero(cj):
                continue
            cij = dom.mul(ci, cj)
            for k in range(data["hh3_dim"]):
                out[k] = dom.add(out[k], dom.mul(cij, B[i][j][k]))
    return out


def _verdict(data):
    """(unobstructed, witness, witness_coords): unobstructed iff every pairwise bracket
    class B[i][j] vanishes (the map alpha|->[alpha,alpha] is 0 on all of HH^2, by
    polarization in char != 2). The witness is a basis-DEPENDENT display -- a single
    basis cocycle if one self-obstructs, else a two-term combination a_i + a_j (whose
    self-bracket is 2 B[i][j] != 0). The whole quadratic form is scanned; we NEVER rely
    on the diagonal alone (a basis diagonal may or may not vanish -- basis-dependent)."""
    dom, B, n = data["dom"], data["B"], data["hh2_dim"]
    unobstructed = all(_is_zero_vec(B[i][j], dom) for i in range(n) for j in range(n))
    if unobstructed:
        return True, None, None
    e = lambda i: [dom.one() if k == i else dom.zero() for k in range(n)]
    for i in range(n):                                  # a basis cocycle self-obstructs?
        if not _is_zero_vec(B[i][i], dom):
            return False, "a_%d" % i, e(i)
    for i in range(n):                                  # else an off-diagonal combination
        for j in range(i + 1, n):
            if not _is_zero_vec(B[i][j], dom):
                coords = [dom.one() if k in (i, j) else dom.zero() for k in range(n)]
                if not _is_zero_vec(self_bracket_class(data, coords), dom):
                    return False, "a_%d + a_%d" % (i, j), coords
    # char-2 pathology (2 B[i][j] = 0): fall back to a full nonzero-B display, honest
    for i in range(n):
        for j in range(i, n):
            if not _is_zero_vec(B[i][j], dom):
                return False, "a_%d, a_%d (bracket nonzero; char-2 form)" % (i, j), None
    return False, None, None


# ---------------------------------------------------------------------------
# public: infinitesimal deformations + the obstruction map
# ---------------------------------------------------------------------------
def infinitesimal_deformations(A, *, engine="auto", max_cells=4_000_000):
    """The tangent space to the formal-deformation functor: dim HH^2(A) (the
    infinitesimal deformations), over any exact field, from the shipped HH engine.
    Returns a dict {hh2_dim, basis, characteristic, char0_note, references}."""
    _require_presented(A, "infinitesimal_deformations")
    _require_budget(A, DEFORM_MAXDIM, "infinitesimal_deformations")
    table = A.hochschild_cohomology(2, engine=engine, max_cells=max_cells)
    return {
        "hh2_dim": table[2],
        "basis": "cs/%s" % A.domain.name,
        "characteristic": A.domain.characteristic,
        "char0_note": _char0_note(A),
        "references": list(_DEFORM_REFERENCES),
    }


def obstruction_map(A, *, engine="auto", max_cells=4_000_000):
    """The primary obstruction alpha |-> [alpha,alpha] in HH^3 -- the targeted (2,2) CS
    Gerstenhaber self-bracket. `unobstructed` is True iff the map vanishes on ALL of HH^2
    (the whole quadratic form Q(c) = sum c_i c_j [alpha_i,alpha_j] is scanned, never the
    diagonal alone). `witness` is a basis-DEPENDENT display of an obstructed direction.
    Field-general: over char p the block is still returned (with a char0_note) -- the
    DEFORMATION reading is char-0 gated, the raw bracket is not (`engine` is accepted for
    API symmetry; the obstruction always uses the CS-native bracket)."""
    _require_presented(A, "obstruction_map")
    _require_budget(A, DEFORM_MAXDIM, "obstruction_map")
    dom = A.domain
    data = _obstruction_data(A, max_cells=max_cells)
    unobstructed, witness, _ = _verdict(data)
    B, n = data["B"], data["hh2_dim"]
    consts = tuple(
        tuple(tuple(str(x) for x in B[i][j]) for j in range(n)) for i in range(n))
    return Obstruction(
        hh2_dim=data["hh2_dim"], hh3_dim=data["hh3_dim"],
        basis="cs/%s" % dom.name, unobstructed=unobstructed, witness=witness,
        obstruction_constants=consts, characteristic=dom.characteristic,
        char0_note=_char0_note(A),
        references=tuple(_DEFORM_REFERENCES + _BRACKET_REFERENCES))


# ---------------------------------------------------------------------------
# the nilpotent-regime gate (MRRS Thm 5.4) + the Maurer-Cartan report (char 0)
# ---------------------------------------------------------------------------
def _has_parallel_arrows(q):
    st = list(q.arrows.values())
    return len(st) != len(set(st))


def _is_acyclic(q):
    """True iff the quiver has no oriented cycle (a loop is a length-1 cycle)."""
    import collections
    adj = collections.defaultdict(list)
    for (s, t) in q.arrows.values():
        adj[s].append(t)
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {v: WHITE for v in q.vertices}
    cyclic = [False]

    def dfs(u):
        color[u] = GRAY
        for w in adj[u]:
            if color[w] == GRAY:
                cyclic[0] = True
            elif color[w] == WHITE:
                dfs(w)
        color[u] = BLACK
    for v in q.vertices:
        if color[v] == WHITE:
            dfs(v)
    return not cyclic[0]


def is_l_infinity_nilpotent(A):
    """The MRRS Thm 5.4 gate: True iff the L-infinity structure on B(A)[1] is certified
    nilpotent (so MC = Z^2). Char-0 conservative reading of R13's stated hypotheses --
    gentle, no parallel arrows, no oriented cycles. Returns None (NOT certified; the full
    MC equation applies) outside these conditions: char != 0, non-gentle, parallel arrows,
    or an oriented cycle. (`# PIN`: the exact Thm-5.4 quiver conditions were not resolvable
    against the paper's section 5 at authoring; R13's conditions are used conservatively --
    a True is always safe, a None never over-claims.)"""
    if A.domain.characteristic != 0:
        return None
    if A.quiver is None:
        return None
    try:
        gentle = A.is_gentle()
    except QuiverlabError:
        return None
    if not gentle:
        return None
    q = A.quiver
    if _has_parallel_arrows(q) or not _is_acyclic(q):
        return None
    return True


def _z2_dim(A, *, max_cells=4_000_000):
    """dim Z^2 = dim ker(delta^2 : C^2 -> C^3) on the CS cochain complex = dim_C(2) -
    rank(delta^2). The Maurer-Cartan set in the nilpotent regime (MRRS: MC = Z^2)."""
    from quiverlab.fields.linalg import rank
    from quiverlab.resolutions_cs.build import reduction_system_of
    from quiverlab.resolutions_cs.homology import _require_admissible
    from quiverlab.resolutions_cs.resolution import ChouhySolotarResolution
    rs = reduction_system_of(A)
    _require_admissible(rs)
    res = ChouhySolotarResolution(A, rs, max_degree=3, max_cells=max_cells)
    return res.dim_C(2, "coh") - rank(res.matrix(2, "coh"), A.domain)


def maurer_cartan(A, *, order=2, engine="auto", max_cells=4_000_000):
    """The Maurer-Cartan / formal-deformation report on the Hochschild DGLA C(A), solved
    order by order to the certified truncation order `order` (the `order` PARAMETER is the
    truncation limit -- the native CS bracket runs at any degree, so no bracket window
    bounds the reachable order). CHAR-0 ONLY (loud QuiverlabError over char p -- the
    deformation interpretation needs char 0; the field-general HH^2/[a,a] block stays
    available via infinitesimal_deformations / obstruction_map).

    In the MRRS nilpotent regime the MC set equals the 2-cocycles Z^2 (every infinitesimal
    integrates, unobstructed) -- reported directly, no bracket. Outside it, the primary
    obstruction [alpha,alpha] gates the order-2 lift (delta mu_2 = -1/2 [mu_1,mu_1] is
    solvable iff the class vanishes); a complete formal solution past the certified order is
    NEVER claimed."""
    _require_presented(A, "maurer_cartan")
    _require_budget(A, DEFORM_MAXDIM, "maurer_cartan")
    dom = A.domain
    if dom.characteristic != 0:
        raise QuiverlabError(
            "formal-deformation / Maurer-Cartan outputs need characteristic 0 (RRB / MRRS "
            "are char-0 theorems; the L-infinity 1/n! and the DGLA 1/2 need char 0 / char "
            "!= 2); this algebra is over %s" % dom.name,
            hint="the field-general HH^2 / [alpha,alpha] block is still available via "
                 "infinitesimal_deformations(A) / obstruction_map(A)")
    if order < 1:
        raise QuiverlabError("maurer_cartan order must be >= 1 (got %d)" % order)

    nilpotent = is_l_infinity_nilpotent(A)
    table = A.hochschild_cohomology(3, engine=engine, max_cells=max_cells)
    hh2_dim, hh3_dim = table[2], table[3]
    window_note = ("truncation order = %d (the `order` parameter); the native CS bracket "
                   "has no degree window -- higher orders are simply not computed here"
                   % order)

    if nilpotent is True:
        z2 = _z2_dim(A, max_cells=max_cells)
        return MCReport(
            characteristic=0, hh2_dim=hh2_dim, hh3_dim=hh3_dim,
            basis="cs/%s" % dom.name, nilpotent_regime=True, unobstructed=True,
            obstruction_witness=None,
            mc_description=("MC = Z^2 (all %d 2-cocycles integrate; nilpotent regime, MRRS "
                            "Thm 5.4 -- every infinitesimal deformation is unobstructed)" % z2),
            mc_order_certified=order, z2_dim=z2, order_requested=order,
            window_note=window_note, char0_note=None,
            note="nilpotent L-infinity regime certified (gentle, no parallel arrows, no "
                 "oriented cycles); the higher brackets vanish on Z^2.",
            references=tuple(_DEFORM_REFERENCES))

    # non-nilpotent (or not-certified): solve order by order on C(A), gate on [alpha,alpha]
    if hh2_dim == 0:
        # no infinitesimal deformations at all -> trivially unobstructed, MC = {0}
        return MCReport(
            characteristic=0, hh2_dim=0, hh3_dim=hh3_dim, basis="cs/%s" % dom.name,
            nilpotent_regime=nilpotent, unobstructed=True, obstruction_witness=None,
            mc_description="MC = {0}: HH^2 = 0, there are no nontrivial infinitesimal "
                           "deformations to integrate.",
            mc_order_certified=order, z2_dim=_z2_dim(A, max_cells=max_cells),
            order_requested=order, window_note=window_note, char0_note=None,
            note="no infinitesimal deformations (HH^2 = 0); nilpotent regime not certified "
                 "but the deformation functor is trivial.",
            references=tuple(_DEFORM_REFERENCES))

    data = _obstruction_data(A, max_cells=max_cells)
    unobstructed, witness, _ = _verdict(data)
    z2 = data["res"].dim_C(2, "coh") - _rank(data["res"].matrix(2, "coh"), dom)
    if unobstructed:
        certified = min(order, 2)
        desc = ("unobstructed at second order: the primary obstruction [alpha,alpha] "
                "vanishes on all of HH^2, so every infinitesimal deformation lifts to "
                "order 2 (delta mu_2 = -1/2 [mu_1,mu_1] is solvable for every direction)")
        note = ("the primary (order-2) obstruction vanishes; orders > 2 were not computed "
                "(default order = 2). No complete formal solution is claimed.")
    else:
        certified = 1
        desc = ("obstructed: some infinitesimal direction has [alpha,alpha] != 0 in HH^3 "
                "and does NOT lift past first order (the order-2 lift delta mu_2 = "
                "-1/2 [mu_1,mu_1] is unsolvable for the witness direction); directions with "
                "a vanishing self-bracket still lift")
        note = ("a genuine second-order obstruction exists; the DGLA recursion is solvable "
                "iff the right-hand-side class vanishes.")
    return MCReport(
        characteristic=0, hh2_dim=data["hh2_dim"], hh3_dim=data["hh3_dim"],
        basis="cs/%s" % dom.name, nilpotent_regime=nilpotent, unobstructed=unobstructed,
        obstruction_witness=witness, mc_description=desc, mc_order_certified=certified,
        z2_dim=z2, order_requested=order, window_note=window_note, char0_note=None,
        note=note, references=tuple(_DEFORM_REFERENCES + _BRACKET_REFERENCES))


def _rank(M, dom):
    from quiverlab.fields.linalg import rank
    return rank(M, dom)


# ---------------------------------------------------------------------------
# the presented deformed algebra A_alpha (RRRV) + the Ext-algebra handoff (RRR)
# ---------------------------------------------------------------------------
def _deformed_relation(r, pert, t):
    """The deformed relation string r - t*pert (the first-order product a*b = ab +
    t*f(a,b), the radical direction f supported on paths of length >= 2). Folds the sign
    of t into the +/- like QuantumCI (no parenthesised coefficient -- the grammar splits
    every +/-). t == '1' / '-1' drop the redundant coefficient (so the unit direction
    yields 'x*x - 1' with the clean 'term 1 has no arrows' RelationError, Pin 1)."""
    tok = str(t)
    if tok == "1":
        return "%s - %s" % (r, pert)
    if tok == "-1":
        return "%s + %s" % (r, pert)
    if tok.startswith("-"):
        return "%s + %s*%s" % (r, tok[1:], pert)
    return "%s - %s*%s" % (r, tok, pert)


def deformed_algebra(A, direction, *, t="1"):
    """Build the presented deformed algebra A_alpha = kQ/I_alpha for a RADICAL 2-cocycle
    `direction` and value `t`, and RE-CERTIFY it (RRRV, algebraically closed).

    `direction` is a dict {base_relation_str: perturbation_path_str}: each named base
    relation r is replaced by r - t*pert (a radical, length->2 perturbation), the rest of
    the presentation unchanged; the algebra is rebuilt on the SAME quiver and re-certified
    admissible + FLAT (dim A_alpha == dim A).

    Loud on the honest-scope boundaries: a UNIT direction (a length-0 term) is refused at
    PARSE level by RelationError ('term ... has no arrows', Pin 1); a parseable but
    non-admissible / infinite radical direction relays AdmissibilityError /
    NotFiniteDimensionalError; a non-flat jump is refused too. All relayed as QuiverlabError
    (the radical-collapse honest scope, section 5)."""
    from quiverlab.errors import (RelationError, AdmissibilityError,
                                  NotFiniteDimensionalError)
    _require_presented(A, "deformed_algebra")
    if not isinstance(direction, dict) or not direction:
        raise QuiverlabError(
            "deformed_algebra: `direction` must be a non-empty dict "
            "{base_relation_str: perturbation_path_str}",
            hint="e.g. {'x*y': 'y*x'} deforms the relation x*y by - t*(y*x)")
    base = [str(r) for r in A.relations]
    unknown = [k for k in direction if k not in base]
    if unknown:
        raise QuiverlabError(
            "deformed_algebra: direction keys %r are not relations of A (relations: %r)"
            % (unknown, base))
    new_rels = [_deformed_relation(r, direction[r], t) if r in direction else r
                for r in base]
    try:
        A_alpha = A.quiver.algebra(relations=new_rels, field=A.domain)
    except (RelationError, AdmissibilityError, NotFiniteDimensionalError) as exc:
        raise QuiverlabError(
            "deformed_algebra: the direction leaves the admissible radical sub-locus -- "
            "%s: %s. A unit (length-0) direction (e.g. k[x]/(x^2) ~> k[x]/(x^2 - t), the "
            "semisimplification) is an HH^2/MC story on C(A), NOT a presented A_alpha; the "
            "presented feedback covers the RADICAL (admissible) directions only (section 5)."
            % (type(exc).__name__, exc)) from exc
    if A_alpha.dim != A.dim:
        raise QuiverlabError(
            "deformed_algebra: the deformation is NOT flat (dim A_alpha = %d != dim A = %d) "
            "-- the radical changed size; this direction is not a flat formal deformation"
            % (A_alpha.dim, A.dim))
    return A_alpha


def _ext_summary(E):
    """A one-line summary of the Ext-algebra YonedaPresentation (Plan 27) of A_alpha."""
    ngen = sum(len(v) for v in E.generators_by_degree.values())
    nrel = sum(len(v) for v in E.relations_by_degree.values())
    return "E(A_alpha): %d generator%s, %d relation%s, koszul=%s" % (
        ngen, "" if ngen == 1 else "s", nrel, "" if nrel == 1 else "s", E.koszul)


def _parallel_paths(A, k, src, tgt, limit=64):
    """Length-k paths src -> tgt in A's quiver, as '*'-joined arrow-name strings."""
    arrows = A.quiver.arrows
    out = []

    def dfs(cur, word):
        if len(out) >= limit:
            return
        if len(word) == k:
            if cur == tgt:
                out.append("*".join(word))
            return
        for name, (s, tt) in arrows.items():
            if s == cur:
                dfs(tt, word + [name])
    dfs(src, [])
    return out


def _relation_endpoints(A, r):
    """(source, target, length) of a monomial relation string r = 'a*b*...'."""
    arrows = A.quiver.arrows
    names = r.split("*")
    src = arrows[names[0]][0]
    tgt = arrows[names[-1]][1]
    return src, tgt, len(names)


def _canonical_radical_direction(A):
    """A canonical radical deformation direction for the report (display-only): perturb
    one monomial relation r by a parallel length-|r| path p (same source/target) so that
    A_alpha stays FLAT. Best-effort over monomial A; returns ({r: p}, A_alpha) or (None,
    None). The flatness check filters bad perturbations (no separate ideal-membership
    test needed)."""
    if A.quiver is None or A.relations is None:
        return None, None
    mono = [str(r) for r in A.relations if getattr(r, "is_monomial", False)]
    for r in mono:
        try:
            src, tgt, k = _relation_endpoints(A, r)
        except (KeyError, IndexError):
            continue
        if k < 2:
            continue
        for p in _parallel_paths(A, k, src, tgt):
            if p == r:
                continue
            try:
                A_alpha = deformed_algebra(A, {r: p}, t="1")
            except QuiverlabError:
                continue
            return {r: p}, A_alpha
    return None, None


# ---------------------------------------------------------------------------
# Task 4: the B(A)[1] L-infinity companion -- rad^2=0 dg-Lie certificate
# (unconditional) + the FEASIBILITY-GATED char-0 Bardzell l_3 (deferred, see below)
# ---------------------------------------------------------------------------
def _is_monomial(A):
    if A.relations is None:
        return False
    return all(getattr(r, "is_monomial", False) for r in A.relations)


def dg_lie_certificate(A):
    """rad^2 = 0 (monomial) => l_{>=3} == 0, i.e. B(A) is a genuine dg-Lie algebra (RRB).
    A DECIDABLE structural check (rad^2 = 0 iff dim A = |Q_0| + |Q_1|), the one exact
    l_{>=3} claim v1 ships -- UNCONDITIONALLY. Returns:
      * True  -- monomial with rad^2 = 0 (kZ_n/J^2): B(A) is dg-Lie, l_{>=3} == 0;
      * False -- monomial with rad^2 != 0 (QuantumCI(0)): honest L-infinity, NOT dg-Lie;
      * None  -- non-monomial (outside RRB's B(A) scope) or presentation-less.
    (The dg-Lie interpretation is a char-0 theorem; the rad^2=0 fact is structural.)"""
    if A.quiver is None or not _is_monomial(A):
        return None
    rsz = _rad_square_zero(A)
    if rsz is None:
        return None
    return bool(rsz)


# --- FEASIBILITY SPIKE OUTCOME (Task 4 Step 1, recorded here + in the plan Change log) ---
# The bounded spike assessed implementing RRB's l_3 homotopy-transfer formula on a char-0
# Bardzell complex (arXiv:2008.08122). The field-free Bardzell substrate IS reusable over
# QQ (confirmed), but RRB's l_3 requires the FULL homotopy-transfer machinery (the explicit
# contracting homotopy of Bardzell's complex + the transferred tree-summed bracket) -- a
# whole plan's worth of research-grade engine work, beyond one plan (the DECISION-2 freeze
# condition). Per the DECISION, v1 ships:
#   (a) dg_lie_certificate (rad^2=0 => l_{>=3}=0) -- UNCONDITIONAL, above;
#   (b) the induced-l_2 == CS-bracket MODEL-INDEPENDENCE THEOREM STATEMENT (section 4,
#       L-infinity-quasi-iso invariance; rrb_linfty_bardzell) -- CITED, not computed (no
#       B(A) l_2 adapter is built, so there is nothing to compute against; a computed
#       crossengine self-cert is explicitly gated behind the un-built spike, M1 ruling);
#   (c) a LEDGERED deferral of the general l_3/l_4 (DEEPER-ENGINES-BACKLOG.md, P80 reconciles).
# l4 is NEVER claimed zero ("l_n=0 for n>=5" is a sufficient collapse condition, not
# automatic -- honest L-infinity).
_L3_DEFERRAL_NOTE = (
    "the char-0 Bardzell l_3 (RRB, arXiv:2008.08122) is DEFERRED: RRB's homotopy-transfer "
    "l_3 formula is beyond one plan (a full contracting-homotopy + tree-summed transfer on "
    "Bardzell's complex). v1 ships the rad^2=0 dg-Lie certificate (l_{>=3}=0, unconditional) "
    "and the induced-l_2 == CS-bracket model-independence theorem statement (section 4, "
    "cited, not computed -- no B(A) l_2 adapter). l_4 is NEVER claimed zero (honest "
    "L-infinity). The general l_3/l_4 is ledgered (DEEPER-ENGINES-BACKLOG.md, P80).")

_MODEL_INDEPENDENCE_NOTE = (
    "B(A)'s transferred l_2 induces the Gerstenhaber bracket on HH (the L-infinity "
    "quasi-isomorphism B(A) ~ C(A) is a theorem, rrb_linfty_bardzell section 4); this is "
    "stated as the model-independence THEOREM, cited -- not a computed check, since no "
    "B(A) l_2 adapter is built in v1.")


def l3_bracket(A, *, budget=DEFORM_MAXDIM):
    """The feasibility-gated char-0 Bardzell l_3 on B(A)[1] (monomial A, RRB). In v1 this
    is DEFERRED (see the spike outcome): returns an L3Report with status='deferred' and the
    honest note (never a fabricated l_3). Loud QuiverlabError on the scope boundary --
    non-monomial A (RRB's B(A) L-infinity is monomial) or characteristic != 0."""
    _require_presented(A, "l3_bracket")
    if not _is_monomial(A):
        raise QuiverlabError(
            "l3_bracket: RRB's B(A) L-infinity structure is defined for MONOMIAL algebras "
            "only; this algebra is non-monomial",
            hint="the deformation FUNCTOR (HH^2 / obstruction / MC) is served over any "
                 "admissible presentation via obstruction_map / maurer_cartan")
    if A.domain.characteristic != 0:
        raise QuiverlabError(
            "l3_bracket: RRB's L-infinity structure is a characteristic-0 theory "
            "(this algebra is over %s)" % A.domain.name)
    dg = dg_lie_certificate(A)
    return L3Report(
        status="deferred", monomial=True, rad_square_zero=_rad_square_zero(A), dg_lie=dg,
        note=_L3_DEFERRAL_NOTE + (
            " (This algebra IS rad^2=0 monomial, so l_{>=3}=0 by the dg-Lie certificate -- "
            "the deferral concerns the GENERAL monomial l_3, not this collapsed case.)"
            if dg else ""),
        references=tuple(["rrb_linfty_bardzell", "mrrs_mc_gentle"]))


def _l3_status(A, dg, char):
    """The l3_status string for the report/block: honest scope + spike outcome."""
    if A.quiver is None or not _is_monomial(A):
        return "n/a (non-monomial)"
    if char != 0:
        return "n/a (char p)"
    if dg is True:
        return "dg_lie"          # rad^2=0 => l_{>=3}=0 (unconditional certificate)
    return "deferred"            # honest L-infinity; the general l_3 is ledgered


# ---------------------------------------------------------------------------
# the full report + the runner block
# ---------------------------------------------------------------------------
_BASE_CHANGE_NOTE = (
    "RRRV's presentation of A_alpha by quiver and relations holds over an algebraically "
    "closed field; the build here is per-instance certified (flat: dim A_alpha = dim A).")


def deformation_structure(A, *, budget=DEFORM_MAXDIM, engine="auto",
                          max_cells=4_000_000):
    """The full formal-deformation report (frozen Deformations). Char 0: infinitesimal
    HH^2, the primary obstruction [alpha,alpha] in HH^3 (+ witness), the MRRS nilpotent
    verdict, the MC description + certified order, the rad^2=0 dg-Lie certificate + l_3
    status, and a canonical radical A_alpha (display-only) with its Ext-algebra summary.
    Char p: the field-general HH^2/HH^3 dims + a char0_note (the deformation interpretation
    is char-0 gated); the expensive obstruction bracket is skipped (its verdict is a char-0
    notion). Loud on oversize / presentation-less."""
    _require_presented(A, "deformation_structure")
    _require_budget(A, budget, "deformation_structure")
    dom = A.domain
    char = dom.characteristic
    window_note = ("MC truncation order is bounded by the `order` parameter of "
                   "maurer_cartan (default 2 -- the primary obstruction); the native CS "
                   "bracket has no degree window.")

    if char != 0:
        table = A.hochschild_cohomology(3, engine=engine, max_cells=max_cells)
        return Deformations(
            characteristic=char, hh2_dim=table[2], hh3_dim=table[3],
            basis="cs/%s" % dom.name, unobstructed=None, obstruction_witness=None,
            nilpotent_regime=None,
            mc_description="the Maurer-Cartan / L-infinity interpretation is "
                           "characteristic-0 gated; only the field-general HH^2 / HH^3 "
                           "dimensions are reported here.",
            mc_order_certified=1, dg_lie=None, l3_status="n/a (char p)",
            a_alpha=None, ext_algebra_summary=None, base_change_note=None,
            char0_note=_char0_note(A), window_note=None, status="complete",
            note="field-general HH block (characteristic %d); "
                 "the deformation theory is char-0." % char,
            references=tuple(_DEFORM_REFERENCES))

    # char 0: the full report
    data = _obstruction_data(A, max_cells=max_cells)
    unobstructed, witness, _ = _verdict(data)
    hh2, hh3 = data["hh2_dim"], data["hh3_dim"]
    nilpotent = is_l_infinity_nilpotent(A)
    if nilpotent is True:
        z2 = data["res"].dim_C(2, "coh") - _rank(data["res"].matrix(2, "coh"), dom)
        mc_desc = ("MC = Z^2 (all %d 2-cocycles integrate; nilpotent regime, MRRS Thm 5.4 "
                   "-- every infinitesimal deformation is unobstructed)" % z2)
        mc_order = 2
    elif hh2 == 0:
        mc_desc = ("MC = {0}: HH^2 = 0, no nontrivial infinitesimal deformations")
        mc_order = 2
    elif unobstructed:
        mc_desc = ("unobstructed at second order: the primary obstruction [alpha,alpha] "
                   "vanishes on all of HH^2 (a second-order lift exists for every direction)")
        mc_order = 2
    else:
        mc_desc = ("obstructed: some direction has [alpha,alpha] != 0 in HH^3 and does not "
                   "lift past first order (witness %s); directions with a vanishing "
                   "self-bracket still lift" % witness)
        mc_order = 1
    dg = dg_lie_certificate(A)
    l3_status = _l3_status(A, dg, char)
    direction, A_alpha = _canonical_radical_direction(A)
    a_alpha = None
    ext_summary = None
    if A_alpha is not None:
        a_alpha = {
            "quiver": {"vertices": list(A_alpha.quiver.vertices),
                       "arrows": {n: [s, t] for n, (s, t) in A_alpha.quiver.arrows.items()}},
            "relations": [str(r) for r in A_alpha.relations],
            "direction": dict(direction), "t": "1", "flat": True}
        try:
            ext_summary = _ext_summary(A_alpha.ext_algebra(3))
        except QuiverlabError:
            ext_summary = None
    return Deformations(
        characteristic=0, hh2_dim=hh2, hh3_dim=hh3, basis="cs/%s" % dom.name,
        unobstructed=unobstructed, obstruction_witness=witness,
        nilpotent_regime=nilpotent, mc_description=mc_desc, mc_order_certified=mc_order,
        dg_lie=dg, l3_status=l3_status, a_alpha=a_alpha, ext_algebra_summary=ext_summary,
        base_change_note=(_BASE_CHANGE_NOTE if a_alpha else None), char0_note=None,
        window_note=window_note, status="complete",
        note=(_MODEL_INDEPENDENCE_NOTE + " " + _L3_DEFERRAL_NOTE
              if l3_status == "deferred" else _MODEL_INDEPENDENCE_NOTE),
        references=tuple(_DEFORM_REFERENCES + _BRACKET_REFERENCES))


def deformations_block(A, *, budget=DEFORM_MAXDIM):
    """The JSON block for the `deformations` compute kind (Plan 78), shared by all three
    tiers (byte-identical). The raw bracket constants are NOT shipped (they explode and are
    basis-dependent); the block reports HH^2/HH^3, the obstruction verdict + witness, the
    nilpotent verdict, the MC description, the dg-Lie certificate + l_3 status, and the
    display-only presented A_alpha. A char-p / oversize / presentation-less refusal is a
    clean {"error": ...} block, never a 500."""
    refs = list(_DEFORM_REFERENCES)
    try:
        D = deformation_structure(A, budget=budget)
    except QuiverlabError as exc:
        return {"kind": "deformations", "status": "budget", "error": str(exc),
                "references": refs}
    return {
        "kind": "deformations", "characteristic": D.characteristic,
        "hh2_dim": D.hh2_dim, "hh3_dim": D.hh3_dim,
        "unobstructed": D.unobstructed, "obstruction_witness": D.obstruction_witness,
        "nilpotent_regime": D.nilpotent_regime, "mc_description": D.mc_description,
        "mc_order_certified": D.mc_order_certified, "dg_lie": D.dg_lie,
        "l3_status": D.l3_status, "a_alpha": D.a_alpha,
        "ext_algebra_summary": D.ext_algebra_summary,
        "base_change_note": D.base_change_note, "char0_note": D.char0_note,
        "window_note": D.window_note, "status": D.status, "note": D.note,
        "references": refs,
    }
