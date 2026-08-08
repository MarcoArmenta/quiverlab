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
    coordinate vector over the recorded HH^3 basis: sum_{i,j} coords_i coords_j B[i][j].
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
