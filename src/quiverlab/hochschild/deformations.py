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
