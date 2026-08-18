"""The Han bounded-extension recognizer + Jacobi-Zariski Han transport (Plan 73, R7).

Cibils-Lanzilotta-Marcos-Solotar, ``2101.02597`` (J. Algebra 2022). An extension
``B subset A`` is **left (resp. right) bounded** iff the ``B``-bimodule ``A/B`` is
(i) ``B``-tensor nilpotent, (ii) of finite projective dimension over ``B^e``, and
(iii) left (resp. right) ``B``-projective. Under a bounded extension Han's
conjecture transports: ``B |= Han  <=>  A |= Han`` (Thm 4.6).

This module houses the three legs of the recognizer and the transport ladder. Leg
(i) -- the load-bearing one -- is an **honest capped semi-decision** on the direct
tensor powers ``(A/B)^{ox_B m}`` PLUS a certificate route (relative cycles +
``J``-interrupters); ``"undecided"`` is a first-class verdict, never a guess.

Composition is left-to-right (quiverlab); the CLMS paper is right-to-left.
"""
from dataclasses import dataclass, field as _field

from quiverlab.errors import QuiverlabError
from quiverlab.families.extension import (
    arrow_removal_subalgebra, enveloping_algebra, _f_length, _tokens)


# --------------------------------------------------------------------------- #
# leg (i): tensor nilpotency (CLMS Def. 2.3, Sec. 5)
# --------------------------------------------------------------------------- #
@dataclass
class TensorNilpotency:
    """A data report for leg (i). Read ``status``; no ``__bool__``."""
    status: str                          # "nilpotent" | "not_nilpotent" | "undecided"
    index: "int | None" = None           # nilpotency index n when "nilpotent"
    route: str = ""                      # "no_relative_cycles"|"J_interrupter"|"direct_cap"|"witness"
    witness: "object | None" = None      # the surviving relative cycle when "not_nilpotent"
    cap: int = 0
    note: str = ""


def is_tensor_nilpotent(ext, *, cap=8, use_certificate=True):
    """Decide whether ``A/B`` is ``B``-tensor nilpotent (CLMS Def. 2.3).

    The verdict is **authoritative from the direct tensor route** and the
    length-index theorem (CLMS Def. 5.19: if nilpotent the index ``<=`` the length
    index ``L``, so ``(A/B)^{ox L} = 0`` iff nilpotent) whenever ``|R|^L`` is within
    the cell guard; the relative cycles + ``J``-interrupters supply the ``route``
    label and the non-nilpotency ``witness``.

    ``use_certificate=False`` forces the **naive capped semi-decision** (the record's
    honest capped route): compute ``(A/B)^{ox m}`` for ``m = 1..cap``; a vanishing
    power ``=> "nilpotent"``; otherwise ``"undecided"`` -- it NEVER concludes
    ``"not_nilpotent"`` (powers can grow before dying).
    """
    L = ext.length_index()

    # ---- naive capped semi-decision (no length-index shortcut) -------------- #
    if not use_certificate:
        return _capped_semidecision(ext, cap)

    cycles = ext.relative_cycles(cap=cap)

    # ---- authoritative direct decision via the length-index theorem --------- #
    if ext._tensor_feasible(L):
        index = None
        for m in range(1, L + 1):
            if ext.tensor_power_dim(m) == 0:
                index = m
                break
        if index is not None:
            route = "no_relative_cycles" if not cycles else "J_interrupter"
            note = (f"no relative cycles => nilpotent (CLMS Thm 5.14)"
                    if not cycles else
                    f"every relative cycle is J-interrupted => nilpotent (CLMS Thm 5.16)")
            return TensorNilpotency("nilpotent", index=index, route=route, cap=cap,
                                    note=f"{note}; length index L={L}")
        # (A/B)^{ox L} != 0 => NOT nilpotent (length-index theorem, rigorous)
        witness = _surviving_cycle(ext, cycles)
        return TensorNilpotency(
            "not_nilpotent", index=None, route="witness", witness=witness, cap=cap,
            note=f"(A/B)^(ox {L}) != 0 at the length index L={L} => not nilpotent "
                 f"(CLMS Def. 5.19); surviving relative cycle {witness!r}")

    # ---- length index infeasible: fall to the honest capped route ----------- #
    tn = _capped_semidecision(ext, cap)
    tn.note = (tn.note + f"; length index L={L} too large for the direct decision "
                         f"(|R|^L past the cell guard)")
    return tn


def _capped_semidecision(ext, cap):
    for m in range(1, cap + 1):
        try:
            d = ext.tensor_power_dim(m)
        except QuiverlabError:
            return TensorNilpotency(
                "undecided", index=None, route="direct_cap", cap=cap,
                note=f"tensor power m={m} exceeds the cell guard; undecided")
        if d == 0:
            return TensorNilpotency("nilpotent", index=m, route="direct_cap", cap=cap,
                                    note=f"(A/B)^(ox {m}) = 0 within the cap")
    return TensorNilpotency(
        "undecided", index=None, route="direct_cap", cap=cap,
        note=f"{cap} tensor powers nonzero; cannot conclude (powers may grow before "
             f"dying) -- undecided (no a-priori index bound)")


def _surviving_cycle(ext, cycles):
    """A relative cycle whose self-tensor survives (no ``J``-interrupter) -- the
    non-nilpotency witness (CLMS Ex. 5.5). Falls back to any relative cycle."""
    for w in cycles:
        if not ext.has_J_interrupter(w):
            return w
    return cycles[0] if cycles else None


# --------------------------------------------------------------------------- #
# leg (ii): finite pd_{B^e}(A/B) -- gl.dim primary (CLMS Ex. 6.1) + fallback
# --------------------------------------------------------------------------- #
def finite_pd_Be(ext, *, pd_cap=16):
    """Is ``pd_{B^e}(A/B) < infinity``? Returns a dict
    ``{route, status, gldim_B|value, note}``.

    **PRIMARY (free) route (CLMS Ex. 6.1):** if ``gl.dim B < infinity`` then
    ``gl.dim B^e < infinity`` and every ``B``-bimodule -- in particular ``A/B`` --
    has finite ``pd_{B^e}``. No enveloping algebra is built. **FALLBACK
    (``gl.dim B = infinity``):** ``A/B`` as a right ``B^e``-module, ``pd`` via the
    shipped module resolution capped at ``pd_cap`` (Task II3); honest ``"undecided"``
    + certified lower bound over the cap (never ``infinity`` unproven).
    """
    B = ext.B
    gd = B.global_dimension()
    if gd.exact:
        v = int(gd)
        return {"route": "gldim", "status": "finite", "gldim_B": v, "value": None,
                "note": f"gl.dim B = {v} < inf => gl.dim B^e < inf => "
                        f"pd_{{B^e}}(A/B) < inf (CLMS Ex. 6.1)"}
    return _finite_pd_Be_envelope(ext, pd_cap=pd_cap)


def _finite_pd_Be_envelope(ext, *, pd_cap=16):
    """The ``gl.dim B = infinity`` fallback: build ``B^e`` and resolve ``A/B`` as a
    right ``B^e``-module, capped at ``pd_cap``. Finite ``pd`` ``=> "finite"``;
    unresolved within the cap ``=> "undecided"`` + the reached-length certified lower
    bound (never ``infinity`` unproven). Loud refusals (heavy product quiver over
    budget) surface as ``"undecided"`` with the reason."""
    from quiverlab.hochschild.coefficients import Bimodule
    B = ext.B
    try:
        Be = enveloping_algebra(B)
        M = ext.quotient_bimodule()
        N = _bimodule_to_Be_module(M, B, Be)
        pd = N.projective_resolution(pd_cap).pd()
    except QuiverlabError as exc:
        return {"route": "enveloping", "status": "undecided", "gldim_B": None,
                "value": None,
                "note": f"gl.dim B = inf; enveloping fallback unavailable: {exc}"}
    if pd is not None:
        return {"route": "enveloping", "status": "finite", "gldim_B": None,
                "value": pd,
                "note": f"pd_{{B^e}}(A/B) = {pd} via B^e = B (x) B^op (gl.dim B = inf)"}
    return {"route": "enveloping", "status": "undecided", "gldim_B": None,
            "value": None,
            "note": f"pd_{{B^e}}(A/B) not resolved within pd_cap={pd_cap}: certified "
                    f"lower bound pd > {pd_cap} (never infinity unproven)"}


# --------------------------------------------------------------------------- #
# module bridges: a B-bimodule / A/B as a one-sided or B^e module
# --------------------------------------------------------------------------- #
def _bimodule_to_Be_module(M, B, Be):
    """A ``B``-bimodule ``M`` as a **right** ``B^e``-module over ``Be = B (x) B^op``
    (``m . (b (x) c^op) = c . m . b``): the ``Be``-vertex ``(u, v)`` is the corner
    ``e_v . M . e_u`` (right ``u``, left ``v``); ``L_{al}_v = al (x) e_v`` acts by
    right-mult ``m . al``; ``R_{be}_u = e_u (x) be^op`` acts by left-mult ``be . m``.
    Validated by ``from_arrow_action`` (the module axioms are the correctness gate)."""
    from quiverlab.modules.module import Module
    dom = B.domain

    def which_vertex(s, apply):
        es = M._unit_vec(s)
        hits = [v for v in B.quiver.vertices
                if apply(B._basis_vec(B.basis_labels.index(f"e_{v}")), es) == es]
        return hits[0] if len(hits) == 1 else None

    # corner (u, v) = e_v . M . e_u  (right vertex u, left vertex v)
    corner = {}
    for s in range(M.dim_M):
        lv = which_vertex(s, M.left_apply)
        rv = which_vertex(s, M.right_apply)
        if lv is None or rv is None:
            raise QuiverlabError(
                "A/B is not B-homogeneous: a basis vector lies in no single corner "
                "e_v M e_u -- cannot present it as a right B^e-module",
                hint="the arrow-extension bimodules are corner-homogeneous")
        corner.setdefault((rv, lv), []).append(s)

    verts = list(Be.quiver.vertices)
    dimvec = {vv: len(corner.get(vv, [])) for vv in verts}
    starts, off = {}, 0
    for vv in verts:
        starts[vv] = off
        off += dimvec[vv]
    n = off
    gidx = {}
    for key, lst in corner.items():
        for p, s in enumerate(lst):
            gidx[s] = starts[key] + p

    def bvec(name):
        return B._basis_vec(B.basis_labels.index(name))

    arrow_action = {}
    for nm, (sv, tv) in Be.quiver.arrows.items():
        mat = [[dom.zero()] * n for _ in range(n)]
        gen = nm[2:].rsplit("_", 1)[0]
        if nm.startswith("L_"):                       # al (x) e_v -> m . al
            action = [M.right_apply(bvec(gen), M._unit_vec(s)) for s in range(M.dim_M)]
        else:                                         # e_u (x) be^op -> be . m
            action = [M.left_apply(bvec(gen), M._unit_vec(s)) for s in range(M.dim_M)]
        for s in corner.get(sv, []):
            res = action[s]
            for t in range(M.dim_M):
                if not dom.is_zero(res[t]) and t in gidx \
                        and starts[tv] <= gidx[t] < starts[tv] + dimvec[tv]:
                    mat[gidx[t]][gidx[s]] = dom.add(mat[gidx[t]][gidx[s]], res[t])
        arrow_action[nm] = mat
    return Module.from_arrow_action(Be, dimvec, arrow_action, name="A/B over B^e")


def _one_sided_module(ext, side):
    """``A/B`` as a right (``side="right"``) or left (``side="left"``) ``B``-module.
    A left ``B``-module is a right ``B^op``-module (built over ``B.opposite()``)."""
    from quiverlab.modules.module import Module
    A = ext.A
    dom = A.domain
    labels = A.basis_labels
    Q = A.quiver
    rel = list(ext.rel_idx)
    nR = len(rel)
    rpos = {ai: k for k, ai in enumerate(rel)}

    def projR(vec):
        return [vec[rel[k]] for k in range(nR)]

    if side == "right":
        rep = ext.B                                   # right B-module
        vfun = lambda ai: Q.word_target(_tokens(labels[ai]))

        def act(ai, arrow):                           # r . arrow
            return projR(A.multiply(A._basis_vec(ai), A._basis_vec(labels.index(arrow))))
    elif side == "left":
        rep = ext.B.opposite()                        # left B-module = right B^op-module
        vfun = lambda ai: Q.word_source(_tokens(labels[ai]))

        def act(ai, arrow):                           # arrow . r  (left-mult in A)
            return projR(A.multiply(A._basis_vec(labels.index(arrow)), A._basis_vec(ai)))
    else:
        raise QuiverlabError(f"side must be 'left' or 'right', got {side!r}")

    verts = list(rep.quiver.vertices)
    corner = {v: [] for v in verts}
    for ai in rel:
        corner[vfun(ai)].append(ai)
    dimvec = {v: len(corner[v]) for v in verts}
    starts, off = {}, 0
    for v in verts:
        starts[v] = off
        off += dimvec[v]
    n = off
    gidx = {}
    for v in verts:
        for p, ai in enumerate(corner[v]):
            gidx[ai] = starts[v] + p
    arrow_action = {}
    for arrow, (sv, tv) in rep.quiver.arrows.items():
        mat = [[dom.zero()] * n for _ in range(n)]
        for ai in corner[sv]:
            res = act(ai, arrow)
            for k in range(nR):
                if not dom.is_zero(res[k]) and vfun(rel[k]) == tv:
                    mat[gidx[rel[k]]][gidx[ai]] = dom.add(mat[gidx[rel[k]]][gidx[ai]], res[k])
        arrow_action[arrow] = mat
    return Module.from_arrow_action(rep, dimvec, arrow_action, name=f"A/B {side}")


# --------------------------------------------------------------------------- #
# leg (iii): one-sided B-projectivity (CLMS Thm 5.20)
# --------------------------------------------------------------------------- #
@dataclass
class OneSidedProjectivity:
    """Leg (iii): is ``A/B`` projective as a one-sided ``B``-module? A data report."""
    side: "str | None"                   # the projective side ("left"|"right") or None
    projective: "bool | None"            # True if projective on `side` (or either, auto)
    pd_left: "int | None" = None
    pd_right: "int | None" = None
    note: str = ""


def one_sided_projective(ext, *, side="auto", bound=16):
    """Decide whether ``A/B`` is projective as a one-sided ``B``-module (leg iii,
    CLMS Def. 2.3 / Thm 5.20). ``pd_B(A/B) = 0`` on the shipped module resolution
    stack is the exact test.

    ``side="auto"`` tries both sides and reports the projective one (CLMS bounds an
    extension on a FIXED side; either suffices for the Han transport, Thm 4.6).
    NOTE the composition-convention translation: quiverlab composes left-to-right,
    the CLMS paper right-to-left, so a CLMS "left"-bounded extension is
    "right"-projective here (Ex. 5.3 is RIGHT ``B``-projective in this convention)."""
    def pd_side(sd):
        M = _one_sided_module(ext, sd)
        return M.projective_resolution(bound).pd()

    if side in ("left", "right"):
        pd = pd_side(side)
        proj = (pd == 0)
        rep = {("left" if side == "left" else "right"): pd}
        return OneSidedProjectivity(
            side=(side if proj else None), projective=proj,
            pd_left=pd if side == "left" else None,
            pd_right=pd if side == "right" else None,
            note=f"pd_B(A/B as a {side} B-module) = {pd}")
    # auto: try right then left
    pdr = pd_side("right")
    pdl = pd_side("left")
    if pdr == 0:
        winner = "right"
    elif pdl == 0:
        winner = "left"
    else:
        winner = None
    return OneSidedProjectivity(
        side=winner, projective=(winner is not None),
        pd_left=pdl, pd_right=pdr,
        note=f"pd_B(A/B) left={pdl}, right={pdr}; projective side = {winner}")


# --------------------------------------------------------------------------- #
# the bounded-extension certificate (three legs)
# --------------------------------------------------------------------------- #
@dataclass
class BoundedCertificate:
    """The three-leg bounded-extension certificate (CLMS Def. 2.3). A data report."""
    bounded: "bool | None"               # True (a side) / False / None (undecided)
    side: "str | None"                   # "left" | "right" | None
    tensor_nilpotent: TensorNilpotency
    one_sided: OneSidedProjectivity
    pd_Be: dict                          # {route, status, gldim_B|value, ...}
    note: str = ""


def bounded_extension(A, new_arrows, *, side="auto", nilp_cap=8, pd_cap=16):
    """Decide whether ``B subset A`` (by removing ``new_arrows``) is a **bounded**
    extension: ``A/B`` (i) ``B``-tensor nilpotent, (ii) finite ``pd_{B^e}``, (iii)
    one-sided ``B``-projective (CLMS Def. 2.3). Bounded ``<=>`` all three on a FIXED
    side. Returns a :class:`BoundedCertificate` -- honest ``None`` when a leg is
    undecided."""
    ext = arrow_removal_subalgebra(A, new_arrows)
    tn = is_tensor_nilpotent(ext, cap=nilp_cap)
    if tn.status == "not_nilpotent":
        # leg (i) fails => not a bounded extension via this route
        osp = OneSidedProjectivity(side=None, projective=False, note="leg (i) False")
        return BoundedCertificate(
            bounded=False, side=None, tensor_nilpotent=tn, one_sided=osp,
            pd_Be={"route": None, "status": "n/a"},
            note="A/B is not tensor-nilpotent (leg i) -- not bounded (CLMS Ex. 5.5)")
    if tn.status == "undecided":
        osp = OneSidedProjectivity(side=None, projective=None, note="leg (i) undecided")
        return BoundedCertificate(
            bounded=None, side=None, tensor_nilpotent=tn, one_sided=osp,
            pd_Be={"route": None, "status": "undecided"},
            note="leg (i) tensor-nilpotency undecided (cap reached)")
    # leg (i) nilpotent -> the other two legs
    pd = finite_pd_Be(ext, pd_cap=pd_cap)
    osp = one_sided_projective(ext, side=side)
    pd_finite = True if pd["status"] == "finite" else None
    if osp.projective and pd_finite:
        return BoundedCertificate(
            bounded=True, side=osp.side, tensor_nilpotent=tn, one_sided=osp, pd_Be=pd,
            note=f"bounded on the {osp.side} side: A/B tensor-nilpotent (index "
                 f"{tn.index}), pd_{{B^e}} finite ({pd['route']}), {osp.side}-projective")
    reason = []
    if not osp.projective:
        reason.append("A/B not one-sided B-projective")
    if not pd_finite:
        reason.append(f"finite pd_{{B^e}} {pd['status']}")
    return BoundedCertificate(
        bounded=(None if (pd_finite is None) else False),
        side=None, tensor_nilpotent=tn, one_sided=osp, pd_Be=pd,
        note="not bounded: " + "; ".join(reason))


# --------------------------------------------------------------------------- #
# Han transport -- the injection/iso LADDER (CLMS Thm 3.1/4.6, H1)
# --------------------------------------------------------------------------- #
def transport_verdict(nilp, pd_finite, one_sided):
    """The pure ladder function (H1). ``nilp``/``pd_finite``/``one_sided`` are
    ``True``/``False``/``None`` (None = undecided). There is NO "iso from leg (i)
    alone" verdict."""
    if nilp is False:
        return "not_bounded"
    if nilp is None:
        return "undecided"
    if pd_finite is True:
        return "bounded" if one_sided is True else "pd_injection"
    return "nilpotent_injection"


@dataclass
class HanTransport:
    """The Han-transport certificate -- the injection/iso LADDER (H1). A data report."""
    transport: str                       # bounded|pd_injection|nilpotent_injection|not_bounded|undecided
    certificate: BoundedCertificate
    han_B: "bool | None" = None
    han_A: "bool | None" = None
    injection_from: "int | None" = None  # CERTIFIED LOWER BOUND (H3), not a precise threshold
    injection_bound_ok: "bool | None" = None
    references: list = _field(default_factory=list)
    note: str = ""


_HAN_REFERENCES = ["clms_bounded_extensions", "clms_jacobi_zariski",
                   "kaygun_jacobi_zariski", "han_conjecture", "assem_book"]


def han_transport(A, new_arrows, *, side="auto", nilp_cap=8, pd_cap=16, hh_top=None):
    """Transport Han's conjecture across ``B subset A`` and label the claim by the
    exact row of the injection/iso ladder (H1). The self-cert gate is the **injection
    bound** ``dim HH_m(B) <= dim HH_m(A)`` on the shipped HH engines, with **equality
    asserted only under ``"bounded"``**."""
    cert = bounded_extension(A, new_arrows, side=side, nilp_cap=nilp_cap, pd_cap=pd_cap)
    tn = cert.tensor_nilpotent
    nilp = {"nilpotent": True, "not_nilpotent": False, "undecided": None}[tn.status]
    # True only when pd_{B^e}(A/B) is CERTIFIED finite; otherwise None = "unknown",
    # deliberately NOT False -- an uncertified or infinite-looking pd must never let
    # `transport_verdict` claim a row, it must fall through to the honest "undecided"
    # end of the ladder. (This was written as `None if nilp else None`, a no-op ternary
    # that read as if one branch were meant to be False; the conservative None is what
    # the ladder actually relies on, so it is now spelled once.)
    pd_finite = True if cert.pd_Be.get("status") == "finite" else None
    one_sided = cert.one_sided.projective if nilp else None
    transport = transport_verdict(nilp, pd_finite, one_sided)

    top = hh_top if hh_top is not None else 6
    engine = "cs" if A.domain.name.startswith("QQ") else "auto"
    ext = arrow_removal_subalgebra(A, new_arrows)
    hh_A = A.hochschild_homology(top, engine=engine, verbose=False).dims
    hh_B = ext.B.hochschild_homology(top, engine=engine, verbose=False).dims
    inj_ok = all(b <= a for a, b in zip(hh_A, hh_B))
    equality = (hh_A == hh_B)
    if transport == "bounded" and not equality:
        # the iso row demands equality for * >> 0; a mismatch is a loud drift
        raise QuiverlabError(
            "han_transport self-cert FAILED: transport='bounded' asserts "
            f"HH_*(B) ~= HH_*(A) for * >> 0, but HH_*(B)={hh_B} != HH_*(A)={hh_A}",
            hint="a bounded extension forces the isomorphism (CLMS Thm 4.6)")

    # Han per side, honestly scaled (claim the iff only under "bounded")
    han_B = True if cert.pd_Be.get("route") == "gldim" else None
    han_A = han_B if transport == "bounded" else None

    injection_from = tn.index          # certified LOWER BOUND (H3), not a precise degree
    claim = {
        "bounded": "HH_*(B) ~= HH_*(A) (iso) => full iff B |= Han <=> A |= Han (Thm 4.6)",
        "pd_injection": "HH_*(B) ↪ HH_*(A) (injection, ordinary coeff)",
        "nilpotent_injection": "H_*(B,A) ↪ H_*(A,A) (injection, coefficients in A)",
        "not_bounded": "not a bounded extension via this route (Ex. 5.5)",
        "undecided": "honest cap reached; per-leg status reported",
    }[transport]
    return HanTransport(
        transport=transport, certificate=cert, han_B=han_B, han_A=han_A,
        injection_from=injection_from, injection_bound_ok=inj_ok,
        references=list(_HAN_REFERENCES),
        note=f"{claim}; injection bound dim HH_m(B) <= dim HH_m(A) held: {inj_ok}"
             + ("; EQUALITY (bounded => iso)" if equality and transport == "bounded" else ""))


# --------------------------------------------------------------------------- #
# the GUI / webapp block (algebra-level certificate carrying the new-arrow subset)
# --------------------------------------------------------------------------- #
def han_transport_block(A, new_arrows, *, hh_top=None):
    """The ``han_transport`` compute-kind block (all three tiers, shared by both
    runners). Loud refusals (structure-constant ``A`` / bad arrow name) become a
    clean ``error`` block -- never a 500."""
    block = {"kind": "han_transport", "new_arrows": list(new_arrows),
             "references": list(_HAN_REFERENCES)}
    try:
        ht = han_transport(A, tuple(new_arrows), hh_top=hh_top)
    except QuiverlabError as exc:
        block.update(status="unsupported", error=str(exc),
                     note="Han transport unavailable for this input")
        return block
    cert = ht.certificate
    tn = cert.tensor_nilpotent
    ext = arrow_removal_subalgebra(A, tuple(new_arrows))    # dims (cheap arrow-removal)
    block.update(
        status="complete",
        dim_A=ext.dim_A, dim_B=ext.dim_B, dim_quotient=ext.dim_quotient,
        transport=ht.transport,
        tensor_nilpotent={"status": tn.status, "index": tn.index, "route": tn.route},
        one_sided={"side": cert.one_sided.side, "projective": cert.one_sided.projective},
        pd_Be={"route": cert.pd_Be.get("route"), "value": cert.pd_Be.get("value"),
               "status": cert.pd_Be.get("status"),
               "gldim_B": cert.pd_Be.get("gldim_B")},
        injection_from=ht.injection_from,
        injection_bound_ok=ht.injection_bound_ok,
        han_B=ht.han_B, han_A=ht.han_A,
        note=ht.note,
    )
    return block
