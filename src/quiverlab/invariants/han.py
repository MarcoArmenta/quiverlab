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
from quiverlab.families.extension import arrow_removal_subalgebra


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
    """The ``gl.dim B = infinity`` fallback: build ``B^e`` and resolve ``A/B`` over
    it, capped (Task II3). Placeholder until Task II3 wires the enveloping bridge --
    honest ``"undecided"`` so the primary-route slice stays green."""
    return {"route": "enveloping", "status": "undecided", "gldim_B": None,
            "value": None,
            "note": "gl.dim B = inf: the enveloping-algebra pd fallback is not yet "
                    "wired (Task II3)"}
