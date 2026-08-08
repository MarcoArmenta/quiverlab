"""Shared no-code block builders for pi1(Q, I) + simple connectivity (Plan 56).

BOTH runners -- the server (quiverlab.hpc.spec) and the Pyodide twin
(docs/gui/runner) -- import these, so the blocks are byte-identical by
construction; each runner only appends the resolved ``citations`` locally.
"""
from quiverlab.invariants.coverings import fundamental_group, is_simply_connected

_FG_REFS = ["assem_delapena", "martinez_villa_delapena", "crs_hurewicz",
            "crs_gradings", "briggs_ryd_tori"]
_SC_REFS = ["assem_book", "skowronski_ssc", "le_meur_pi1"]

_INTRINSIC_NOTE = (
    "The INTRINSIC fundamental group (the inverse limit over connected gradings) is "
    "not bounded-computable and is NOT reported here: e.g. the presentation "
    "pi1(k[x]/(x^p)) = Z differs from the intrinsic pi1 = Z x C_p in characteristic p "
    "(Cibils-Redondo-Solotar). quiverlab reports the PRESENTATION group pi1(Q, I).")


def _abelianization_latex(g) -> str:
    parts = []
    if g.free_rank == 1:
        parts.append(r"\mathbb{Z}")
    elif g.free_rank > 1:
        parts.append(r"\mathbb{Z}^{%d}" % g.free_rank)
    parts += [r"\mathbb{Z}/%d" % d for d in g.invariant_factors]
    return r" \oplus ".join(parts) if parts else "0"


def _jsonable(x):
    """Recursively turn sets/frozensets into sorted lists so the block is JSON-safe
    and deterministic (sorted by str to tolerate mixed vertex types)."""
    if isinstance(x, (set, frozenset)):
        return sorted((_jsonable(v) for v in x), key=str)
    if isinstance(x, dict):
        return {k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    return x


def fundamental_group_block(A) -> dict:
    """The pi1(Q, I) block: generators, relators, abelianization (free rank +
    invariant factors), Hom(pi1, k+) dim, components, and the honest intrinsic note.
    The HH^1 bound is NOT computed here (it stays a verification-side cross-check)."""
    g = fundamental_group(A)
    char = A.domain.characteristic
    return {
        "kind": "fundamental_group",
        "generators": list(g.generators),
        "relators": list(g.relators),
        "abelianization": {"free_rank": g.free_rank,
                           "invariant_factors": list(g.invariant_factors)},
        "hom_to_additive_dim": g.hom_to_additive_dim(char),
        "components": g.components,
        "presentation_note": g.presentation_note,
        "latex": r"\pi_1^{\mathrm{ab}} = " + _abelianization_latex(g),
        "intrinsic_note": _INTRINSIC_NOTE,
        "references": list(_FG_REFS),
    }


def simply_connected_block(A) -> dict:
    """The three-valued simple-connectivity block WITH the strongly-simply-connected
    certificate (strong=True, so the P62 consumable is always present). A char/budget
    event surfaces as strongly.verdict = null with its reason, never an error."""
    r = is_simply_connected(A, strong=True)
    strong = r.strongly
    strong_blk = None
    if strong is not None:
        strong_blk = {
            "verdict": strong.verdict,
            "witness": _jsonable(strong.witness),
            "reason": strong.reason,
            "checked_convex": strong.checked_convex,
        }
    return {
        "kind": "simply_connected",
        "verdict": r.verdict,
        "reason": r.reason,
        "witness": _jsonable(r.witness),
        "strongly": strong_blk,
        "references": list(_SC_REFS),
    }
