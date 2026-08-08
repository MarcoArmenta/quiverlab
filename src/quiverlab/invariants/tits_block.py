"""Shared no-code block builder for the Tits-form tame/wild certificate (Plan 62 /
R19).

BOTH runners -- the server (quiverlab.hpc.spec) and the Pyodide twin
(docs/gui/runner) -- import this, so the block is byte-identical by construction;
each runner only appends the resolved ``citations`` locally.

A presentation-less algebra or a non-triangular (loop / oriented-cycle) quiver has
no combinatorial Tits unit form: the builder returns ``{"error": <loud message>}``
(the Plan-30 honest-per-entry precedent), never a 500 and never a silent default.
"""
from quiverlab.errors import QuiverlabError

_REFS = ["bongartz_criterion", "bdps_tame_tits", "kasjan_skowronski",
         "ovsienko_forms", "vonhohne_wnn", "bjp_quadratic_forms"]


def _form_latex(gram, labels) -> str:
    """q_A(x) written out from the integer Gram matrix (diagonal 2): the squares plus
    the off-diagonal cross terms G_ij x_i x_j (i < j)."""
    def sym(v):
        return "x_{%s}" % v
    n = len(gram)
    terms = [sym(labels[i]) + "^{2}" for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            c = gram[i][j]
            if c == 0:
                continue
            mag = abs(c)
            coeff = "" if mag == 1 else str(mag)
            terms.append(("-" if c < 0 else "+", coeff + sym(labels[i]) + sym(labels[j])))
    out = terms[0]
    for t in terms[1:]:
        if isinstance(t, tuple):
            out += " " + t[0] + " " + t[1]
        else:
            out += " + " + t
    return "q_A(x) = " + out


def tame_wild_block(A) -> dict:
    """The tame_wild block: the combinatorial Tits Gram matrix, the minimal-relation
    counts r_ij, the two form booleans (weakly positive / weakly nonnegative) with the
    exact witness, the P56 simple/strong-simple-connectivity trail, and the
    theorem-gated rep-finite/tame/wild verdict (None off scope) with the honest
    scope note. Loud inputs (presentation-less / non-triangular) -> {"error": ...}."""
    try:
        from quiverlab.invariants.tits import tame_wild_certificate
        c = tame_wild_certificate(A)
    except QuiverlabError as exc:
        return {"kind": "tame_wild", "error": str(exc), "references": list(_REFS)}
    labels = list(A.quiver.vertices)
    mrc = {f"{s},{t}": v for (s, t), v in c.minimal_relation_counts.items()}
    strong = c.strong_certificate
    strong_blk = None
    if strong is not None:
        w = strong.witness
        strong_blk = {
            "verdict": strong.verdict,
            "reason": strong.reason,
            "checked_convex": strong.checked_convex,
        }
    return {
        "kind": "tame_wild",
        "gram": [list(row) for row in c.gram],
        "minimal_relation_counts": mrc,
        "is_unit_form": c.is_unit_form,
        "weakly_positive": c.weakly_positive,
        "weakly_nonnegative": c.weakly_nonnegative,
        "witness": list(c.witness) if c.witness is not None else None,
        "witness_value": c.witness_value,
        "simply_connected": c.simply_connected,
        "strongly_simply_connected": c.strongly_simply_connected,
        "strong_certificate": strong_blk,
        "field_alg_closed": c.field_alg_closed,
        "rep_type": c.rep_type,
        "certified": c.certified,
        "reason": c.reason,
        "scope_note": c.scope_note,
        "labels": [str(v) for v in labels],
        "latex": _form_latex(c.gram, labels),
        "references": list(_REFS),
    }
