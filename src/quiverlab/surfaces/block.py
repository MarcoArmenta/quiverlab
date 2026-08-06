"""surface_block: the algebra-adjacent descriptor of a triangulated surface (Plan 48).
Packages the surface invariants + the resulting gentle algebra's presentation, gentle
verdict, and AG invariant (P46) for the report and preset provenance. Not a webapp compute
kind -- surfaces are an INPUT method; the produced algebra runs through every existing
compute kind. Float-free."""
from __future__ import annotations

from quiverlab.invariants.recognizers import is_gentle
from quiverlab.surfaces.qp import jacobian_of, quiver_of

_REFS = ("fomin_shapiro_thurston", "labardini", "abcp")


def surface_block(T, field=None) -> dict:
    """A descriptor dict of the triangulated surface ``T``: the surface invariants, the
    resulting Jacobian's quiver/relations + dimension, the P46 gentle verdict, and (when
    gentle) the AG derived invariant. ``kind == "surface"``."""
    S = T.surface
    Q = quiver_of(T)
    J = jacobian_of(T, field=field)
    gentle = is_gentle(J)
    block = {
        "kind": "surface",
        "genus": S.genus,
        "boundary_marked": list(S.boundary_marked),
        "punctures": S.punctures,
        "arc_count": S.arc_count(),
        "triangle_count": S.triangle_count(),
        "euler_characteristic": S.euler_characteristic(),
        "vertices": list(Q.vertices),
        "arrows": {a: list(Q.arrows[a]) for a in Q.arrows},
        "relations": [str(r) for r in (J.relations or ())],
        "dim": J.dim,
        "is_gentle": gentle,
        "references": list(_REFS),
        "citations": list(_REFS),
    }
    if gentle:
        from quiverlab.strings.ag import ag_invariant
        block["ag_invariant"] = [list(pair) for pair in ag_invariant(J).pairs]
    return block
