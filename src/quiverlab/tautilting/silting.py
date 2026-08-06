"""2-term silting bridge (Plan 45 / C4). A support tau-tilting pair (M, P) IS a 2-term
silting object ``T(M,P) = M (+) P[1]`` in K^b(proj A): each module summand M_i presents as
``[P_1(M_i) -> P_0(M_i)]`` and each killed projective P_v as ``P_v[1] = [P_v -> 0]`` (AIR
2014 Thm 3.2). The count of distinct 2-term silting complexes is the fourth leg of the AIR
four-way identity ``#s-tau-tilt = #f.f. torsion = #2-term silting = #semibricks``.

Soft P43 dependency: the raw ``(P1_vertices, P0_vertices, d1)`` tuples are emitted ALWAYS
(the payload the report/GUI render); the :class:`ChainComplex` wrapper is built only when
``quiverlab.modules.complexes`` imports (a convenience for the derived-category surface)."""
from __future__ import annotations

from quiverlab.modules.resolution import minimal_resolution

try:                                                  # SOFT P43 dependency
    from quiverlab.modules.complexes import ChainComplex
except ImportError:                                   # pragma: no cover
    ChainComplex = None


def two_term_silting(pair):
    """The 2-term silting object of ``pair`` as ``{"summands": [{"P1": [v...], "P0": [v...],
    "d1": [[...]]}, ...], "complex": <ChainComplex | None>}``. Module summands present as
    ``P_1 -> P_0``; killed projectives as the shift ``P_v[1]``. The ChainComplex is built
    iff ``modules.complexes`` imports."""
    A = pair.algebra
    verts = list(A.quiver.vertices)
    summands = []
    for Mi in pair.summands:
        terms, dmats = minimal_resolution(Mi, 1)
        p0 = list(terms[0].vertices)
        p1 = list(terms[1].vertices) if len(terms) > 1 else []
        d1 = dmats[1] if (len(dmats) > 1 and p1) else []
        summands.append({"P1": p1, "P0": p0, "d1": d1})
    for v in sorted(pair.support, key=verts.index):
        summands.append({"P1": [v], "P0": [], "d1": []})     # P_v[1]
    complex_obj = _build_complex(pair) if ChainComplex is not None else None
    return {"summands": summands, "complex": complex_obj}


def _build_complex(pair):
    """The whole 2-term silting object as a :class:`ChainComplex` (degree 1 = the shifted /
    P_1 part, degree 0 = the P_0 part). Returns None if the build cannot be certified (the
    wrapper is a soft convenience -- never fatal)."""
    from quiverlab.tautilting import _twoterm as tt
    try:
        comps = tt.to_complexes(pair)
        if not comps:
            return None
        C, _ = tt.direct_sum_complexes(comps)
        P0 = tt._proj_sum(C.algebra, C.terms.get(0, []))
        P1 = tt._proj_sum(C.algebra, C.terms.get(-1, []))
        terms = {}
        if P0.dim:
            terms[0] = P0
        if P1.dim:
            terms[1] = P1
        if not terms:
            return None
        dmats = {}
        d = C.diffs.get(-1)
        if d and d[0] and P1.dim and P0.dim:
            dmats[1] = d
        return ChainComplex(terms, dmats, check=True)
    except Exception:
        return None


def silting_count(A, budget=512):
    """The number of distinct 2-term silting complexes over ``A`` (the fourth leg of the
    AIR four-way identity), via the exchange graph. By the AIR bijection this equals the
    number of support tau-tilting pairs; honest complete-iff contract."""
    from quiverlab.tautilting.mutation import exchange_graph
    eg = exchange_graph(A, budget_pairs=budget)
    if not eg.is_complete:
        return {"count": 0, "complete": False, "status": eg.status}
    return {"count": len(eg.vertices), "complete": True, "status": "complete"}
