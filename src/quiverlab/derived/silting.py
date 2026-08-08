"""Silting objects in K^b(proj A): verifier, single mutation, bounded exploration,
and the co-t-structure dictionary (Plan 67 / Aihara-Iyama arXiv:1009.3370, J. LMS 85
(2012) 633-668).

A silting object T satisfies Hom_{D^b}(T, T[n]) = 0 for all n > 0 (PRESILTING) and
thick(T) = K^b(proj A) (GENERATION). This is WEAKER than tilting (which also needs the
n < 0 vanishing); tilting => silting (AI Def 2.1). Presilting is DECIDED on the exact
positive window (perfect => bounded => finite window; outside it hyper-Hom is provably
the zero cochain group). Generation is three-valued: certified True on the tilting /
2-term / local classes (where a completion theorem applies), else 'unknown' -- det(g_proj)
= +-1 is NECESSARY (AI Thm 2.27) but NOT sufficient in general (thick subcategories are
not K0-classified; Krah phantom, arXiv:2302.12502). The g-matrix is computed in the
PROJECTIVE basis K0(K^b proj A) = (+)_v Z[P_v] via ``derived.tilting.g_proj`` (Plan 67
Task 0), NOT the composition-factor basis -- ``_chi`` gives det(Cartan . g) and is wrong
on non-unimodular Cartan (self-injective/symmetric; Example 2.47).

Single mutation (AI Def 2.30/2.34) is one approximation triangle: at an indecomposable
summand X, the cone of the minimal left add(T/X)-approximation (left mutation mu^+) or the
cocone of the right one (right mutation mu^-). Every mutant is RE-VERIFIED silting, shares
exactly n-1 summands with its input, and mu^- o mu^+ = id (AI Thm 2.31 / Prop 2.33). The
bounded-radius exploration truncates loudly (the silting quiver can be infinite -- kA2
already is; transitivity is proven only for local/hereditary/canonical, AI Thm 1.2, so
there is NO general BFS/enumeration claim). Float-free; exact linear algebra only."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.modules.complexes import hyper_hom_dims
from quiverlab.derived.tilting import (is_tilting_complex, _direct_sum_complex, _span,
                                       g_proj)


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


def _is_two_term(summands):
    """Every summand concentrated in two consecutive degrees, all up to ONE common shift
    (2-term silting is defined up to shift). True for the empty/stalk cases."""
    spans = [_span(T) for T in summands if T.degrees()]
    if not spans:
        return True
    widths = {hi - lo for lo, hi in spans}
    if not (widths <= {0, 1}):
        return False
    tops = {hi for _, hi in spans}
    return max(tops) - min(tops) <= 1          # a common 2-window covers them


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
    if is_tilting_complex(summands).is_tilting:     # is_tilting_complex now g_proj-routed
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
