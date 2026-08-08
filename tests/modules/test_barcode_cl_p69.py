"""Commutative-ladder generalized persistence diagrams (Plan 69 / R33). Self-cert: the
diagram's summand dim-vectors sum to M; every summand matched to an AR-quiver vertex.
Cross-engine: the diagram's distinct-indecomposable count is bounded by the knit's
indecomposable count. Honest scope: CL(n<=4) only (char 0 / char > dim); CL(n>=5) refused
at construction; a specific CL(3) worked example is # PIN vs the Escolar-Hiraoka figure."""
import pytest
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ, GF
from quiverlab.families.commutative_ladder import CommutativeLadder
from quiverlab.modules.barcode import barcode

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@selfcert
def test_cl3_diagram_sums_to_module_and_indexes_ar():
    A = CommutativeLadder(3, field=QQ)
    M = A.projective(A.quiver.vertices[0])          # a concrete CL(3) module (P of a source)
    bc = barcode(M)
    assert bc.kind == "commutative_ladder" and bc.ar_status == "complete"
    acc = {}
    for entry in bc.diagram:                        # (ar_index, ar_name, dimvec, mult, is_interval)
        assert entry.ar_index is not None          # matched to an AR vertex
        for v, d in entry.dimvec.items():
            acc[v] = acc.get(v, 0) + entry.multiplicity * d
    assert acc == {v: d for v, d in M.dimension_vector().items() if d}


@xeng
def test_cl3_diagram_indecs_are_knit_vertices():
    from quiverlab.modules.ar import knit_ar_quiver
    A = CommutativeLadder(3, field=QQ)
    ar = knit_ar_quiver(A, budget_modules=400, budget_dim=8000)
    M = A.injective(A.quiver.vertices[-1])
    bc = barcode(M)
    assert len({e.ar_index for e in bc.diagram}) <= len(ar.vertices)   # subset of the 29


def test_cl_char_scope_refused_off_char0():
    # CL indecomposables need not be bricks -> AR knit / decompose char-scoped.
    A = CommutativeLadder(3, field=GF(2))            # char 2 <= dim -> knit/decompose refuse
    with pytest.raises(QuiverlabError):
        barcode(A.projective(A.quiver.vertices[0]))


def test_cl_char_scope_knit_surfaces_error_not_silent_complete():
    # M6: verify the refusal CHAIN, not just that barcode raises. Over GF(2) the knit
    # CATCHES the internal decompose char-refusal and returns status="error",
    # is_complete=False -- it must NEVER silently drain to a truncated "complete".
    # (LIVE-VERIFIED 2026-08-08: CL(3)/GF(2) -> status="error", is_complete=False, 8 verts.)
    from quiverlab.modules.ar import knit_ar_quiver
    A = CommutativeLadder(3, field=GF(2))
    ar = knit_ar_quiver(A, budget_modules=400, budget_dim=8000)
    assert ar.status == "error" and ar.is_complete is False
