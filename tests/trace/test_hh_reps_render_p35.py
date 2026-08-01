"""Plan 35 wave 3d -- the PLAIN hh_cohomology / hh_homology blocks now RENDER the
element-wise classical dictionary (central elements / derivations / deformation cochain
/ commutator residues) AND the per-degree explicit representatives, in the report
(Computed-results HTML) and both GUIs. Rendering-only: every number is the capture
layer's (hochschild.hh_reps); the read-offs are the ONE shared builder
``trace.interpretations`` both surfaces read. Element-wise ONLY where reps are present;
the framing sentence covers every other degree.
"""
import pathlib

import pytest

import quiverlab as ql
from quiverlab.fields.primefield import PrimeField
from quiverlab.hochschild.hh_reps import hh_reps_blocks
from quiverlab.trace.results_html import _block_html


def _hh_block(a, p, kind, top):
    A = ql.truncated_polynomial(a, field=PrimeField(p))
    tbl = (A.hochschild_cohomology if kind == "hh_cohomology"
           else A.hochschild_homology)(top, verbose=False)
    block = {"kind": tbl.kind, "top": top, "dims": list(tbl.dims), "engine": tbl.engine}
    block.update(hh_reps_blocks(A, kind, top, list(tbl.dims), tbl.engine) or {})
    return block


# --------------------------------------------------------------------------- #
# (1) HH^1 renders D(arrow)=value + the inner-derivation subspace dimension.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_HH1_renders_derivation_and_inner_dim():
    html = "".join(_block_html("hh_cohomology", _hh_block(2, 7, "hh_cohomology", 3)))
    assert "D(x) = x" in html
    assert "OUTER DERIVATIONS" in html                    # framing still present
    assert "inner derivations" in html and "rank δ⁰" in html
    # the element-wise read-off is under the explicit-class name alpha^1_1
    assert "\\alpha^{1}_{1}" in html


@pytest.mark.oracle_selfcert
def test_HH0_renders_central_elements():
    html = "".join(_block_html("hh_cohomology", _hh_block(2, 7, "hh_cohomology", 2)))
    assert "central element" in html
    # both central elements of the dual numbers are named
    assert "e_1" in html and "x" in html


@pytest.mark.oracle_selfcert
def test_HH2_renders_deformation_heading():
    html = "".join(_block_html("hh_cohomology", _hh_block(3, 5, "hh_cohomology", 3)))
    assert "2-cocycle μ(a, b)" in html


@pytest.mark.oracle_selfcert
def test_HH0_homology_renders_commutator_residues():
    html = "".join(_block_html("hh_homology", _hh_block(2, 7, "hh_homology", 2)))
    assert "residue of an element modulo the commutators" in html


# --------------------------------------------------------------------------- #
# (2) The per-degree explicit-reps sections render with stable anchors.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_per_degree_reps_sections_and_anchors():
    html = "".join(_block_html("hh_cohomology", _hh_block(2, 7, "hh_cohomology", 3)))
    assert "Explicit representatives by degree" in html
    # product_degree_sections anchors, namespaced by the block kind + side
    assert "id='cr-hh_cohomology-hh-coh-deg-1'" in html
    assert "Basis classes" in html and "Verification" in html
    # the term-sum -> coordinate-vector view of the HH^1 class
    assert "[x → x]" in html


@pytest.mark.oracle_selfcert
def test_homology_reps_sections_render():
    html = "".join(_block_html("hh_homology", _hh_block(2, 7, "hh_homology", 3)))
    assert "id='cr-hh_homology-hh-hom-deg-0'" in html
    assert "every 0-chain is a cycle" in html             # b_0 = 0 verification


# --------------------------------------------------------------------------- #
# (3) Tolerance: an old-cache block with no reps renders just the dims + framing.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_tolerates_block_without_reps():
    html = "".join(_block_html("hh_cohomology",
                               {"kind": "HH^", "dims": [2, 1, 1], "engine": "x"}))
    assert "classical dictionary" in html                 # framing still renders
    assert "Explicit representatives by degree" not in html
    assert "rank δ⁰" not in html                           # the element-wise read-off is absent


# --------------------------------------------------------------------------- #
# (4) Both gui.js copies mirror the HH element read-offs + reps; byte-identical.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_gui_js_wires_hh_reps():
    root = pathlib.Path(__file__).resolve().parents[2]
    a = (root / "docs/gui/gui.js").read_text(encoding="utf-8")
    b = (root / "webapp/static/gui/gui.js").read_text(encoding="utf-8")
    assert a == b, "gui.js copies must be byte-identical"
    for needle in ("appendHHInterpretation", "appendHHReps", "derivationValues",
                   "deformationCochain", "HH_INTERP", "rank δ⁰"):
        assert needle in a
