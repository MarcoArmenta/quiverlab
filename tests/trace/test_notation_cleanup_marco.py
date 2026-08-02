"""Marco 2026-07-31 notation-cleanup wave (render/prose only): the report now

  * writes classes directly over the ordered basis -- NO ``= e_k`` coordinate inline
    (the coordinate vectors stay in the JSON, said once per section);
  * opens each Hochschild (co)homology section with a TYPING statement (what the
    engine computes + what the bar-bracket / tensor notation means), bar vs CS;
  * uses UNIFORM product notation alpha^n_j / z^n_j (no beta / gamma / w);
  * caps displays at 50 -- a matrix with >=50 rows/cols and a basis listing past 50
    state the size and point at the JSON;
  * carries a fine-grained, nested table of contents (per section, per degree);
  * states a ZERO (co)homology degree in ONE line and keeps its anchor (never a
    silent omission, never an empty enumeration/classes scaffold);
  * explains the Ext/Tor class-label syntax (P_v#k, n_{v,j}).

These are contract tests over the RENDERED HTML; the block data is unchanged.
"""
import pathlib

import pytest

import quiverlab as ql
from quiverlab.trace.render_html import (
    render_html, matrix_grid, product_degree_sections, module_reps_sections,
    hh_typing_html)
from quiverlab.trace.results_html import results_section

ROOT = pathlib.Path(__file__).resolve().parents[2]


# --------------------------------------------------------------------------- #
# (1) No coordinate-vector inline; classes written over the ordered basis.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_no_coordinate_vector_inline_in_classes():
    block = {"kind": "HH^", "dims": [2, 2], "engine": "bar/GF(2)",
             "basis_classes": {"1": [
                 {"terms": [["1", ["x"], "x"]], "vector": [[1, "1"]],
                  "kind": "cochain"}]},
             "chain_basis": {"1": ["[x -> e_1]", "[x -> x]"]},
             "differentials": {}}
    html = "".join(results_section({"hh_cohomology": block}))
    assert "[x → x]" in html
    # the class line never carries a "= e_1"/"= e_2" coordinate tail
    assert "= e_1</li>" not in html and "= e_2</li>" not in html
    assert "written over the ordered basis" in html


# --------------------------------------------------------------------------- #
# (2) Typing statement at the top of each (co)homology section, bar vs CS.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_typing_statement_bar_route():
    html = hh_typing_html("hh_cohomology", "bar")
    assert "What the engine computes" in html
    assert "Hom_k(Ā" in html and "Hom_{A^e}" in html
    assert "tensor products over k" in html and "composition of arrows" in html


@pytest.mark.oracle_selfcert
def test_typing_statement_cs_route_and_homology():
    coh = hh_typing_html("HH^", "cs")
    assert "Chouhy–Solotar" in coh and "e_{o(σ)} A e_{t(σ)}" in coh
    hom = hh_typing_html("HH_", "bar")
    assert "A ⊗_k Ā" in hom


@pytest.mark.oracle_selfcert
def test_hh_block_opens_with_typing_statement():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    block = {"kind": "HH^", "dims": [2, 1, 1], "engine": "bar/GF(7)"}
    html = "".join(results_section({"hh_cohomology": block}))
    assert "What the engine computes" in html
    # the typing precedes the dims table
    assert html.index("What the engine computes") < html.index("dim HH^n")


# --------------------------------------------------------------------------- #
# (3) Uniform product notation: no beta / gamma / w in the products source.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_products_source_has_no_beta_gamma_w():
    src = (ROOT / "src" / "quiverlab" / "trace" / "products.py").read_text(
        encoding="utf-8")
    # the symbol dicts and legend use only alpha / z now
    assert r"\beta" not in src and r"\gamma" not in src
    # the legend text uses alpha^n_j / z^n_j
    assert "α^n_j" in src and "z^n_j" in src


@pytest.mark.oracle_selfcert
def test_cup_equation_uses_alpha_throughout():
    A = ql.truncated_polynomial(2, field=ql.GF(2))
    from quiverlab.trace.products import products_chapter
    from tests.trace._matrix_grid import cayley_headers
    html = render_html(products_chapter(A, "cup", A.cup_products(1)),
                       title="cup", algebra=A)
    # Marco 2026-08-01: the Cayley grid's headers are the alpha classes and the corner
    # is the cup operator; the cells (in the alpha output basis) carry no beta/gamma.
    cols, rows = cayley_headers(html)
    assert r"\alpha^{0}_{1}" in cols and r"\alpha^{0}_{1}" in rows
    assert r"\cup" in html
    assert r"\beta" not in html and r"\gamma" not in html


# --------------------------------------------------------------------------- #
# (4) Display caps: 50x50 matrices, 50-element basis listings.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_matrix_over_50_columns_elides_to_json_pointer():
    big = [[str(j) for j in range(51)] for _ in range(2)]
    html = matrix_grid(big, label="M")
    assert "ql-matrix" not in html                     # NOT rendered as a grid
    assert "display cap" in html and "JSON record" in html
    assert "2×51 matrix" in html


@pytest.mark.oracle_selfcert
def test_matrix_under_50_still_renders_grid():
    small = [["1", "0"], ["0", "1"]]
    assert "ql-matrix" in matrix_grid(small)


@pytest.mark.oracle_selfcert
def test_basis_listing_capped_at_50_with_pointer():
    enum = ["u_%d" % i for i in range(60)]
    block_reps = {"0": [{"terms": [["1", "P_1", "n_1,1"]], "vector": [[0, "1"]]}]}
    html = "".join(module_reps_sections(block_reps, {"0": enum}, None, "ext", "cr"))
    assert "u_0" in html and "u_49" in html            # first 50 shown
    assert "u_59" not in html                          # the 60th is not
    assert "10 more" in html and "machine record" in html


# --------------------------------------------------------------------------- #
# (5) Fine-grained nested table of contents.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_toc_is_nested_with_degree_subsections():
    # Marco 2026-08-02: products no longer carry per-degree sub-sections, so the nested
    # ToC is exercised on a plain HH block (which KEEPS its per-degree structure).
    from quiverlab.hochschild.hh_reps import hh_reps_blocks
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    tbl = A.hochschild_cohomology(2, verbose=False)
    block = {"kind": tbl.kind, "top": 2, "dims": list(tbl.dims), "engine": tbl.engine}
    block.update(hh_reps_blocks(A, "hh_cohomology", 2, list(tbl.dims), tbl.engine) or {})
    html = render_html([], title="hh", algebra=A, results={"hh_cohomology": block})
    assert "<h2>Contents</h2>" in html
    assert "ql-toc-sub" in html                        # nested list present
    # the Computed-results section lists its per-degree subsections, linked
    assert "href='#cr-hh_cohomology-hh-coh-deg-1'" in html
    assert "id='cr-hh_cohomology-hh-coh-deg-1'" in html   # the target exists


# --------------------------------------------------------------------------- #
# (6) Zero degree: one-line statement, anchor kept, no empty scaffold.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_zero_cohomology_degree_states_one_line_keeps_anchor():
    block = {"kind": "HH^", "dims": [1, 0, 1], "engine": "bar/GF(7)",
             "basis_classes": {
                 "0": [{"terms": [["1", [], "e_1"]], "vector": [[0, "1"]],
                        "kind": "cochain"}],
                 "1": [],
                 "2": [{"terms": [["1", [], "e_1"]], "vector": [[0, "1"]],
                        "kind": "cochain"}]},
             "chain_basis": {"0": ["e_1"], "2": ["e_1"]}, "differentials": {}}
    html = "".join(results_section({"hh_cohomology": block}))
    # the zero degree keeps its anchor and states the vanishing in one line
    assert "id='cr-hh_cohomology-hh-coh-deg-1'" in html
    assert "HH^{1} = 0" in html
    # ...with NO enumeration / classes scaffold for that degree
    assert "degree-1 cohomology space" not in html
    # the nonzero degrees keep the full structure
    assert "degree-0 cohomology space" in html
    assert "degree-2 cohomology space" in html


@pytest.mark.oracle_selfcert
def test_zero_ext_degree_states_one_line():
    reps = {"0": [{"terms": [["1", "P_1", "n_1,1"]], "vector": [[0, "1"]]}],
            "1": []}
    html = "".join(module_reps_sections(reps, {"0": ["[P_1 -> n_1,1]"]}, None,
                                        "ext", "cr"))
    assert "id='cr-ext-deg-1'" in html
    assert "Ext^{1} = 0" in html


# --------------------------------------------------------------------------- #
# (7) Product sections drop the annihilating differential; HH sections keep it.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_product_section_drops_differential_hh_section_keeps_it():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    reps = A.cup_products(1).blocks().get("basis_classes")
    cb = A.cup_products(1).blocks().get("chain_basis")
    diffs = A.cup_products(1).blocks().get("differentials")
    prod = "".join(product_degree_sections(reps, cb, diffs, "p",
                                           show_differential=False))
    hh = "".join(product_degree_sections(reps, cb, diffs, "h",
                                         show_differential=True))
    assert "Verification" not in prod                  # product: no differential
    assert "Verification" in hh                         # HH: differential kept
