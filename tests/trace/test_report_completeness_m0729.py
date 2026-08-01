"""The report is the session's complete, keepable record (Marco, 2026-07-29).

His two artifacts -- the draw page after Compute (example-a) and the downloaded
worked-steps HTML (example-b) -- disagreed: the page showed rad/top/soc, the AR
translates, the resolutions and the Ext/Tor tables; the report showed only the
worked steps of the one traced computation, so closing the tab lost everything
else. Plus four presentation defects, all pinned here:

  * matrices behind a horizontal SCROLLBAR (rad/top/soc, projectives/injectives);
  * an arrow acting as the exact ZERO map printed as a zero block (soc M, arrow d);
  * a Result section that was one long equation instead of a degree table;
  * differentials repeated verbatim across degrees in a periodic resolution.

And one honesty fix: the resolution terms are NAMED as projective bimodules
instead of quoted as a bare generator count.
"""
import pytest

from quiverlab import CC, GF, Quiver
from quiverlab.trace.events import ResolutionTerm
from quiverlab.trace.recorder import Trace
from quiverlab.trace.render_html import _fit_pct, render_html
from quiverlab.trace.results_html import normalize, results_section

from tests.trace._matrix_grid import grid_indices, grids, has_grid
from tests.trace._result_table import result_dims

pytestmark = [pytest.mark.oracle_selfcert]


# A 1-dimensional-per-vertex module over that algebra: a and c act by 1; the loops
# b and d MUST act as zero (b^2 = d^2 = 0 forces a vanishing scalar in dimension 1).
# Basis order is by vertex, so index 0 <-> v1, 1 <-> v2, 2 <-> v3.
_Z = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
_MAPS = {"a": [[0, 1, 0], [0, 0, 0], [0, 0, 0]],
         "b": _Z,
         "c": [[0, 0, 0], [0, 0, 0], [0, 1, 0]],
         "d": _Z}


# Marco's example-a algebra: three vertices, a loop at 2 and at 3, b^2 = d^2 = 0.
def _marco_algebra(field=CC):
    Q = Quiver(vertices=[1, 2, 3],
               arrows={"a": (2, 1), "b": (2, 2), "c": (2, 3), "d": (3, 3)})
    return Q.algebra(relations=["b*b", "d*d"], field=field)


# --------------------------------------------------------------------------- #
# 1. No scrollbars: nothing in the page clips, wide matrices shrink instead.
# --------------------------------------------------------------------------- #
def test_no_overflow_clipping_anywhere_in_the_page():
    A = _marco_algebra()
    tr = Trace()
    A.hochschild_homology(3, engine="cs", trace=tr)
    html = render_html(list(tr), title="t", algebra=A)
    assert "overflow" not in html            # no scroll box, no clip, no ellipsis


def test_wide_matrices_shrink_to_fit_and_never_magnify():
    assert _fit_pct(None) is None            # unknown width -> full size
    assert _fit_pct(4) is None               # comfortably fits
    assert _fit_pct(16) is None              # the boundary still fits
    assert _fit_pct(20) == 80                # past it: proportional shrink
    assert _fit_pct(1000) == 50              # ...floored, never vanishing
    assert all(_fit_pct(c) <= 100 for c in range(17, 200))   # shrink-ONLY


# --------------------------------------------------------------------------- #
# 2. An arrow acting as zero is named, not printed.
# --------------------------------------------------------------------------- #
def test_zero_arrow_is_named_not_printed():
    block = {"kind": "rad_top_soc",
             "radical": {"dims": {"1": 0, "2": 0, "3": 2},
                         "maps": {"d": [[0, 0], [0, 0]]}},
             "top": {"dims": {"1": 1}, "maps": {}},
             "socle": {"dims": {"3": 2}, "maps": {"d": [["0", "0"], ["0", "0"]],
                                                  "c": [["1", "0"], ["0", "1"]]}}}
    html = "".join(results_section({"rad_top_soc": block}))
    # the all-zero views say so once; the zero arrow of a MIXED view is named...
    assert "every arrow acts as zero" in html
    assert "arrow d acts as zero" in html
    # ...and its zero block is nowhere to be seen, while the live arrow's is.
    assert grids(html) == [[["1", "0"], ["0", "1"]]]   # exactly one: arrow c's


def test_a_zero_entry_inside_a_live_matrix_is_still_printed():
    """The filter is per-ARROW, not per-entry: a matrix with some zeros is live."""
    block = {"kind": "rad_top_soc",
             "radical": {"dims": {"1": 2}, "maps": {"a": [["0", "1"], ["0", "0"]]}},
             "top": {"dims": {}, "maps": {}},
             "socle": {"dims": {}, "maps": {}}}
    html = "".join(results_section({"rad_top_soc": block}))
    assert "rad M, arrow a" in html
    assert has_grid(html, [["0", "1"], ["0", "0"]])


# --------------------------------------------------------------------------- #
# 3. The Result is a degree table.
# --------------------------------------------------------------------------- #
def test_result_is_a_degree_table():
    A = _marco_algebra()
    tr = Trace()
    table = A.hochschild_homology(4, engine="cs", trace=tr)
    from quiverlab.trace.writer import _authoritative_result
    html = render_html(_authoritative_result(list(tr), table), title="t", algebra=A)
    assert 'class="ql-dims"' in html
    assert result_dims(html, "dim HH_n") == list(table.dims)


# --------------------------------------------------------------------------- #
# 4. A repeated differential is referenced, not reprinted.
# --------------------------------------------------------------------------- #
def test_repeated_differentials_reference_the_first_occurrence():
    """Marco's example-b: b_2 = b_4 = b_6 = ... in the periodic resolution, each
    printed in full. Now the repeats say 'the same matrix as above'."""
    A = _marco_algebra()
    tr = Trace()
    A.hochschild_homology(8, engine="cs", trace=tr)
    html = render_html(list(tr), title="t", algebra=A)
    assert "the same matrix as above" in html
    # the two distinct matrices are each shown once; the 7 repeats are references
    assert html.count("the same matrix as above") >= 5


def test_repeated_resolution_differentials_are_referenced():
    same = [["1", "0"], ["0", "1"]]
    block = {"kind": "projective_resolution", "terms": [{}, {}, {}],
             "summands": ["P_{1}", "P_{1}", "P_{1}"], "pd": None,
             "differentials": [{"rows": 2, "cols": 2, "matrix": same},
                               {"rows": 2, "cols": 2, "matrix": [["0", "1"], ["0", "0"]]},
                               {"rows": 2, "cols": 2, "matrix": same}]}
    html = "".join(results_section({"projective_resolution": block}))
    assert "the same matrix as above" in html
    assert html.count("the same matrix as above") == 1     # only d_2 repeats d_0


def test_an_elided_differential_is_never_matched_as_a_repeat():
    """Two elided differentials have no bodies, so claiming they are equal would be
    a fabrication -- each is stated by shape instead."""
    block = {"kind": "projective_resolution", "terms": [{}, {}],
             "summands": ["P_{1}", "P_{1}"], "pd": None,
             "differentials": [{"rows": 900, "cols": 900, "elided": True},
                               {"rows": 900, "cols": 900, "elided": True}]}
    html = "".join(results_section({"projective_resolution": block}))
    assert "the same matrix as above" not in html
    assert html.count("body not recorded") == 2


# --------------------------------------------------------------------------- #
# 5. The resolution terms are NAMED as projective bimodules.
# --------------------------------------------------------------------------- #
def test_resolution_terms_name_their_projective_summands():
    A = _marco_algebra()
    tr = Trace()
    A.hochschild_homology(2, engine="cs", trace=tr)
    terms = [e for e in tr if isinstance(e, ResolutionTerm)]
    assert terms and all(t.corners for t in terms)
    # degree 0: one P(v,v) per vertex; degree 1: one P(s,t) per arrow.
    assert sorted(terms[0].corners) == [(1, 1), (2, 2), (3, 3)]
    assert sorted(terms[1].corners) == [(2, 1), (2, 2), (2, 3), (3, 3)]
    html = render_html(list(tr), title="t", algebra=A)
    assert "P(1,1)" in html and "P(2,1)" in html
    assert r"A e_{v} \otimes e_{w} A" in html          # the notation is defined


def test_a_term_without_recorded_corners_claims_nothing():
    """The bar resolution over a structure-constants algebra is not vertex-graded,
    so no decomposition is recorded -- and none is invented."""
    ev = [ResolutionTerm(degree=0, n_generators=3, collapsed_dim=3)]
    html = render_html(ev, title="t")
    assert "P(" not in html


# --------------------------------------------------------------------------- #
# 6. The report saves EVERY computed block.
# --------------------------------------------------------------------------- #
def test_normalize_accepts_both_runner_shapes():
    server = {"cartan": {"latex": "X"}, "center": {"dim": 1}}
    assert normalize(server) == [("cartan", {"latex": "X"}),
                                 ("center", {"dim": 1})]
    pyodide = [{"invariant": "ext:0..4", "dims": [1]},
               {"invariant": "cartan", "latex": "X"}]
    assert [k for k, _ in normalize(pyodide)] == ["ext", "cartan"]
    assert normalize(None) == [] and normalize("nonsense") == []


def test_every_computed_block_reaches_the_report():
    blocks = {
        "hh_cohomology": {"kind": "HH^", "dims": [5, 2, 2], "engine": "bar"},
        "cartan": {"matrix": [[1, 2], [0, 1]],
                   "latex": r"\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}"},
        "coxeter_polynomial": {"latex": "t^{2} + t + 1"},
        "global_dimension": {"text": "gl.dim = 1 (exact)"},
        "center": {"dim": 3},
        "dimension": {"value": 15},
        "dimension_vector": {"latex": r"\underline{\dim}\, M = (3,\, 4,\, 3)"},
        "decompose": {"iso_classes": 1,
                      "summands": [{"multiplicity": 1, "dim_vector": {"1": 3}}]},
        "ext": {"dims": [1, 7, 3], "target": {"dimvec": {"1": 2}}},
        "tor": {"dims": [0, 6, 3], "target": {"dimvec": {"1": 2}}},
        "injective_dimension": {"value": None, "bound": 32,
                                "latex": r"\operatorname{id} M > 32",
                                "note": "certified lower bound"},
    }
    html = "".join(results_section(blocks))
    for expect in ("Hochschild cohomology", "Cartan matrix", "Coxeter polynomial",
                   "Global dimension", "Centre", "Dimension",
                   "Krull–Schmidt", "Ext", "Tor", "Injective dimension"):
        assert expect in html, expect
    assert result_dims(html, "dim HH^n") == [5, 2, 2]
    assert result_dims(html, "dim Ext^n") == [1, 7, 3]
    assert result_dims(html, "dim Tor_n") == [0, 6, 3]
    assert "certified lower bound" in html          # never a bare, unproven ∞
    assert "∞" not in html and r"\infty" not in html


def test_a_failed_computation_is_recorded_honestly():
    results = [{"invariant": "tau", "error": {"type": "RelationError",
                                              "message": "M violates b*b"}}]
    html = "".join(results_section(results))
    assert "not computed" in html and "M violates b*b" in html


def test_results_are_optional_and_absent_changes_nothing():
    A = _marco_algebra(field=GF(5))
    tr = Trace()
    A.hochschild_cohomology(2, engine="cs", trace=tr)
    ev = list(tr)
    assert render_html(ev, title="t", algebra=A) == \
        render_html(ev, title="t", algebra=A, results=None)
    assert "Computed results" not in render_html(ev, title="t", algebra=A)


# --------------------------------------------------------------------------- #
# 7. The modules themselves: maps and Loewy series, not just dimension vectors.
#    (Marco, second pass over example-b-NEW.)
# --------------------------------------------------------------------------- #
def test_the_modules_section_shows_maps_and_loewy_layers():
    A = _marco_algebra(field=GF(7))
    M = A.module({1: 1, 2: 1, 3: 1}, _MAPS, name="M")
    N = A.simple(2)
    html = render_html([], title="t", algebra=A, modules=[("M", M), ("N", N)])
    assert "<h2 id='modules'>The modules</h2>" in html
    # the dimension vector AND the Loewy series AND top/socle
    assert r"\underline{\dim} = (1{:}1,\ 2{:}1,\ 3{:}1)" in html
    assert r"\operatorname{top} M" in html and r"\operatorname{soc} M" in html
    assert r"\begin{matrix} S_{2} \end{matrix}" in html      # a Loewy layer
    # ...and the actual arrow matrices, which the dimension vector cannot give
    assert "M, arrow a:" in html and "M, arrow c:" in html
    assert "M: arrows b, d act as zero." in html             # zero arrows named once
    assert "<h3>N</h3>" in html


def test_modules_section_is_absent_and_harmless_without_modules():
    A = _marco_algebra(field=GF(7))
    assert "The modules" not in render_html([], title="t", algebra=A)
    # a module that cannot be described must not sink the report
    class _Broken:
        def __getattr__(self, name):
            raise RuntimeError("boom")
    html = render_html([], title="t", algebra=A, modules=[("M", _Broken())])
    assert "The modules" not in html and "<h1>" in html


def test_left_module_is_flagged_as_such():
    A = _marco_algebra(field=GF(7))
    html = render_html([], title="t", algebra=A,
                       modules=[("N", A.simple(2, side="left"))])
    assert "a LEFT module" in html


# --------------------------------------------------------------------------- #
# 8. Krull-Schmidt: standard summands NAMED, the rest shown in full.
# --------------------------------------------------------------------------- #
def test_standard_summands_are_named_and_need_no_matrices():
    from quiverlab.modules.qpa_module import summand_blocks
    Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)})
    A = Q.algebra(relations=[], field=GF(7))
    for mod, want in ((A.simple(1), ("simple", "1")),
                      (A.projective(2), ("projective", "2"))):
        blk = summand_blocks(mod)
        assert blk["standard"] == {"kind": want[0], "vertex": want[1]}
        assert "maps" not in blk               # the name IS the module
    html = "".join(results_section({"decompose": {
        "iso_classes": 2,
        "summands": [summand_blocks(A.projective(2)), summand_blocks(A.simple(1))]}}))
    assert "P_2" in html and "S_1" in html and "M_1" not in html
    assert grids(html) == []                   # no matrices for standard summands


def test_a_nonstandard_summand_is_shown_in_full():
    from quiverlab.modules.qpa_module import summand_blocks
    A = _marco_algebra(field=GF(7))
    M = A.module({1: 1, 2: 1, 3: 1}, _MAPS)
    blk = summand_blocks(M)
    assert "standard" not in blk and blk["maps"]          # unnamed -> full action
    html = "".join(results_section({"decompose": {"iso_classes": 1,
                                                  "summands": [blk]}}))
    assert "M_1, arrow a:" in html and grids(html)


def test_identify_standard_never_guesses():
    from quiverlab.modules.hom import identify_standard
    A = _marco_algebra(field=GF(7))
    M = A.module({1: 1, 2: 1, 3: 1}, _MAPS)
    assert identify_standard(M) is None            # not standard -> no name invented
    assert identify_standard(A.simple(3)) == ("simple", 3)


# --------------------------------------------------------------------------- #
# 9. No "Result" heading, and no HH table printed twice.
# --------------------------------------------------------------------------- #
def test_homology_is_named_and_never_printed_twice():
    A = _marco_algebra()
    tr = Trace()
    table = A.hochschild_homology(3, engine="cs", trace=tr)
    from quiverlab.trace.writer import _authoritative_result
    ev = _authoritative_result(list(tr), table)
    results = {"hh_homology": {"kind": "HH_", "dims": list(table.dims),
                               "engine": "cs"}}

    both = render_html(ev, title="t", algebra=A, results=results)
    assert "Hochschild homology" in both
    assert both.count('class="ql-dims"') == 1        # ONE table, not two
    assert ">Result<" not in both                    # never that empty heading

    alone = render_html(ev, title="t", algebra=A)    # no results passed
    assert "<h2 id='result'>Hochschild homology</h2>" in alone
    assert ">Result<" not in alone


def test_cohomology_gets_its_own_name():
    A = _marco_algebra()
    tr = Trace()
    table = A.hochschild_cohomology(2, engine="cs", trace=tr)
    from quiverlab.trace.writer import _authoritative_result
    html = render_html(_authoritative_result(list(tr), table), title="t", algebra=A)
    assert "<h2 id='result'>Hochschild cohomology</h2>" in html
    assert "Hochschild homology" not in html         # never the wrong variance


# --------------------------------------------------------------------------- #
# 10. Every matrix is an INDEXED GRID: an extra index row and column, and a rule
#     between cells (Marco 2026-07-29, third pass).
# --------------------------------------------------------------------------- #
def test_matrices_carry_row_and_column_indices():
    html = render_html([], title="t", results={"cartan": {
        "matrix": [[1, 2, 4], [0, 2, 4]],
        "latex": r"\begin{pmatrix} 1 & 2 & 4 \\ 0 & 2 & 4 \end{pmatrix}"}})
    assert has_grid(html, [[1, 2, 4], [0, 2, 4]])       # entries, verbatim
    cols, rows = grid_indices(html)
    assert cols == ["1", "2", "3"]                      # an extra header ROW
    assert rows == ["1", "2"]                           # an extra header COLUMN
    assert 'class="ql-matrix"' in html and "ql-corner" in html
    assert "border:1px solid #d0d0d0" in html           # the grey grid rule


def test_a_zero_dimensional_matrix_is_the_symbol_zero_not_an_empty_grid():
    from quiverlab.trace.render_html import matrix_grid
    assert "ql-matrix" not in matrix_grid([], label="d_{1}")
    assert "ql-matrix" not in matrix_grid([[]], label="d_{1}")
    assert "d_{1} = 0" in matrix_grid([], label="d_{1}")


# --------------------------------------------------------------------------- #
# 11. Plan 35: the HH product blocks (cup/cap/bracket/connes_b) render as tables
#     in the Computed-results section, not the generic "see the JSON" stub.
# --------------------------------------------------------------------------- #
def test_product_blocks_render_in_the_computed_results_section():
    """The DEFAULT report case: an HH computation claims the worked-steps bundle and
    the products ride along as result blocks. The block renderer used to have no
    cup/cap/bracket/connes_b branch, so the section printed a bare stub; now it
    renders the map labels + structure-constant equations (cup/cap/bracket) and the
    induced-B grids (connes_b)."""
    import quiverlab as ql
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    tr = Trace()
    table = A.hochschild_cohomology(2, engine="bar", trace=tr)
    results = {"hh_cohomology": {"kind": "HH^", "dims": list(table.dims),
                                 "engine": "bar"},
               "cup": A.cup_products(2).blocks(),
               "connes_b": A.connes_differentials(2).blocks()}
    html = render_html(list(tr), title="t", algebra=A, results=results)
    # the sections carry the gui.js i18n titles (was the bare "cup"/"connes_b"); the
    # h3 now carries a ToC anchor id (Marco 2026-08-01 fine-grained contents).
    assert "<h3 id='cr-cup'>Cup product tables</h3>" in html
    assert "<h3 id='cr-connes_b'>Connes differentials</h3>" in html
    assert ">cup</h3>" not in html and ">connes_b</h3>" not in html
    # ...the cup equations really render (the ∪ operator appears in the map label
    # AND the structure-constant equations), NOT the missing-branch stub...
    assert r"\cup" in html
    assert "computed; see the JSON record for its data" not in html
    # ...and the Connes differential is a matrix GRID with its induced-rank line.
    assert has_grid(html, A.connes_differentials(2).blocks()["matrices"]["0"])
    assert "induced rank B_0" in html


def test_product_tables_carry_a_notation_legend_with_the_basis():
    """Plan-35 follow-up (Marco): the Computed-results product tables are preceded by
    a legend defining the symbols (alpha/beta/gamma/z/w) and naming the concrete
    recorded basis, so the reader knows what the structure constants refer to."""
    import quiverlab as ql
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    block = A.cup_products(2).blocks()
    html = "".join(results_section({"cup": block}))
    assert "recorded basis" in html
    assert block["basis"] in html                  # e.g. "bar/GF(7)"
    # the legend appears BEFORE the first cup map label / equation.
    assert html.index("recorded basis") < html.index(r"\cup")


def test_a_fully_vanishing_product_bidegree_states_so_not_a_stub():
    block = {"kind": "cup", "engine": "bar/GF(7)",
             "tables": [{"degrees": [1, 1], "out_degree": 2, "dims": [1, 1, 1],
                         "constants": [[["0"]]]}]}
    html = "".join(results_section({"cup": block}))
    assert "every product in this bidegree vanishes" in html
    assert "computed; see the JSON record for its data" not in html


def test_grid_entries_are_copied_verbatim_and_escaped():
    from quiverlab.trace.render_html import matrix_grid
    got = matrix_grid([["1/2", "-3", "x^1"]])
    assert "<td>1/2</td>" in got and "<td>-3</td>" in got and "<td>x^1</td>" in got
    assert "<td>&lt;b&gt;</td>" in matrix_grid([["<b>"]])   # never raw markup
