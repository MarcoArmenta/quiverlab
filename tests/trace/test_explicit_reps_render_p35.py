"""Plan 35 UNIT 2 -- the explicit representatives are now RENDERED (not just captured):

  * the products worked-steps chapter and the Computed-results product block both lay
    out one per-degree sub-section (ordered basis enumeration -> explicit classes as
    term-sum + coordinate vector -> annihilating differential + a one-line verification
    sentence) under stable, referenceable anchors, then the structure-constant tables;
  * the HH worked-steps resolution narration gains each bar (co)chain term's ordered
    basis (reconstructed with UNIT 1's enumeration builders, shown only when the
    reconstructed length matches the recorded term dim -- never mislabelled);
  * projective/injective resolution blocks carry `term_basis` (the ordered concatenated
    path bases of the summands), rendered next to the differentials, with
    len(term_basis[n]) == the differential's column (projective) / row (injective) dim;
  * every surface tolerates a block WITHOUT the new fields (an old cache);
  * both gui.js copies wire the same degree-grouped rendering.

Rendering-only: the numbers are UNIT 1's, captured at build time; nothing recomputes.
"""
import pathlib

import pytest

import quiverlab as ql
from quiverlab.fields import QQ
from quiverlab.trace.products import products_chapter, notation_legend
from quiverlab.trace.render_html import (
    render_html, product_degree_sections)
from quiverlab.trace.results_html import results_section

ROOT = pathlib.Path(__file__).resolve().parents[2]
GUI_DOCS = ROOT / "docs" / "gui" / "gui.js"
GUI_WEBAPP = ROOT / "webapp" / "static" / "gui" / "gui.js"


def _order(html, *needles):
    """The indices of ``needles`` in ``html``; each must be present."""
    pos = []
    for s in needles:
        i = html.find(s)
        assert i >= 0, "missing %r" % s
        pos.append(i)
    return pos


# --------------------------------------------------------------------------- #
# (1) Marco 2026-08-02: the product section is a FLAT class list, then the table.
#     No per-degree sub-sections (chain enumeration / per-degree headings /
#     differentials) INSIDE the product section -- those live in the HH sections.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_products_chapter_flat_class_list_then_table():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    html = render_html(products_chapter(A, "cup", A.cup_products(2)),
                       title="cup", algebra=A)
    prod = html[html.index("id='product-steps'"):]     # the product section onward
    # the flat class list is present...
    assert "Hochschild (co)homology basis classes" in prod
    # ...with NO per-degree sub-section markers inside the product section:
    assert "id='ws-cup-hh-coh-deg-1'" not in html       # no per-degree anchor
    assert "Ordered basis of the degree" not in prod    # chain enumeration not repeated
    assert "Verification:" not in prod                  # differential dropped
    assert "Explicit representatives by degree" not in prod
    assert "Structure-constant tables" not in prod      # no intermediate heading
    # the flat list precedes the Cayley table (table follows the list, immediately)
    assert (prod.index("Hochschild (co)homology basis classes")
            < prod.index("ql-cayley"))


@pytest.mark.oracle_literature
def test_class_written_over_ordered_basis_kx2_gf7():
    """k[x]/(x^2) over GF(7): HH^1 is 1-dim, class alpha^1_1 = [x -> x], written
    directly over the ordered basis (Marco 2026-07-31: no e_k coordinate reduction --
    the coordinate vector lives in the JSON)."""
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    html = render_html(products_chapter(A, "cup", A.cup_products(1)),
                       title="cup", algebra=A)
    assert "[x → x]" in html
    assert "[x → x] = e_2" not in html          # the coordinate inline is gone


@pytest.mark.oracle_literature
def test_classes_written_over_ordered_basis_kx2_gf2():
    """Over GF(2), HH^1 is 2-dim: the classes are [x -> e_1] and [x -> x], each written
    over the ordered basis -- no ``= e_k`` coordinate tail (that is JSON-only)."""
    A = ql.truncated_polynomial(2, field=ql.GF(2))
    html = render_html(products_chapter(A, "cup", A.cup_products(1)),
                       title="cup", algebra=A)
    assert "[x → e_1]" in html
    assert "[x → x]" in html
    assert "[x → e_1] = e_1" not in html
    assert "[x → x] = e_2" not in html


@pytest.mark.oracle_selfcert
def test_computed_results_product_flat_list_and_hh_sections_unchanged():
    """The Computed-results surface shows the product family as a FLAT class list then the
    Cayley table (no per-degree product sections) -- while a plain HH block in the SAME
    report KEEPS its per-degree sections (guards against over-deleting)."""
    from quiverlab.hochschild.hh_reps import hh_reps_blocks
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    tbl = A.hochschild_cohomology(2, verbose=False)
    hh = {"kind": tbl.kind, "top": 2, "dims": list(tbl.dims), "engine": tbl.engine}
    hh.update(hh_reps_blocks(A, "hh_cohomology", 2, list(tbl.dims), tbl.engine) or {})
    html = "\n".join(results_section({"hh_cohomology": hh,
                                      "cup": A.cup_products(1).blocks()}))
    # cup: flat list, NO per-degree product anchors
    assert "Hochschild (co)homology basis classes" in html
    assert "id='cr-cup-hh-coh-deg-1'" not in html
    assert "ql-cayley" in html                               # the product table renders
    # the HH block KEEPS its per-degree structure (unchanged by this change)
    assert "id='cr-hh_cohomology-hh-coh-deg-1'" in html
    assert "Explicit representatives by degree" in html      # from the HH block


@pytest.mark.oracle_selfcert
def test_both_surfaces_render_the_flat_list_once_each():
    """A generated products report carries the worked-steps chapter AND the
    Computed-results block; each renders the flat class list ONCE, and neither emits the
    old per-degree product anchors."""
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    hp = A.cup_products(1)
    html = render_html(products_chapter(A, "cup", hp), title="cup", algebra=A,
                       results={"cup": hp.blocks()})
    assert html.count("Hochschild (co)homology basis classes") == 2   # ws + cr surfaces
    assert "hh-coh-deg-1" not in html                                 # no per-degree product anchor


@pytest.mark.oracle_selfcert
def test_connes_flat_homology_class_list():
    """connes_b (Marco 2026-08-02): ONE flat homology (z^n_j) class list, all degrees at
    once -- NO per-degree sub-sections, NO chain enumerations (those live in the HH_
    sections). The cycle class z^1_1 = e_1 (x) x is written over the ordered basis."""
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    cb = A.connes_differentials(2)
    html = render_html(products_chapter(A, "connes_b", cb), title="B", algebra=A)
    assert "Hochschild (co)homology basis classes" in html   # the flat list heading
    assert "id='ws-connes_b-hh-hom-deg-1'" not in html        # no per-degree anchor
    assert "Ordered basis of" not in html                     # no chain enumeration
    # the cycle class z^1_1 = e_1 (x) x, written over the ordered basis (no coord tail)
    assert "e_1 ⊗ x" in html
    assert "e_1 ⊗ x = e_1" not in html
    # the induced B matrices still render with rank lines
    assert "induced rank" in html


@pytest.mark.oracle_selfcert
def test_cap_flat_list_carries_cohomology_and_homology_classes():
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    html = render_html(products_chapter(A, "cap", A.cap_products(2)),
                       title="cap", algebra=A)
    # the flat list precedes the first Cayley table; bound to it so we read the LIST only
    flat = html[html.index("Hochschild (co)homology basis classes"):
                html.index('<table class="ql-matrix ql-cayley"')]
    assert r"\alpha^{1}_{1}" in flat               # a cohomology class in the flat list
    assert r"z^{1}_{1}" in flat                    # a homology class in the flat list
    assert "id='ws-cap-hh-coh-deg-1'" not in html  # no per-degree sub-sections
    assert "id='ws-cap-hh-hom-deg-1'" not in html


# --------------------------------------------------------------------------- #
# (1b) Legend now points at the explicit listings (not naming-only).
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
@pytest.mark.parametrize("kind", ["cup", "cap", "bracket", "connes_b"])
def test_legend_points_at_explicit_listings(kind):
    legend = notation_legend(kind, "", "bar/GF(7)" if kind != "connes_b" else None)
    assert "explicitly by degree" in legend


# --------------------------------------------------------------------------- #
# (1c) Tolerance: a block WITHOUT the reps fields falls back cleanly.
# --------------------------------------------------------------------------- #
def test_product_degree_sections_tolerates_missing_reps():
    assert product_degree_sections(None, None, None, "x") == []
    assert product_degree_sections({}, {}, {}, "x") == []


def test_results_product_block_without_reps_renders_tables_only():
    """A legacy/old-cache product block (no basis_classes) still renders its tables --
    no crash, no empty 'Explicit representatives' heading."""
    from quiverlab.hochschild.products import HHProducts, ProductTable
    t = ProductTable(kind="cup", degrees=(0, 0), out_degree=0, dims=(1, 1, 1),
                     constants=((("1",),),))         # [dout=1][dl=1][dr=1]
    b = HHProducts(kind="cup", top=0, tables={(0, 0): t}, engine="x",
                   basis="bar/GF(7)", window=None, references=["cup"]).blocks()
    assert "basis_classes" not in b
    html = "\n".join(results_section({"cup": b}))
    assert "Explicit representatives by degree" not in html
    assert "Hochschild (co)homology basis classes" not in html   # no flat list either
    # the big degree-major table still renders (family heading + grid)
    assert "HH^{*} \\cup HH^{*}" in html and "ql-cayley" in html


def test_products_chapter_without_productbasis_event_is_tolerated():
    """A chapter stream lacking a ProductBasis event (a legacy object) renders the
    tables without degree sections -- the ALL_EVENTS gate and _products_html tolerate
    its absence."""
    from quiverlab.trace.events import ProductStep, StepNote
    events = [StepNote("The cup product", "def", heading=True),
              ProductStep(kind="cup", degrees=(0, 0), heading=r"HH^{0} \otimes HH^{0}",
                          lines=[r"\alpha^{0}_{1} \cup \beta^{0}_{1} = \gamma^{0}_{1}"])]
    html = render_html(events, title="cup")
    assert "Explicit representatives by degree" not in html
    assert "Hochschild (co)homology basis classes" not in html   # no flat list either
    assert "gamma" in html


# --------------------------------------------------------------------------- #
# (2) HH worked-steps resolution narration: the bar (co)chain term basis.
# --------------------------------------------------------------------------- #
@pytest.mark.oracle_selfcert
def test_hh_worked_steps_carry_bar_cochain_term_basis():
    """The bar-route HH cohomology worked steps gain each term's ordered cochain basis
    (reconstructed with UNIT 1's builders; the k[x]/(x^2) degree-1 cochain space is
    {[x -> e_1], [x -> x]})."""
    A = ql.truncated_polynomial(2, field=QQ)          # non-GF(p) -> bar route fills trace
    events = []
    A.hochschild_cohomology(2, verbose=False, trace=events)
    html = render_html(events, title="HH", algebra=A)
    assert "Worked resolution steps" in html
    assert "Ordered basis of the degree-1 cochain space" in html
    # the two degree-1 cochains are listed (order-matched to the shown differential)
    p = _order(html, "Ordered basis of the degree-1 cochain space", "[x → e_1]", "[x → x]")
    assert p == sorted(p)


@pytest.mark.oracle_selfcert
def test_hh_homology_worked_steps_carry_chain_term_basis():
    A = ql.truncated_polynomial(2, field=QQ)
    events = []
    A.hochschild_homology(2, verbose=False, trace=events)
    html = render_html(events, title="HH", algebra=A)
    assert "Ordered basis of the degree-1 chain space" in html
    assert "e_1 ⊗ x" in html                          # chain label, kind="chain"


def test_bar_term_basis_omitted_without_algebra():
    """Without the algebra the worked steps cannot reconstruct a basis -- and must not
    fabricate one (no cochain-space line)."""
    A = ql.truncated_polynomial(2, field=QQ)
    events = []
    A.hochschild_cohomology(2, verbose=False, trace=events)
    html = render_html(events, title="HH", algebra=None)
    assert "Ordered basis of the degree-1 cochain space" not in html


@pytest.mark.oracle_selfcert
def test_cs_worked_steps_carry_cs_term_basis():
    """Review fix (IMPORTANT): an engine='cs' HH report used to show the differential
    grid with corner naming only and NO way to read its columns. The CS resolution term
    basis is now enumerated (rebuilt from the deterministic CS resolution, length-guarded
    against the recorded term dim). k[x]/(x^3): C^1 = {[x -> e_1], [x -> x], [x -> x*x]}."""
    A = ql.truncated_polynomial(3, field=QQ)
    events = []
    A.hochschild_cohomology(2, verbose=False, trace=events, engine="cs")
    html = render_html(events, title="HH cs", algebra=A)
    assert "Worked resolution steps" in html
    assert "Ordered basis of the degree-1 cochain space" in html
    p = _order(html, "Ordered basis of the degree-1 cochain space",
               "[x → e_1]", "[x → x]", "[x → x*x]")
    assert p == sorted(p), p


# --------------------------------------------------------------------------- #
# (3) Module resolutions: term_basis interpretability + rendering + tolerance.
# --------------------------------------------------------------------------- #
def _kA3():
    return ql.Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(field=ql.GF(7))


@pytest.mark.oracle_selfcert
def test_term_basis_lengths_match_projective_differential_cols():
    from quiverlab.hpc.spec import _term_basis_blocks
    A = _kA3()
    M = A.simple(1)
    res = M.projective_resolution(3)
    tb = _term_basis_blocks(res, "projective_resolution", M)
    assert tb is not None
    dvs = res.dimension_vectors()
    for n in range(len(res.terms)):
        D = res.differential(n)
        cols = len(D[0]) if D and D[0] else 0
        assert len(tb[n]) == cols, (n, len(tb[n]), cols)          # cols index the term
        assert len(tb[n]) == sum(dvs[n].values())                 # == term dimension


@pytest.mark.oracle_selfcert
def test_term_basis_lengths_match_injective_differential_rows():
    from quiverlab.hpc.spec import _term_basis_blocks
    A = _kA3()
    M = A.simple(3)
    res = M.injective_resolution(3)
    tb = _term_basis_blocks(res, "injective_resolution", M)
    assert tb is not None
    for n in range(len(res.terms)):
        D = res.differential(n)
        rows = len(D) if D else 0
        assert len(tb[n]) == rows, (n, len(tb[n]), rows)          # rows index E^n


@pytest.mark.oracle_selfcert
def test_injective_term_basis_content_order_pinned():
    """Review MINOR 2: pin the injective term_basis CONTENT ORDER (labels order ==
    differential row order), not only its length. E^n's dual basis IS the projective
    resolution of DM (over A^op), and the injective differential is that resolution's
    transpose -- so the injective term basis paths equal the proj-resolution-of-DM term
    basis paths IN ORDER, and each differential's ROW count is the term size (rows index
    E^n). Plus an exact label-order literal for a known injective resolution."""
    from quiverlab.hpc.spec import _term_basis_blocks
    from quiverlab.modules.duality import dualize
    A = _kA3()
    M = A.simple(3)
    inj = M.injective_resolution(3)
    tb = _term_basis_blocks(inj, "injective_resolution", M)
    assert tb[0] == ["I_3: e_3", "I_3: b", "I_3: b*a"]     # exact order, not just length
    assert tb[1] == ["I_2: e_2", "I_2: a"]
    # independent order pin: the injective differential is the transpose of the proj.
    # resolution of DM, so its rows index that resolution's term basis, in order.
    DM = dualize(M)
    ptb = _term_basis_blocks(DM.projective_resolution(3), "projective_resolution", DM)
    for n in range(len(inj.terms)):
        assert [x.split(": ", 1)[1] for x in tb[n]] == [x.split(": ", 1)[1] for x in ptb[n]]
        D = inj.differential(n)
        assert len(tb[n]) == (len(D) if D else 0)          # rows index E^n, this order


def test_products_and_connes_emit_no_per_degree_anchors():
    """Marco 2026-08-02: cup/cap/bracket AND connes_b all show a flat class list then the
    tables/matrices -- no per-degree product/Connes sub-sections, so no per-degree
    anchors (or links to them) in either surface."""
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    cb = A.connes_differentials(2)
    h2 = render_html(products_chapter(A, "connes_b", cb), title="B", algebra=A,
                     results={"connes_b": cb.blocks()})
    for token in ("ws-connes_b-hh-hom-deg-1", "cr-connes_b-hh-hom-deg-1"):
        assert token not in h2                             # no per-degree anchors/links
    hp = A.cup_products(1)
    hc = render_html(products_chapter(A, "cup", hp), title="cup", algebra=A,
                     results={"cup": hp.blocks()})
    assert "ws-cup-hh-coh-deg-1" not in hc
    assert "cr-cup-hh-coh-deg-1" not in hc


@pytest.mark.oracle_crossengine
def test_term_basis_shared_serializer_matches_across_runners(tmp_path):
    """The wheel spec and the Pyodide runner build term_basis with the SAME shape (they
    cannot import each other) -- pin key-for-key equality for a resolution block."""
    import importlib.util
    from quiverlab.hpc.spec import _term_basis_blocks as spec_tb
    runner_path = ROOT / "docs" / "gui" / "runner.py"
    spec = importlib.util.spec_from_file_location("gui_runner_tb", runner_path)
    gui = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gui)
    A = _kA3()
    M = A.simple(1)
    res = M.projective_resolution(3)
    assert spec_tb(res, "projective_resolution", M) == \
        gui._term_basis_blocks(res, "projective_resolution", M)


def test_term_basis_rendered_in_results_block():
    from quiverlab.hpc.spec import _dispatch_module, ComputeItem
    A = _kA3()
    M = A.simple(1)
    block = _dispatch_module(A, ComputeItem("projective_resolution", 0, 3),
                             M, None, None)
    assert "term_basis" in block
    html = "\n".join(results_section({"projective_resolution": block}))
    assert "Ordered basis of each resolution term" in html
    assert "P_1: a" in html                                # a path label of P_1


def test_resolution_block_without_term_basis_is_tolerated():
    """An older cached resolution block (no term_basis) renders its differentials
    without the term-basis section -- no crash."""
    block = {"kind": "projective_resolution", "top": 1,
             "terms": [{"1": 1}, {}], "summands": ["P_1", "0"], "betti": [1, 0],
             "differentials": [{"rows": 1, "cols": 1, "matrix": [["1"]]},
                               {"rows": 1, "cols": 0, "matrix": []}],
             "pd": 0}
    html = "\n".join(results_section({"projective_resolution": block}))
    assert "Ordered basis of each resolution term" not in html
    assert "d_0" in html or "d_{0}" in html                # the differential still shows


def test_term_basis_blocks_returns_none_without_path_basis():
    """A term whose summand projective has no path basis -> the guard returns None so
    the field is omitted (structure-constants algebras carry no path basis)."""
    from quiverlab.hpc.spec import _term_basis_blocks

    class _NoPathAlg:
        pass

    class _FakeM:
        algebra = _NoPathAlg()

    class _FakeRes:
        terms = [0]

        def term(self, n):
            return [1]

    assert _term_basis_blocks(_FakeRes(), "projective_resolution", _FakeM()) is None


# --------------------------------------------------------------------------- #
# (4) GUI wiring (source checks, the tests/gui convention) -- both copies mirror.
# --------------------------------------------------------------------------- #
def test_gui_js_copies_byte_identical():
    assert GUI_DOCS.read_bytes() == GUI_WEBAPP.read_bytes()


def test_gui_js_wires_product_flat_list_and_connes_flat_list():
    src = GUI_DOCS.read_text(encoding="utf-8")
    for fn in ("function appendProductReps", "function appendProductFlatClasses",
               "function appendRepsClasses", "function appendRepsDifferential",
               "function termSumText"):
        assert fn in src, fn
    # Marco 2026-07-31: the coordinate-vector inline is gone -- coordVectorText removed.
    assert "function coordVectorText" not in src
    # Marco 2026-08-02: cup/cap/bracket AND connes_b render the FLAT class list; no
    # per-degree product/Connes sub-sections (the per-degree call is gone).
    assert "appendProductFlatClasses(div, b)" in src
    assert "appendProductReps(div, b, name, false)" not in src
    assert 'appendProductReps(div, b, "connes_b", false)' not in src   # connes now flat
    # the Connes degree-link helpers are gone (no per-degree anchors to link to)
    assert "function degreeLink" not in src and "function productTableLinks" not in src


def test_gui_js_wires_term_basis_next_to_differentials():
    src = GUI_DOCS.read_text(encoding="utf-8")
    assert "function appendTermBasis" in src
    # appendTermBasis is called immediately before appendDifferentials
    i = src.index("appendTermBasis(div, b, proj);")
    j = src.index("appendDifferentials(div, b, proj);", i)
    assert 0 <= i < j
