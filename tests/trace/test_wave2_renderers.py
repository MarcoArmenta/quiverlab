"""Wave-2 report renderers: radical_filtration_ss, ar_quiver, derived_compare, and
the quasi_hereditary / orbit_geometry enrichment fields flow through
``results_html.results_section`` and land in the page -- non-empty chunks, never a
``None`` (which would kill the report on ``out.extend(None)``, the m0806 regression).

Renderer-contract smoke over REAL library blocks (base-dep builders); no engine
values are pinned, so it is unmarked infrastructure (not an oracle class)."""
import ast
import inspect

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.trace import results_html

# NB: do NOT touch quiverlab.verbose here -- the tests/conftest.py autouse fixture
# sets it False for these (non-verbose_default) tests and restores it, so a
# module-level write would leak into the `verbose_default` test at collection time.

_kA3 = linear_path_algebra(3, field=GF(32003))
_kA4 = linear_path_algebra(4, field=GF(32003))
_kxx = Quiver([1], {"x": (1, 1)}).algebra(relations=["x*x"], field=GF(5))


def _render(kind, block):
    chunks = results_html.results_section({kind: block})
    assert chunks, f"{kind}: no chunks"
    assert all(c is not None for c in chunks), f"{kind}: a None chunk (report would die)"
    return "".join(chunks)


def test_radical_filtration_ss_renders():
    from quiverlab.specseq.block import radical_filtration_ss_block
    page = _render("radical_filtration_ss", radical_filtration_ss_block(_kA3, 4))
    assert "Radical-filtration spectral sequence" in page
    assert "E_inf" in page or "abutment" in page.lower() or "H_n" in page


def test_radical_filtration_ss_error_renders():
    page = _render("radical_filtration_ss",
                   {"kind": "radical_filtration_ss", "top": 4,
                    "error": "needs a quiver-presented algebra", "references": []})
    assert "quiver-presented" in page


def test_ar_quiver_complete_renders():
    from quiverlab.modules.ar import ar_quiver_block
    page = _render("ar_quiver", ar_quiver_block(_kA3, budget=64))
    assert "Auslander" in page and "representation-finite" in page
    assert "dim vector" in page


def test_ar_quiver_self_injective_renders():
    # a block carrying an `error` field is caught by results_section's generic
    # "not computed — <message>" path (before the per-kind renderer), so the loud
    # self-injective refusal surfaces there.
    from quiverlab.modules.ar import ar_quiver_block
    page = _render("ar_quiver", ar_quiver_block(_kxx, budget=64))
    assert "not computed" in page and "self-injective" in page


def test_derived_compare_distinguished_renders():
    from quiverlab.derived.block import derived_compare_block
    page = _render("derived_compare", derived_compare_block(_kA3, _kA4, 4))
    assert "Derived fingerprint comparison" in page
    assert "Verdict" in page and "distinguished by" in page


def test_derived_compare_not_distinguished_renders():
    from quiverlab.derived.block import derived_compare_block
    page = _render("derived_compare", derived_compare_block(_kA3, _kA3, 4))
    assert "not distinguished by these invariants" in page
    assert "equivalent" not in page.lower()


def test_quasi_hereditary_enrichment_renders():
    from quiverlab.modules.quasihereditary import quasi_hereditary_block
    page = _render("quasi_hereditary", quasi_hereditary_block(_kA3))
    assert "Characteristic tilting module" in page
    assert "Ringel dual" in page


def test_orbit_geometry_degeneration_renders():
    from quiverlab.invariants.geometry import orbit_geometry_block
    from quiverlab.modules.morphism import direct_sum
    M = direct_sum(_kA3.simple(1), _kA3.simple(2))[0]   # dim (1,1,0), degenerate
    page = _render("orbit_geometry", orbit_geometry_block(M))
    assert "Degeneration" in page and "Hasse" in page


def test_orbit_geometry_trivial_has_no_degeneration_section():
    from quiverlab.invariants.geometry import orbit_geometry_block
    page = _render("orbit_geometry", orbit_geometry_block(_kA3.simple(2)))
    assert "Degeneration" not in page               # unique module: no order to show


def test_new_renderers_end_in_return():
    # the AST gate (test_results_html_returns.py) covers ALL helpers; pin the three new
    # ones explicitly so a future edit that drops a return fails here too.
    tree = ast.parse(inspect.getsource(results_html))
    byname = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
    for name in ("_radical_filtration_ss_html", "_ar_quiver_html", "_derived_compare_html"):
        assert name in byname, name
        assert isinstance(byname[name].body[-1], ast.Return), name
