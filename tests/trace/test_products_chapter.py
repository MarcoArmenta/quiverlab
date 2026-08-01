"""The products chapter: events render, dims drift-gate, full record (Plan 35).

`quiverlab.trace.products.products_chapter(A, kind, obj)` turns a Task-1 product
result object into the typed trace-event stream the HTML/JSON renderers turn into a
homework-grade chapter: a prose definition (StepNote), the authoritative HH
dimensions (ResultDims, which also DRIFT-GATES the narrated table dims against a
fresh Hochschild (co)homology), and one ProductStep per bidegree -- the nonzero
structure-constant equations (cup/cap/bracket) or the induced Connes matrix
(connes_b). `spec.run` backs the worked-steps bundle with the first product kind
when no HH/module trace claimed it.
"""
import pytest

import quiverlab as ql

pytestmark = [pytest.mark.oracle_selfcert]


def test_chapter_events_render_and_carry_tables(tmp_path):
    from quiverlab.trace.products import products_chapter
    from quiverlab.trace.render_html import render_html
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    hp = A.cup_products(2)
    events = products_chapter(A, "cup", hp)
    assert events, "chapter must not be empty"
    html = render_html(events, title="cup", algebra=A)
    assert "cup product" in html.lower()
    assert "HH" in html
    # every nonzero constant appears in the page
    for t in hp.tables.values():
        for mat in t.constants:
            for row in mat:
                for c in row:
                    if c != "0":
                        assert c in html


def test_drift_gate_fires_on_dim_mismatch():
    from quiverlab.trace.products import products_chapter
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    hp = A.cup_products(1)
    # sabotage a dim -- the builder must refuse to narrate inconsistent data
    hp.tables[(0, 0)] = hp.tables[(0, 0)].__class__(
        kind="cup", degrees=(0, 0), out_degree=0, dims=(99, 1, 1),
        constants=hp.tables[(0, 0)].constants)
    with pytest.raises(Exception):
        products_chapter(A, "cup", hp)


def test_cap_and_bracket_chapters_render():
    from quiverlab.trace.products import products_chapter
    from quiverlab.trace.render_html import render_html
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    for kind, obj in (("cap", A.cap_products(2)),
                      ("bracket", A.gerstenhaber_brackets(2))):
        events = products_chapter(A, kind, obj)
        assert events
        html = render_html(events, title=kind, algebra=A)
        assert "HH" in html
        for t in obj.tables.values():
            for mat in t.constants:
                for row in mat:
                    for c in row:
                        if c != "0":
                            assert c in html


def test_connes_chapter_carries_the_matrix():
    from quiverlab.trace.products import products_chapter
    from quiverlab.trace.render_html import render_html
    from tests.trace._matrix_grid import grids
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    cb = A.connes_differentials(2)
    events = products_chapter(A, "connes_b", cb)
    assert events
    html = render_html(events, title="B", algebra=A)
    # the induced B matrices are rendered as indexed grids
    got = grids(html)
    assert got, "the Connes differential matrices must render as grids"
    # B_0 = [['0', '1']] must be among them
    assert [["0", "1"]] in got


def test_connes_drift_gate_fires_on_dim_mismatch():
    from quiverlab.trace.products import products_chapter
    A = ql.truncated_polynomial(2, field=ql.GF(7))
    cb = A.connes_differentials(2)
    cb.hh_dims[0] = 99            # sabotage the recorded HH_0 dimension
    with pytest.raises(Exception):
        products_chapter(A, "connes_b", cb)


def test_products_back_the_worked_steps_bundle(tmp_path):
    """A products request with artifacts.pdf and no HH/module trace: the first
    product kind backs the worked-steps HTML+JSON bundle (spec.run selection)."""
    from quiverlab.hpc.spec import run as spec_run, parse_request
    req = parse_request({
        "schema": 2,
        "algebra": {"kind": "quiver", "vertices": [1],
                    "arrows": {"x": [1, 1]}, "relations": ["x*x*x"],
                    "field": {"kind": "GF", "p": 2, "n": 1}},
        "compute": ["cup:0..2"],
        "artifacts": {"pdf": True, "tikz": False}})
    res = spec_run(req, tmp_path)
    assert res["meta"]["pdf"] == "worked steps in trace_steps.html"
    html = (tmp_path / "trace_steps.html").read_text(encoding="utf-8")
    assert "cup product" in html.lower()
    assert (tmp_path / "trace.json").exists()
