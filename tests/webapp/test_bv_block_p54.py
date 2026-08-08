"""Plan 54 Task H -- the bv_operator compute kind, cross-runner + report contract.

The BV operator Delta is a member of ``PRODUCT_KINDS`` (spec.py): its ``.blocks()``
IS the served block. This pins (1) the GUI (Pyodide) runner block shape, (2) that the
GUI runner and the server core produce a BYTE-IDENTICAL block (the two runners share
``Algebra.bv_operator``), and (3) that the worked-steps report renders the Delta
matrices as indexed grids readable back out. Unmarked (extras-gated dir, the cross-
runner-pair convention).
"""
import importlib.util
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
RUNNER_PATH = ROOT / "docs" / "gui" / "runner.py"

# k[x]/(x^3) over GF(7): symmetric (commutative Frobenius) -> the Tradler route.
_KXX3 = {"schema": 2,
         "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                     "relations": ["x*x*x"], "field": {"kind": "GF", "p": 7, "n": 1}},
         "compute": ["bv_operator:0..3"], "artifacts": {"pdf": False, "tikz": False}}


def _gui_block(body, spec_str):
    spec = importlib.util.spec_from_file_location("gui_runner_bv_twin", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one(spec_str))       # full "bv_operator:0..3" spec
    assert out["ok"], out
    return out["block"]


def test_gui_bv_block_shape():
    b = _gui_block(_KXX3, "bv_operator:0..3")
    assert b["kind"] == "bv_operator"
    assert b["hh_dims"] == [3, 2, 2, 2]
    assert b["hypothesis"] == "symmetric (Tradler AIF 2008)"
    assert set(b["matrices"]) == {"1", "2", "3"}
    assert b["nakayama"]["inner"] is True
    assert b["bracket_check"]["agrees"] is True
    assert "bv_tradler" in b["references"]


def test_gui_matches_server(tmp_path):
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    server = run_spec(ComputeRequest.model_validate(_KXX3), tmp_path)["results"]["bv_operator"]
    gui = _gui_block(_KXX3, "bv_operator:0..3")
    assert json.dumps(server, sort_keys=True) == json.dumps(gui, sort_keys=True)


def test_report_renders_delta_grids():
    # The worked-steps report shows each Delta_n matrix as an indexed grid.
    import quiverlab as ql
    from quiverlab.hpc.spec import _product_object
    from quiverlab.trace.products import products_chapter
    from quiverlab.trace.render_html import render_html
    from tests.trace._matrix_grid import grids
    A = ql.truncated_polynomial(3, field=ql.GF(7))
    obj = _product_object(A, "bv_operator", 3)
    events = products_chapter(A, "bv_operator", obj)
    html = render_html(events, algebra=A)
    gs = grids(html)
    # every NON-ZERO Delta_n matrix appears as an indexed grid; an exactly-zero
    # Delta_n is STATED (d_n = 0), never drawn (the Marco report convention), so it
    # is not among the grids.
    checked = 0
    for n in sorted(obj.matrices):
        mat = [[str(x) for x in row] for row in obj.matrices[n]]
        if any(x != "0" for row in mat for x in row):
            assert mat in gs, n
            checked += 1
    assert checked >= 1
