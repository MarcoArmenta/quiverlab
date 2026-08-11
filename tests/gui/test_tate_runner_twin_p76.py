"""Plan 76 Task G1 -- the BROWSER half of the `tate_hochschild` kind.

`tests/webapp/test_tate_kind_p76.py` pins the shared block builder; here we pin what only
the browser tier can get wrong: that `docs/gui/runner.py`'s OWN dispatch branch is reached
(not merely that the library function exists), that it produces the same JSON the server
tier does, that the `gui.js` wiring is complete in BOTH copies, and that a refusal comes
back as an error block rather than an exception escaping into Pyodide.

Unmarked per the Plan-32 extras-gated ruling.
"""
import importlib.util
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
RUNNER_PATH = ROOT / "docs" / "gui" / "runner.py"
GUI_DOCS = ROOT / "docs" / "gui" / "gui.js"
GUI_WEB = ROOT / "webapp" / "static" / "gui" / "gui.js"


def _runner():
    spec = importlib.util.spec_from_file_location("gui_runner_p76", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# kZ_2/J^3 spelled out as a DRAWN QUIVER -- which is what the GUI actually submits (the
# Pyodide twin refuses `family` blocks by design; those are the server tier). Two vertices,
# a 2-cycle, every path of length 3 killed: dim 6, self-injective and symmetric.
_KZ2J3 = {"kind": "quiver", "vertices": [1, 2],
          "arrows": {"a": [1, 2], "b": [2, 1]},
          "relations": ["a*b*a", "b*a*b"],
          "field": {"kind": "GF", "p": 32003}}

# 1 -> 2 -> 3 with rad^2 = 0: Gorenstein but NOT self-injective, so the native ring is the
# DEFERRED case and the request must come back as a block, not an exception.
_LINE_RADSQ = {"kind": "quiver", "vertices": [1, 2, 3],
               "arrows": {"a": [1, 2], "b": [2, 3]},
               "relations": ["a*b"],
               "field": {"kind": "GF", "p": 32003}}


def _body(top=2, algebra=None):
    return {"schema": 1, "algebra": algebra or _KZ2J3,
            "compute": ["tate_hochschild:0..%d" % top],
            "artifacts": {"pdf": False, "tikz": False}}


def _compute(mod, body):
    assert json.loads(mod.run_build(json.dumps(body)))["ok"]
    out = json.loads(mod.compute_one(body["compute"][0]))
    assert out["ok"], out
    return out["block"]


# --------------------------------------------------------------------------- #
# the twin's own dispatch branch
# --------------------------------------------------------------------------- #

def test_twin_dispatch_branch_is_reached_and_matches_the_server():
    """Goes through `compute_one`, i.e. the branch the browser actually executes."""
    from quiverlab.hpc.spec import _dispatch, parse_compute_item
    from quiverlab import NakayamaAlgebra
    from quiverlab.fields import GF

    twin = _compute(_runner(), _body(top=2))
    A = NakayamaAlgebra(n=2, l=3, cyclic=True, field=GF(32003))
    server, _ = _dispatch(A, parse_compute_item("tate_hochschild:0..2"), None, {})
    twin.pop("invariant", None)
    twin.pop("citations", None)
    server.pop("citations", None)
    assert json.dumps(twin, sort_keys=True) == json.dumps(server, sort_keys=True)


def test_twin_reports_both_halves():
    b = _compute(_runner(), _body(top=3))
    assert b["engine"] == "native" and b["scope"] == "self_injective"
    assert b["pos_dims"] == [1, 1, 1, 1]        # kZ2/J3 is symmetric: 1 in every degree
    assert b["neg_dims"] == [1, 1, 1]
    assert b["hat_hh0"] == 1 and b["agrees_from"] == 1
    assert b["ordinary_pos"][0] == 3            # HH^0 = dim Z(A), NOT the Tate degree 0


def test_twin_refusal_is_a_block_not_an_exception():
    """A degree the route cannot know is null; a scope failure is an `error` block. Either
    way `compute_one` returns ok -- an exception here would kill the Pyodide worker."""
    b = _compute(_runner(), _body(top=2, algebra=_LINE_RADSQ))
    assert b["kind"] == "tate_hochschild"
    # Not self-injective => the native ring is DEFERRED. Whatever the honest outcome is,
    # it must be DATA: either a typed error block or a route that claims nothing.
    assert "error" in b or b.get("pos_dims", [None])[0] is None


# --------------------------------------------------------------------------- #
# gui.js wiring, both copies
# --------------------------------------------------------------------------- #

def test_gui_js_wiring_complete_in_both_copies():
    for path in (GUI_DOCS, GUI_WEB):
        js = path.read_text(encoding="utf-8")
        assert 'id="qlgui-tate_hochschild"' in js, path
        assert 'id="qlgui-tate_hochschild-top"' in js, path
        assert '"tate_hochschild:0.." + el["tate_hochschild-top"].value' in js, path
        assert 'name === "tate_hochschild"' in js, path          # renderBlock branch
        assert 'tate_hochschild: { cb: "tate_hochschild"' in js, path   # probe sizing
        assert '"tate_hochschild"]}' in js, path                 # picker theme


def test_vendored_gui_js_is_byte_identical():
    assert GUI_DOCS.read_bytes() == GUI_WEB.read_bytes()


def test_renderer_prints_a_dash_for_an_absent_degree():
    """The null cells are the honest part of the block -- the browser must not print
    'null' or '0' for a degree the route does not know."""
    js = GUI_DOCS.read_text(encoding="utf-8")
    i = js.index('name === "tate_hochschild"')
    branch = js[i:i + 4000]
    assert 'v == null ? "\\u2014"' in branch


def test_i18n_keys_present_in_every_locale():
    import json as _json
    keys = ["pick.kind.tate_hochschild", "block.tate_hochschild.title",
            "block.tate_hochschild.positive", "block.tate_hochschild.negative",
            "block.tate_hochschild.agrees", "block.tate_hochschild.period"]
    for loc in ("en", "es", "fr", "zh"):
        cat = _json.loads((ROOT / "webapp" / "server" / "i18n" / f"{loc}.json")
                          .read_text(encoding="utf-8"))
        for k in keys:
            assert k in cat and cat[k].strip(), (loc, k)
