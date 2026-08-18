"""Plan 77 Task 7 -- the browser half of the `koszul` kind: the Pyodide twin's OWN
dispatch, byte-parity with the server tier, and the `gui.js` wiring (checkbox, degree
picker, renderer, theme, id registry, probe map, catalog entry) in BOTH vendored copies.

The organizing invariant is that both runners route through the SINGLE shared builder
`modules.nkoszul.koszul_profile_block`, so their blocks are byte-identical by
construction; this file proves it by running the twin's own `compute_one`, not by
re-calling the library.

Unmarked per the Plan-32 extras-gated ruling (tests/gui may carry no `oracle_*` marker).
"""
import importlib.util
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
RUNNER_PATH = ROOT / "docs" / "gui" / "runner.py"
GUI_DOCS = ROOT / "docs" / "gui" / "gui.js"
GUI_WEB = ROOT / "webapp" / "static" / "gui" / "gui.js"

_KX3 = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
        "relations": ["x*x*x"], "field": {"kind": "QQ"}}


def _runner():
    spec = importlib.util.spec_from_file_location("gui_runner_p77", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _twin_block(algebra=None, compute="koszul:0..8"):
    """The block the BROWSER would render: built through the twin's own run_build +
    compute_one, never by re-calling the library (which would test nothing)."""
    m = _runner()
    body = {"schema": 1, "algebra": algebra or _KX3, "compute": [compute],
            "artifacts": {"pdf": False, "tikz": False}}
    built = json.loads(m.run_build(json.dumps(body)))
    assert built["ok"], built
    out = json.loads(m.compute_one(compute))
    assert out["ok"], out
    return out["block"]


def _server_block(algebra=None, compute="koszul:0..8"):
    import tempfile

    from quiverlab.hpc.spec import parse_request
    from quiverlab.hpc.spec import run as spec_run
    body = {"schema": 1, "algebra": algebra or _KX3, "compute": [compute],
            "artifacts": {"pdf": False, "tikz": False}}
    with tempfile.TemporaryDirectory() as td:
        return spec_run(parse_request(body), pathlib.Path(td))["results"]["koszul"]


# --------------------------------------------------------------------------- #
# the twin's OWN dispatch
# --------------------------------------------------------------------------- #

def test_twin_serves_koszul_through_its_own_dispatch():
    block = _twin_block()
    assert block["kind"] == "koszul"
    assert block["n_homogeneous"] == 3
    assert block["generation_degrees"]["1"][:4] == [0, 1, 3, 4]


def test_twin_parity_is_byte_identical_to_the_server():
    twin_block = _twin_block()
    server = _server_block()
    assert (json.dumps(twin_block, sort_keys=True, default=str)
            == json.dumps(server, sort_keys=True, default=str))


# Pi(A3) SPELLED AS A QUIVER, because the Pyodide twin accepts kind 'quiver' only (the
# Plan-09 tier boundary: families are the server tier, bar the three P75 added
# client-side). A browser user draws this by hand, so this IS the reachable path -- and
# it keeps the almost-Koszul branch inside the twin's coverage. Relations transcribed
# from the shipped PreprojectiveAlgebra('A3') presentation, verbatim.
_PI_A3 = {"kind": "quiver", "vertices": [1, 2, 3],
          "arrows": {"e12": [1, 2], "e12s": [2, 1], "e23": [2, 3], "e23s": [3, 2]},
          "relations": ["e12*e12s", "e23*e23s - e12s*e12", "-e23s*e23"],
          "field": {"kind": "QQ"}}


def test_the_hand_drawn_pi_a3_is_the_shipped_family():
    """The quiver spelling above must BE Pi(A3), not merely resemble it -- otherwise the
    almost-Koszul parity below would be testing an algebra nobody named."""
    from quiverlab import PreprojectiveAlgebra
    from quiverlab.fields import QQ
    fam = PreprojectiveAlgebra("A3", field=QQ)
    assert list(fam.quiver.vertices) == _PI_A3["vertices"]
    assert {a: [fam.quiver.source(a), fam.quiver.target(a)]
            for a in fam.quiver.arrows} == _PI_A3["arrows"]
    assert [str(r) for r in fam.relations] == _PI_A3["relations"]


def test_twin_parity_on_a_preprojective_almost_koszul():
    twin_block = _twin_block(_PI_A3, "koszul:0..6")
    server = _server_block(_PI_A3, "koszul:0..6")
    assert twin_block["almost_koszul"]["verdict"] is True
    assert (json.dumps(twin_block, sort_keys=True, default=str)
            == json.dumps(server, sort_keys=True, default=str))


def test_twin_default_top_matches_the_server_default():
    twin_block = _twin_block(compute="koszul")
    assert twin_block["certified_through_degree"] == 8
    assert (json.dumps(twin_block, sort_keys=True, default=str)
            == json.dumps(_server_block(compute="koszul"), sort_keys=True, default=str))


def test_twin_reproduce_snippet_names_the_public_method():
    m = _runner()
    assert m.ETA_MODEL["scalars"]["koszul"] >= m.ETA_MODEL["scalars"]["ext_algebra"]


# --------------------------------------------------------------------------- #
# the gui.js wiring -- both copies, byte-identical
# --------------------------------------------------------------------------- #

def test_both_gui_copies_are_byte_identical():
    assert GUI_DOCS.read_bytes() == GUI_WEB.read_bytes()


def test_gui_wires_every_koszul_touchpoint():
    js = GUI_DOCS.read_text(encoding="utf-8")
    # the checkbox + degree picker, in the SAME panel that hosts ext_algebra
    assert 'id="qlgui-koszul"' in js and 'id="qlgui-koszul-top"' in js
    # the id registry and the element list
    assert '"koszul", "koszul-top",' in js
    assert 'el.koszul, el["koszul-top"],' in js
    # the request push-list uses the ext_algebra range grammar
    assert 'compute.push("koszul:0.." + el["koszul-top"].value)' in js
    # the renderer branch and the probe/control map
    assert 'name === "koszul"' in js
    assert 'koszul: { cb: "koszul", top: "koszul-top" }' in js
    # the theme places it beside ext_algebra in the structure theme
    assert '"recognizers", "ext_algebra", "koszul", "strings"' in js
    # a catalog example that is genuinely N-Koszul (k[x]/x^3, Berger N=3)
    assert '"id": "koszul"' in js and '"koszul:0..8"' in js


def test_gui_renderer_never_prints_a_window_bounded_claim_bare():
    js = GUI_DOCS.read_text(encoding="utf-8")
    start = js.index('name === "koszul"')
    branch = js[start:start + 4200]
    # the renderer must SAY "through degree" when the window is not complete
    assert "through degree" in branch
    assert "not decided" in branch
    # ... and it must render the error block rather than a traceback
    assert "b.error" in branch


def test_i18n_has_the_picker_label_in_all_four_locales():
    for lang in ("en", "es", "fr", "zh"):
        path = ROOT / "webapp" / "server" / "i18n" / f"{lang}.json"
        cat = json.loads(path.read_text(encoding="utf-8"))
        assert cat["pick.kind.koszul"].strip()
    draw = (ROOT / "webapp" / "templates" / "draw.html").read_text(encoding="utf-8")
    assert "data-pick-kind-koszul" in draw
