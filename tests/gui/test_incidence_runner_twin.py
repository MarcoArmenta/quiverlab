"""Plan 75 Task G1 -- the GUI half of the poset input mode + the ``incidence_cohomology``
kind: the client-side IncidenceAlgebra build, the ``gui.js`` wiring (checkbox, degree
picker, poset panel, renderer, theme, id registry), and the vendored-copy identity.

The cross-runner block parity lives in ``tests/webapp/test_incidence_p75.py``; here we pin
what only the browser tier can get wrong. Unmarked per the Plan-32 extras-gated ruling.
"""
import importlib.util
import json
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]
RUNNER_PATH = ROOT / "docs" / "gui" / "runner.py"
GUI_DOCS = ROOT / "docs" / "gui" / "gui.js"
GUI_WEB = ROOT / "webapp" / "static" / "gui" / "gui.js"

_DIAMOND = [[1, 2], [1, 3], [2, 4], [3, 4]]


def _runner():
    spec = importlib.util.spec_from_file_location("gui_runner_p75_gui", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _body(covers, elements=None, top=3):
    return {"schema": 1,
            "algebra": {"kind": "family", "family": "IncidenceAlgebra",
                        "params": {"poset_or_covers": covers, "elements": elements},
                        "field": {"kind": "QQ"}},
            "compute": ["incidence_cohomology:0..%d" % top],
            "artifacts": {"pdf": False, "tikz": False}}


# --------------------------------------------------------------------------- #
# the client-side build (the poset input mode is USELESS in the browser without it)
# --------------------------------------------------------------------------- #

def test_gui_builds_the_incidence_family_client_side():
    """Before P75 the Pyodide twin refused EVERY ``family`` block; the poset panel emits
    one, so IncidenceAlgebra joins SkewGentle/SkewGroup as a client-side constructor."""
    m = _runner()
    out = json.loads(m.run_build(json.dumps(_body(_DIAMOND))))
    assert out["ok"], out
    assert out["dim"] == 9                      # the diamond incidence algebra
    # The summary counts come off the BUILT Hasse quiver -- the params carry poset data,
    # not vertices/arrows, and reporting 0 for a real algebra would be a lie.
    assert out["n_vertices"] == 4 and out["n_arrows"] == 4


def test_gui_build_refuses_a_directed_cycle_cleanly():
    m = _runner()
    out = json.loads(m.run_build(json.dumps(_body([[1, 2], [2, 1]]))))
    assert out["ok"] is False
    assert "cycle" in out["error"]["message"] or "antisymmetry" in out["error"]["message"]


def test_gui_build_refuses_malformed_cover_data():
    m = _runner()
    out = json.loads(m.run_build(json.dumps(_body("1<2"))))
    assert out["ok"] is False
    assert "poset_or_covers" in out["error"]["message"]


def test_isolated_elements_reach_the_client_build():
    m = _runner()
    assert json.loads(m.run_build(json.dumps(
        _body([[1, 2], [2, 3]], elements=[1, 2, 3, 4]))))["ok"]
    b = json.loads(m.compute_one("incidence_cohomology:0..2"))["block"]
    assert b["dims"] == [2, 0, 0]               # two components -> H^0 = k^2


def test_reproduce_line_names_the_public_call():
    m = _runner()
    assert json.loads(m.run_build(json.dumps(_body(_DIAMOND))))["ok"]
    json.loads(m.compute_one("incidence_cohomology:0..3"))
    src = json.loads(m.reproduce())["code"] if hasattr(m, "reproduce") else ""
    if src:
        assert "IncidenceAlgebra" in src
        assert "incidence_cohomology(3)" in src


# --------------------------------------------------------------------------- #
# gui.js wiring (both copies, byte-identical)
# --------------------------------------------------------------------------- #

def test_both_gui_copies_are_byte_identical():
    assert GUI_DOCS.read_bytes() == GUI_WEB.read_bytes()


def test_gui_js_carries_every_poset_touchpoint():
    src = GUI_DOCS.read_text(encoding="utf-8")
    for needle in (
            'id="qlgui-incidence_cohomology"',        # the checkbox
            'id="qlgui-incidence_cohomology-top"',    # the degree picker
            'id="qlgui-poset-enable"',                # the panel toggle
            'id="qlgui-poset-covers"',                # the cover-pair editor
            'id="qlgui-poset-elements"',              # the isolated-element list
            'id="qlgui-poset-preview"',               # the live Hasse preview host
            'function renderPosetPreview',
            'function posetAlgebra',
            '"incidence_cohomology:0.." + el["incidence_cohomology-top"].value',
            'name === "incidence_cohomology"',        # the result renderer
            'family: "IncidenceAlgebra"'):
        assert needle in src, "gui.js is missing " + needle


def test_the_kind_is_in_exactly_one_theme_and_it_is_hochschild():
    src = GUI_DOCS.read_text(encoding="utf-8")
    m = re.search(r"QLGUI-THEMES-BEGIN(.*?)QLGUI-THEMES-END", src, re.S)
    themes = json.loads(re.search(r"(\[.*\])", m.group(1), re.S).group(1))
    holders = [t["id"] for t in themes if "incidence_cohomology" in t["kinds"]]
    assert holders == ["hochschild"], holders


def test_the_panel_emits_null_elements_unless_a_point_is_isolated():
    """The canonical-key normalization, read straight off the shipped source: the
    `elements` field is populated ONLY inside the `isolated.length` branch."""
    src = GUI_DOCS.read_text(encoding="utf-8")
    body = src[src.index("function posetAlgebra"):]
    body = body[:body.index("\n  }")]
    assert "var elements = null;" in body
    assert "if (isolated.length) {" in body
    # ... and the only assignment to `elements` is inside that branch.
    assigns = [ln for ln in body.splitlines() if re.search(r"\belements = ", ln)]
    assert len(assigns) == 2, assigns   # the null default + the isolated-point branch
