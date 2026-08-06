"""Plan 45 / C4: static wiring of the tau_tilting compute kind + the LIVE wall-and-chamber
fan into the draw-page GUI. Asserts on the JS + i18n SOURCES (the tests/gui string-check
pattern -- there is no browser here): the checkbox + budget picker exist with the expected
ids, buildRequest() pushes 'tau_tilting:<budget>', renderBlock is wired, the SVG fan
renderer exists, both gui.js copies are byte-identical, the runner carries the ETA scalar +
snippet + compute_one branch, and inv.tau_tilting is present in BOTH locales."""
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
GUI_DOCS = ROOT / "docs" / "gui" / "gui.js"
GUI_WEBAPP = ROOT / "webapp" / "static" / "gui" / "gui.js"
RUNNER = ROOT / "docs" / "gui" / "runner.py"
EN = ROOT / "webapp" / "server" / "i18n" / "en.json"
ES = ROOT / "webapp" / "server" / "i18n" / "es.json"


def test_gui_js_copies_byte_identical():
    assert GUI_DOCS.read_bytes() == GUI_WEBAPP.read_bytes()


def test_tau_tilting_checkbox_ids_present():
    src = GUI_DOCS.read_text(encoding="utf-8")
    assert 'id="qlgui-tau_tilting"' in src                 # the checkbox
    assert 'id="qlgui-tau_tilting-budget"' in src          # its budget picker
    assert '"tau_tilting"' in src                          # registered in `el`
    assert '"tau_tilting-budget"' in src


def test_compute_push_is_a_budget_not_a_range():
    src = GUI_DOCS.read_text(encoding="utf-8")
    assert 'compute.push("tau_tilting:" + el["tau_tilting-budget"].value)' in src


def test_renderblock_and_fan_wired_in_gui_js():
    src = GUI_DOCS.read_text(encoding="utf-8")
    assert 'name === "tau_tilting"' in src                 # renderBlock branch
    assert "function renderTauTilting" in src
    assert "function renderWallAndChamber" in src          # the LIVE SVG fan


def test_runner_wiring():
    src = RUNNER.read_text(encoding="utf-8")
    assert '"tau_tilting": 2.0' in src                     # ETA_MODEL scalar
    assert '"tau_tilting": "A.exchange_graph(budget_pairs=%d)"' in src   # snippet
    assert 'name == "tau_tilting"' in src                  # compute_one branch


def test_tau_tilting_i18n_present_in_both_locales():
    en = json.loads(EN.read_text(encoding="utf-8"))
    es = json.loads(ES.read_text(encoding="utf-8"))
    for key in ("inv.tau_tilting", "tt.title", "tt.four_way", "tt.budget"):
        assert key in en and en[key], key
        assert key in es and es[key], key
