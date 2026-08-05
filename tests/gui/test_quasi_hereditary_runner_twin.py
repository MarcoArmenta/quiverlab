"""Plan 47: static wiring of the quasi_hereditary compute kind into the draw-page GUI.

Asserts on the JS + i18n + runner SOURCES (the ``tests/gui`` string-check pattern) -- there
is no browser here. The checkbox exists with the expected id, ``buildRequest()`` pushes
``quasi_hereditary``, ``renderBlock`` has its branch, ``inv.quasi_hereditary`` is present in
BOTH locales, and the runner carries the ETA + snippet + compute_one branch. Byte-equality
of the two ``gui.js`` copies is re-asserted so a drift is caught even here."""
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


def test_quasi_hereditary_checkbox_and_registry():
    src = GUI_DOCS.read_text(encoding="utf-8")
    assert 'id="qlgui-quasi_hereditary"' in src                # the checkbox
    assert '"quasi_hereditary"' in src                         # registered in `el`


def test_compute_push_and_render_branch():
    src = GUI_DOCS.read_text(encoding="utf-8")
    assert '"strings", "quasi_hereditary"' in src              # in the plain-kind push-list
    assert 'name === "quasi_hereditary"' in src                # renderBlock branch


def test_runner_eta_snippet_and_dispatch():
    src = RUNNER.read_text(encoding="utf-8")
    assert '"quasi_hereditary": 0.5' in src                    # ETA_MODEL scalars
    assert '"quasi_hereditary": "A.is_quasi_hereditary()"' in src   # python_snippet call
    assert 'name == "quasi_hereditary"' in src                 # compute_one branch


def test_i18n_present_in_both_locales():
    en = json.loads(EN.read_text(encoding="utf-8"))
    es = json.loads(ES.read_text(encoding="utf-8"))
    assert en.get("inv.quasi_hereditary") and es.get("inv.quasi_hereditary")
    assert en.get("block.quasi_hereditary.title") and es.get("block.quasi_hereditary.title")
