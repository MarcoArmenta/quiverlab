"""Plan 42: static wiring of the ss_hochschild compute kind into the draw-page GUI.

Asserts on the JS + i18n SOURCES (the ``tests/gui`` string-check pattern) -- there
is no browser here. The ``(b, B)`` spectral-sequence checkbox exists with the
expected ids, ``buildRequest()`` pushes ``ss_hochschild`` immediately AFTER the
``cyclic_homology`` push, and ``inv.ss_hochschild`` is present in BOTH locales.
Byte-equality of the two ``gui.js`` copies is re-asserted here so a drift is caught
even when only this file is run."""
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


def test_ss_hochschild_checkbox_ids_present():
    src = GUI_DOCS.read_text(encoding="utf-8")
    assert 'id="qlgui-ss_hochschild"' in src               # the checkbox
    assert 'id="qlgui-ss_hochschild-top"' in src           # its degree picker
    assert '"ss_hochschild"' in src                        # registered in `el`
    assert '"ss_hochschild-top"' in src


def test_compute_push_order_after_cyclic_homology():
    src = GUI_DOCS.read_text(encoding="utf-8")
    hc = ('compute.push("cyclic_homology:0.." + '
          'el["cyclic_homology-top"].value)')
    ss = ('compute.push("ss_hochschild:0.." + '
          'el["ss_hochschild-top"].value)')
    assert hc in src and ss in src
    assert src.index(hc) < src.index(ss)                   # ss follows cyclic homology


def test_renderblock_wired_in_gui_js():
    src = GUI_DOCS.read_text(encoding="utf-8")
    assert 'name === "ss_hochschild"' in src               # renderBlock branch


def test_eta_and_snippet_wired_in_runner():
    src = RUNNER.read_text(encoding="utf-8")
    assert '"ss_hochschild": 2.0' in src                   # ETA_MODEL scalars
    assert '"ss_hochschild": "A.hochschild_bB_ss(%d)"' in src   # python_snippet call
    assert 'name == "ss_hochschild"' in src                # compute_one branch


def test_ss_hochschild_i18n_present_in_both_locales():
    en = json.loads(EN.read_text(encoding="utf-8"))
    es = json.loads(ES.read_text(encoding="utf-8"))
    assert "inv.ss_hochschild" in en and en["inv.ss_hochschild"]
    assert "inv.ss_hochschild" in es and es["inv.ss_hochschild"]
