"""Plan 35: static wiring of the HH product surface into the draw-page GUI.

Asserts on the JS + i18n SOURCES (the pattern of the other ``tests/gui`` string
checks) -- there is no browser here. The four product checkboxes exist with the
expected ids, ``buildRequest()`` pushes ``cup``/``cap``/``bracket``/``connes_b``
IN THAT ORDER immediately after the ``hh_homology`` push (the Task-12
curated-request order rule), and every i18n key the surface introduces is present
in BOTH locales. Byte-equality of the two ``gui.js`` copies is gated by
``tests/webapp/test_draw_page.py``; re-asserted here so a drift is caught even
when only this file is run.
"""
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
GUI_DOCS = ROOT / "docs" / "gui" / "gui.js"
GUI_WEBAPP = ROOT / "webapp" / "static" / "gui" / "gui.js"
EN = ROOT / "webapp" / "server" / "i18n" / "en.json"
ES = ROOT / "webapp" / "server" / "i18n" / "es.json"

PRODUCT_KINDS = ["cup", "cap", "bracket", "connes_b"]


def test_gui_js_copies_byte_identical():
    # The docs canvas and the vendored webapp copy must not drift (the edit
    # workflow is: edit docs/gui/gui.js, then cp it over the webapp copy).
    assert GUI_DOCS.read_bytes() == GUI_WEBAPP.read_bytes()


def test_product_checkbox_ids_present():
    src = GUI_DOCS.read_text(encoding="utf-8")
    for kind in PRODUCT_KINDS:
        assert 'id="qlgui-%s"' % kind in src, kind          # the checkbox
        assert 'id="qlgui-%s-top"' % kind in src, kind      # its degree picker
        assert '"%s"' % kind in src, kind                   # registered in `el`


def test_compute_push_order_after_hh_homology():
    src = GUI_DOCS.read_text(encoding="utf-8")
    anchor = 'compute.push("hh_homology:0.."'
    pushes = [
        'compute.push("cup:0.." + el["cup-top"].value)',
        'compute.push("cap:0.." + el["cap-top"].value)',
        'compute.push("bracket:0.." + el["bracket-top"].value)',
        'compute.push("connes_b:0.." + el["connes_b-top"].value)',
    ]
    assert anchor in src
    for p in pushes:
        assert p in src, p
    positions = [src.index(anchor)] + [src.index(p) for p in pushes]
    # strictly increasing: hh_homology first, then cup, cap, bracket, connes_b
    assert positions == sorted(positions)
    assert len(set(positions)) == len(positions)


def test_cyclic_homology_checkbox_and_push_order():
    # Plan-35 follow-up: HC checkbox + degree picker exist, are registered in `el`,
    # and buildRequest() pushes cyclic_homology immediately AFTER the connes_b push.
    src = GUI_DOCS.read_text(encoding="utf-8")
    assert 'id="qlgui-cyclic_homology"' in src
    assert 'id="qlgui-cyclic_homology-top"' in src
    assert '"cyclic_homology"' in src                    # registered in `el`
    connes = 'compute.push("connes_b:0.." + el["connes_b-top"].value)'
    hc = ('compute.push("cyclic_homology:0.." + '
          'el["cyclic_homology-top"].value)')
    assert connes in src and hc in src
    assert src.index(connes) < src.index(hc)             # HC follows connes_b


def test_cyclic_homology_i18n_present_in_both_locales():
    for path in (EN, ES):
        data = json.loads(path.read_text(encoding="utf-8"))
        assert "inv.cyclic_homology" in data, path.name


def test_i18n_keys_present_in_both_locales():
    keys = (["inv.%s" % k for k in PRODUCT_KINDS]
            + ["block.cup.title", "block.cap.title", "block.bracket.title",
               "block.connes_b.title", "block.bracket.window",
               "block.products.zero"])
    for path in (EN, ES):
        data = json.loads(path.read_text(encoding="utf-8"))
        missing = [k for k in keys if k not in data]
        assert not missing, (path.name, missing)
