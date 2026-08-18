"""Plan 57: static wiring of the radical_filtration + ar_invariants compute kinds
into the draw-page GUI, and the runner-twin parity contract. Asserts on the JS +
runner SOURCES (the tests/gui string-check pattern -- there is no browser here): the
checkboxes + budget pickers, the budget push (not a degree range), the renderBlock
branches, the picker/THEMES/SEARCH wiring, both gui.js copies byte-identical, and the
Pyodide runner's parser + dispatch + snippet branches; plus i18n in all four locales."""
import importlib.util
import json
import os
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
GUI_DOCS = ROOT / "docs" / "gui" / "gui.js"
GUI_WEBAPP = ROOT / "webapp" / "static" / "gui" / "gui.js"
RUNNER = ROOT / "docs" / "gui" / "runner.py"
I18N = ROOT / "webapp" / "server" / "i18n"


def test_gui_js_copies_byte_identical():
    assert GUI_DOCS.read_bytes() == GUI_WEBAPP.read_bytes()


def test_checkbox_ids_present():
    src = GUI_DOCS.read_text(encoding="utf-8")
    for kind in ("radical_filtration", "ar_invariants"):
        assert 'id="qlgui-%s"' % kind in src
        assert 'id="qlgui-%s-budget"' % kind in src
        assert '"%s"' % kind in src
        assert '"%s-budget"' % kind in src


def test_compute_push_is_a_budget_not_a_range():
    src = GUI_DOCS.read_text(encoding="utf-8")
    assert 'compute.push("radical_filtration:" + el["radical_filtration-budget"].value)' in src
    assert 'compute.push("ar_invariants:" + el["ar_invariants-budget"].value)' in src


def test_renderblock_branches_wired():
    src = GUI_DOCS.read_text(encoding="utf-8")
    assert 'name === "radical_filtration"' in src
    assert 'name === "ar_invariants"' in src
    # the disambiguation vs the Loewy radical series is stated in the code comment
    assert "radical_filtration_ss" in src


def test_picker_and_theme_wiring():
    src = GUI_DOCS.read_text(encoding="utf-8")
    assert 'radical_filtration: { cb: "radical_filtration"' in src
    assert 'ar_invariants: { cb: "ar_invariants"' in src
    # both kinds are in the module_ar theme + a SEARCH card exists
    assert '"radical_filtration", "ar_invariants"' in src


def test_runner_wiring():
    src = RUNNER.read_text(encoding="utf-8")
    assert 'name in ("radical_filtration", "ar_invariants")' in src   # parser branch
    assert 'name == "radical_filtration"' in src                      # dispatch branch
    assert 'name == "ar_invariants"' in src
    assert '"A.radical_filtration(budget_modules=%d)"' in src         # snippet
    assert '"A.ar_invariants(budget_modules=%d)"' in src


def test_i18n_pick_labels_in_all_locales():
    for lang in ("en", "es", "fr", "zh"):
        cat = json.loads((I18N / (lang + ".json")).read_text(encoding="utf-8"))
        for key in ("pick.kind.radical_filtration", "pick.kind.ar_invariants"):
            assert key in cat and cat[key].strip(), (lang, key)


def _runner():
    s = importlib.util.spec_from_file_location(
        "gui_twin_p57_wiring", os.path.join(ROOT, "docs", "gui", "runner.py"))
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    return mod


def test_runner_twin_parses_and_computes():
    mod = _runner()
    req = {"schema": 1,
           "algebra": {"kind": "quiver", "vertices": [1, 2, 3],
                       "arrows": {"a1": [1, 2], "a2": [2, 3]}, "relations": [],
                       "field": {"kind": "QQ"}},
           "compute": ["radical_filtration", "ar_invariants"]}
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    rf = json.loads(mod.compute_one("radical_filtration"))
    assert rf["ok"] and rf["block"]["nilpotency_index"] == 3
    ai = json.loads(mod.compute_one("ar_invariants"))
    assert ai["ok"] and ai["block"]["representation_directed"] is True
