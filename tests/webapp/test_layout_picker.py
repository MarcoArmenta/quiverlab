"""Layout C — the two-pane theme picker + cart (Marco 2026-08-06).

The picker is a VIEW over the real hidden checkbox grid, so most of its contract
is exercised by the existing draw-page / delegation / search tests staying green.
Here we pin the parts unique to the picker: the THEMES taxonomy (every compute
kind lives in exactly one theme, with Marco's explicit placement of the Hochschild
products and the cyclic block), that every theme/kind label is translated in all
four catalogs, and that the localized picker strings reach the draw page under
every language prefix.
"""
import json
import re
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from webapp.server.app import create_app
from webapp.server.config import Config
from webapp.server.i18n import LANGS, catalog

ROOT = Path(__file__).resolve().parents[2]
GUI_JS = ROOT / "webapp" / "static" / "gui" / "gui.js"

# The full compute-kind surface the picker must cover (19 algebra + 14 module).
ALL_KINDS = {
    "hh_cohomology", "hh_homology", "cup", "cap", "bracket", "cyclic_homology",
    "connes_b", "ss_hochschild", "cartan", "coxeter_polynomial",
    "global_dimension", "homological_profile", "center", "recognizers",
    "ext_algebra", "strings", "quasi_hereditary", "derived_fingerprint",
    "tau_tilting", "dimension_vector", "rad_top_soc", "tau", "tau_minus",
    "projective_dimension", "injective_dimension", "projective_resolution",
    "injective_resolution", "decompose", "almost_split", "tilting_check",
    "orbit_geometry", "ext", "tor",
    # Wave-2 surface expansion (2026-08-06): the three new compute kinds.
    "radical_filtration_ss", "ar_quiver", "derived_compare",
}


def _themes():
    src = GUI_JS.read_text(encoding="utf-8")
    m = re.search(r"QLGUI-THEMES-BEGIN(.*?)QLGUI-THEMES-END", src, re.S)
    assert m, "THEMES sentinels not found in gui.js"
    arr = re.search(r"(\[.*\])", m.group(1), re.S)
    assert arr, "THEMES array literal not found"
    return json.loads(arr.group(1))


def test_themes_cover_every_kind_exactly_once():
    themes = _themes()
    seen = []
    for t in themes:
        assert t["id"] and isinstance(t["kinds"], list) and t["kinds"]
        seen.extend(t["kinds"])
    assert len(seen) == len(set(seen)), "a kind appears in two themes: " + str(
        sorted(k for k in seen if seen.count(k) > 1))
    assert set(seen) == ALL_KINDS, (
        "picker kinds != the compute surface; missing "
        + str(sorted(ALL_KINDS - set(seen)))
        + " extra " + str(sorted(set(seen) - ALL_KINDS)))


def test_marco_taxonomy_pins():
    by_id = {t["id"]: t["kinds"] for t in _themes()}
    # Marco: cup, cap, bracket belong to Hochschild (not the cyclic block).
    for k in ("hh_cohomology", "hh_homology", "cup", "cap", "bracket"):
        assert k in by_id["hochschild"], k
    # Marco: keep the rest in cyclic; the radical-filtration SS joins it (wave 2).
    assert set(by_id["cyclic"]) == {
        "cyclic_homology", "connes_b", "ss_hochschild", "radical_filtration_ss"}
    # ar_quiver is algebra-level but sits in the AR theme; derived_compare in structure.
    assert "ar_quiver" in by_id["module_ar"]
    assert "derived_compare" in by_id["structure"]


def test_every_theme_and_kind_label_is_translated_in_all_langs():
    themes = _themes()
    keys = ["pick.theme." + t["id"] for t in themes]
    for t in themes:
        keys.extend("pick.kind." + k for k in t["kinds"])
    for lang in LANGS:
        cat = catalog(lang)
        for key in keys:
            assert key in cat and cat[key].strip(), f"{lang} missing {key}"


def _client(tmp_path):
    return TestClient(create_app(Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})))


@pytest.mark.parametrize("lang", LANGS)
def test_draw_page_carries_picker_attrs_in_every_language(lang, tmp_path):
    prefix = "" if lang == "en" else "/" + lang
    html = _client(tmp_path).get(prefix + "/draw").text
    # a theme name, a kind name, and a chrome string -- all localized data-*.
    for attr in ("data-pick-theme-hochschild", "data-pick-kind-tau_tilting",
                 "data-pick-cart-title", "data-pick-cost-fast", "data-pick-filter"):
        assert attr in html, f"{attr} absent under {prefix}/draw"
    # the localized value actually rendered (not the raw key)
    assert "pick.theme.hochschild" not in html


# --------------------------------------------------------------------------- #
# Wave-2 surface expansion (2026-08-06): QQ field, potential -> Jacobian, the
# three new kinds + the algebra-B input for derived_compare.
# --------------------------------------------------------------------------- #
APP_JS = ROOT / "webapp" / "static" / "app.js"


@pytest.mark.parametrize("lang", LANGS)
def test_draw_page_carries_wave2_attrs_in_every_language(lang, tmp_path):
    prefix = "" if lang == "en" else "/" + lang
    html = _client(tmp_path).get(prefix + "/draw").text
    for attr in ("data-pick-kind-radical_filtration_ss", "data-pick-kind-ar_quiver",
                 "data-pick-kind-derived_compare", "data-potential-label",
                 "data-potential-ph", "data-algb-legend", "data-algb-hint"):
        assert attr in html, f"{attr} absent under {prefix}/draw"
    # localized, not the raw key
    for key in ("pick.kind.ar_quiver", "draw.algb_legend", "draw.potential_label"):
        assert key not in html


def test_new_kinds_labelled_in_all_catalogs():
    for lang in LANGS:
        cat = catalog(lang)
        for key in ("pick.kind.radical_filtration_ss", "pick.kind.ar_quiver",
                    "pick.kind.derived_compare", "draw.potential_label",
                    "draw.potential_ph", "draw.potential_conflict", "draw.algb_legend",
                    "draw.algb_mode_label", "draw.algb_type_label",
                    "draw.algb_preset_label", "draw.algb_hint", "form.field_ph"):
            assert key in cat and cat[key].strip(), f"{lang} missing {key}"


def test_gui_js_wires_qq_potential_and_algebra_b():
    src = GUI_JS.read_text(encoding="utf-8")
    # QQ is a first-class field option and currentField() emits it.
    assert '<option value="QQ">' in src
    assert 'kind: "QQ"' in src
    # potential -> algebra.potential, only when non-empty
    assert "req.algebra.potential = pot" in src
    # derived_compare's second algebra
    assert "req.algebra_b = algb" in src and "function buildAlgebraB" in src


def test_app_js_accepts_qq_field_not_only_cc_gf():
    # readComputeBody must recognise QQ, not coerce every non-CC field to GF.
    src = APP_JS.read_text(encoding="utf-8")
    m = re.search(r"function readComputeBody\(\)\s*\{(.*?)\n\}", src, re.S)
    assert m, "readComputeBody not found"
    body = m.group(1)
    assert "QQ" in body, "readComputeBody still hardcodes only the CC/GF pair"
