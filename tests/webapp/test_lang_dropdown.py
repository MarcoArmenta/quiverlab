"""Header language dropdown (Feature 1).

Marco: "make a language dropdown menu because the menu at the top is getting
crowded." The inline per-language links in base.html are replaced by a compact
menu. Because the strict CSP (script-src 'self', no 'unsafe-inline'; see
security.py) forbids inline <script> and inline onchange handlers -- and a
<select> cannot navigate without one, which test_pages.py's no-inline-script
guard would also catch -- the menu is a native <details>/<summary> disclosure:
zero JS, opens on click, lists the four languages as links, degrades gracefully.

These assertions are template-level (TestClient render), mirroring
tests/webapp/test_pages.py. tests/webapp/test_i18n.py (orchestrator-owned) still
pins i18n.lang_links itself.
"""
import re

import pytest
from fastapi.testclient import TestClient

from webapp.server.app import create_app
from webapp.server.config import Config
from webapp.server.i18n import LANG_NAMES, LANGS, lang_links


def _client(tmp_path):
    return TestClient(create_app(Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})))


def _header(html: str) -> str:
    """The <header>...</header> region only -- so 'no <select>' is asserted against
    the menu, not against the draw page's many compute <select>s below it."""
    m = re.search(r"<header>(.*?)</header>", html, re.DOTALL)
    assert m, "no <header> in the rendered page"
    return m.group(1)


# Every language mount, plus the draw page under two prefixes (the canvas page is
# where Feature 2's data-offline lives, so the header must render there too).
_PAGES = [("/", "en"), ("/es", "es"), ("/fr", "fr"), ("/zh", "zh"),
          ("/about", "en"), ("/es/about", "es"),
          ("/draw", "en"), ("/es/draw", "es"), ("/fr/draw", "fr"), ("/zh/draw", "zh")]


@pytest.mark.parametrize("path,lang", _PAGES)
def test_menu_is_a_details_disclosure_not_a_select(tmp_path, path, lang):
    header = _header(_client(tmp_path).get(path).text)
    assert '<details class="lang-menu">' in header
    # A native disclosure, never a <select> that would need CSP-forbidden JS.
    assert "<select" not in header


@pytest.mark.parametrize("path,lang", _PAGES)
def test_all_four_languages_present_current_selected(tmp_path, path, lang):
    header = _header(_client(tmp_path).get(path).text)
    # The CURRENT language's native name shows in the summary and as the marked
    # current entry -- this also pins base.html's mirrored native-name map to
    # i18n.LANG_NAMES (the source), so the two cannot drift silently.
    summary = re.search(r"<summary[^>]*>([^<]+)</summary>", header)
    assert summary, "no <summary> in the language menu"
    assert summary.group(1).strip() == LANG_NAMES[lang]
    # The current language is a marked, non-navigating entry (not a switch link).
    assert 'class="lang-current"' in header
    assert 'aria-current="true"' in header
    # All four native names appear in the menu (current + the three switch links).
    for code in LANGS:
        assert LANG_NAMES[code] in header, f"{LANG_NAMES[code]} missing on {path}"


@pytest.mark.parametrize("path,lang", _PAGES)
def test_other_language_links_target_the_same_page(tmp_path, path, lang):
    header = _header(_client(tmp_path).get(path).text)
    others = lang_links(path, lang)          # the three OTHER languages (code,name,url)
    assert [c for c, _n, _u in others] == [c for c in LANGS if c != lang]
    for code, name, url in others:
        assert f'href="{url}"' in header, f"missing switch link {url} on {path}"
        assert f'hreflang="{code}"' in header
        assert f'lang="{code}"' in header


def test_index_menu_roots_are_correct(tmp_path):
    # A concrete, non-parametrised anchor: the English index links to each other
    # language's ROOT (not "/en/..."), exactly as before the dropdown.
    header = _header(_client(tmp_path).get("/").text)
    for url in ('href="/es"', 'href="/fr"', 'href="/zh"'):
        assert url in header
    # ...and the Spanish index links back to the English root.
    es_header = _header(_client(tmp_path).get("/es").text)
    assert 'href="/"' in es_header
    assert 'lang="en"' in es_header


def test_no_untranslated_key_leaks_via_the_menu(tmp_path):
    # The menu adds no catalog keys (native names + lang_links only), so the /es
    # page must still be free of any leaked catalog key (mirrors test_pages.py).
    from webapp.server.i18n import all_keys
    html = _client(tmp_path).get("/es").text
    for key in all_keys():
        assert key not in html, f"key leaked into /es via the header: {key}"
