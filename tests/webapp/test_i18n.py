import logging
import re

import pytest

from webapp.server.i18n import (LANG_NAMES, LANGS, MOUNTS, PREFIXES, all_keys,
                                bare_path, catalog, lang_links, t)


def test_langs():
    assert LANGS == ("en", "es", "fr", "zh")


def test_lang_names_and_mounts_cover_every_lang():
    assert set(LANG_NAMES) == set(LANGS) == set(PREFIXES)
    assert MOUNTS == (("", "en"), ("/es", "es"), ("/fr", "fr"), ("/zh", "zh"))


@pytest.mark.parametrize("lang", LANGS[1:])
def test_key_parity_with_english(lang):
    en, other = set(catalog("en")), set(catalog(lang))
    assert en == other, f"en/{lang} key mismatch: " + str(sorted(en ^ other))


def test_translates_known_keys():
    assert t("inv.coxeter", "en") == "Coxeter polynomial"
    assert t("inv.coxeter", "es") == "Polinomio de Coxeter"
    assert t("form.pdf", "es") == "reporte de pasos detallados"
    # fr/zh keep proper names in Latin script (terminology convention).
    assert "Coxeter" in t("inv.coxeter", "fr")
    assert "Coxeter" in t("inv.coxeter", "zh")


def test_unknown_lang_uses_english():
    assert t("form.compute", "de") == "Compute"


def test_missing_key_returns_key_and_warns(caplog):
    with caplog.at_level(logging.WARNING):
        assert t("does.not.exist", "es") == "does.not.exist"
    assert any("does.not.exist" in r.message for r in caplog.records)


@pytest.mark.parametrize("lang", LANGS[1:])
def test_no_placeholder_translations(lang):
    # NOTE: markers are case-sensitive where the lowercase form is a real word
    # ("todo" is Spanish for "all"; only the shouty TODO is a placeholder).
    for k, v in catalog(lang).items():
        assert v.strip(), f"empty {lang} for {k}"
        for marker in ("TODO", "FIXME", "XXX"):
            assert marker not in v, f"placeholder in {lang} {k}: {v!r}"
        for marker in ("traducir", "traduire"):
            assert not re.search(r"\b" + marker + r"\b", v.lower()), \
                f"placeholder in {lang} {k}: {v!r}"


@pytest.mark.parametrize("lang", LANGS[1:])
def test_format_slots_survive_translation(lang):
    """Every {slot} in an English value must survive verbatim in the
    translation — the callers .replace() on those exact slots."""
    slots = lambda s: set(re.findall(r"\{[a-z_]+\}", s))
    en, other = catalog("en"), catalog(lang)
    for k, v in en.items():
        assert slots(v) == slots(other[k]), f"format slots drifted in {lang} {k}"


def test_bare_path_strips_every_prefix():
    for prefix, lang in MOUNTS:
        assert bare_path(prefix + "/draw") == "/draw"
        assert bare_path(prefix or "/") == "/"
    # a path that merely STARTS like a prefix is not one
    assert bare_path("/esoteric") == "/esoteric"


def test_lang_links_point_at_the_same_page():
    links = lang_links("/es/draw", "es")
    assert [code for code, _n, _u in links] == ["en", "fr", "zh"]
    assert dict((code, url) for code, _n, url in links) == {
        "en": "/draw", "fr": "/fr/draw", "zh": "/zh/draw"}
    names = [name for _c, name, _u in links]
    assert names == [LANG_NAMES["en"], LANG_NAMES["fr"], LANG_NAMES["zh"]]
    # the index page round-trips to each language's root
    assert dict((code, url) for code, _n, url in lang_links("/", "en")) == {
        "es": "/es", "fr": "/fr", "zh": "/zh"}


def test_all_keys_is_english_keyset():
    assert all_keys() == set(catalog("en"))
