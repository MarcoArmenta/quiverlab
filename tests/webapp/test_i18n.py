import logging

from webapp.server.i18n import LANGS, all_keys, catalog, t


def test_langs():
    assert LANGS == ("en", "es")


def test_en_es_key_parity():
    en, es = set(catalog("en")), set(catalog("es"))
    assert en == es, "en/es key mismatch: " + str(sorted(en ^ es))


def test_translates_known_keys():
    assert t("inv.coxeter", "en") == "Coxeter polynomial"
    assert t("inv.coxeter", "es") == "Polinomio de Coxeter"
    assert t("form.pdf", "es") == "reporte de pasos detallados"


def test_unknown_lang_uses_english():
    assert t("form.compute", "fr") == "Compute"


def test_missing_key_returns_key_and_warns(caplog):
    with caplog.at_level(logging.WARNING):
        assert t("does.not.exist", "es") == "does.not.exist"
    assert any("does.not.exist" in r.message for r in caplog.records)


def test_no_placeholder_spanish():
    for k, v in catalog("es").items():
        assert v.strip(), f"empty Spanish for {k}"
        assert "TODO" not in v and "traducir" not in v.lower(), f"placeholder in {k}"
