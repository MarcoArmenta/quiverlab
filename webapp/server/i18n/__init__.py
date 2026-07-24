"""EN/ES string catalogs and the t() helper. Catalogs load once at import; a
missing key renders English (logged) rather than a blank page. The JSON files
ship inside the package — the deploy copies the whole webapp/ tree, so no
extra package-data wiring is needed for the container."""
from __future__ import annotations

import json
import logging
from pathlib import Path

_log = logging.getLogger("quiverlab_web.i18n")
LANGS = ("en", "es")
_DIR = Path(__file__).resolve().parent
_CATALOGS = {lang: json.loads((_DIR / f"{lang}.json").read_text(encoding="utf-8"))
             for lang in LANGS}


def catalog(lang: str) -> dict:
    return _CATALOGS.get(lang, _CATALOGS["en"])


def all_keys() -> set:
    return set(_CATALOGS["en"])


def t(key: str, lang: str = "en") -> str:
    cat = _CATALOGS.get(lang, _CATALOGS["en"])
    if key in cat:
        return cat[key]
    if key in _CATALOGS["en"]:
        _log.warning("i18n: key %r missing for lang %r; using English", key, lang)
        return _CATALOGS["en"][key]
    _log.warning("i18n: unknown key %r", key)
    return key
