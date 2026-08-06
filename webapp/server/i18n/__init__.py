"""EN/ES/FR/ZH string catalogs and the t() helper. Catalogs load once at import;
a missing key renders English (logged) rather than a blank page. The JSON files
ship inside the package — the deploy copies the whole webapp/ tree, so no
extra package-data wiring is needed for the container."""
from __future__ import annotations

import json
import logging
from pathlib import Path

_log = logging.getLogger("quiverlab_web.i18n")
LANGS = ("en", "es", "fr", "zh")
# A language's name is written in itself — it is the one string every reader of
# that language recognises, so it is static data, not a catalog key.
LANG_NAMES = {"en": "English", "es": "Español", "fr": "Français", "zh": "中文"}
# (url prefix, language) — the mount points every page gets. English owns the
# bare path; every other language lives under its /xx prefix.
MOUNTS = tuple(("" if lang == "en" else "/" + lang, lang) for lang in LANGS)
PREFIXES = {lang: prefix for prefix, lang in MOUNTS}
_DIR = Path(__file__).resolve().parent
_CATALOGS = {lang: json.loads((_DIR / f"{lang}.json").read_text(encoding="utf-8"))
             for lang in LANGS}


def bare_path(path: str) -> str:
    """The page path with any language prefix removed ("/es/draw" → "/draw")."""
    for prefix, _lang in MOUNTS:
        if prefix and (path == prefix or path.startswith(prefix + "/")):
            return path[len(prefix):] or "/"
    return path or "/"


def lang_links(path: str, lang: str):
    """``(code, native name, url)`` for every OTHER language's version of
    ``path`` — the header language menu."""
    page = bare_path(path)
    return tuple((code, LANG_NAMES[code],
                  (prefix + ("" if page == "/" else page)) or "/")
                 for prefix, code in MOUNTS if code != lang)


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
