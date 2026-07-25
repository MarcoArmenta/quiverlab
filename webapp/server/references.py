"""Bibliography adapter. Everything here is SOURCED from the library
(``quiverlab.bibliography()`` + ``table.references`` / ``A.citations()``); the
webapp only normalizes, links, and groups -- it never authors citation content
(spec 16).

``bibliography()`` returns a ``Bibliography`` (amendment pt 8) that is directly
iterable over ``Entry`` records (``key``, ``bibtex_key``, ``formatted``,
``annotation``, ``doi``, ``arxiv``, ``topic``); it also exposes ``.groups`` and
``.bibtex()``. ``entry_view`` reads either a dict or an object defensively, so a
future shape change degrades gracefully instead of crashing."""
from __future__ import annotations

from functools import lru_cache

import quiverlab as ql


def _get(entry, name: str, default=None):
    if isinstance(entry, dict):
        return entry.get(name, default)
    return getattr(entry, name, default)


def entry_view(entry) -> dict:
    key = _get(entry, "key") or str(entry)
    formatted = _get(entry, "formatted") or _get(entry, "text") or key
    doi = _get(entry, "doi")
    arxiv = _get(entry, "arxiv") or _get(entry, "arxiv_id")
    return {
        "key": key,
        "bibtex_key": _get(entry, "bibtex_key"),
        "formatted": formatted,
        "doi": doi,
        "arxiv": arxiv,
        "doi_url": (f"https://doi.org/{doi}" if doi else None),
        "arxiv_url": (f"https://arxiv.org/abs/{arxiv}" if arxiv else None),
        "topic": _get(entry, "topic") or "other",
        "annotation": _get(entry, "annotation") or _get(entry, "note") or "",
    }


@lru_cache(maxsize=1)
def _index() -> dict:
    return {v["key"]: v for v in (entry_view(e) for e in ql.bibliography())}


def resolve_references(keys) -> list:
    """Map used BibTeX keys to library entries, preserving order. An unknown
    key (should not happen if the library is consistent) degrades to a minimal
    entry rather than dropping the citation."""
    idx = _index()
    out = []
    for k in keys:
        out.append(idx.get(k, {"key": k, "bibtex_key": None, "formatted": k,
                               "doi": None, "arxiv": None, "doi_url": None,
                               "arxiv_url": None, "topic": "other",
                               "annotation": ""}))
    return out


@lru_cache(maxsize=1)
def grouped_bibliography() -> list:
    """The full bibliography grouped by topic, for the /literature page."""
    groups: dict = {}
    for e in (entry_view(x) for x in ql.bibliography()):
        groups.setdefault(e["topic"], []).append(e)
    return [{"topic": topic, "entries": sorted(entries, key=lambda x: x["key"])}
            for topic, entries in sorted(groups.items())]
