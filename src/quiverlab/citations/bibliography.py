"""Grouped, annotated, ITERABLE bibliography (spec §3.9). Plan 09's
references.entry_view() iterates bibliography() and reads e.key / e.formatted /
e.doi / e.arxiv / e.topic / e.annotation off each entry -- those attribute names
are the frozen contract. `formatted`, `doi`, and `arxiv` are built from
references.bib; `topic` is mapped from the registry `kind`; annotations come from
the registry."""
import re
from dataclasses import dataclass

from quiverlab.citations.registry import (
    REGISTRY, bibtex as _bibtex, reference, references_bib_path,
)

# topic = the /literature grouping label; mapped from the registry kind.
_KIND_TOPIC = {"algorithm": "Algorithms", "family": "Families",
               "field": "Finite fields", "foundation": "Foundations"}
_TOPIC_ORDER = ("Algorithms", "Families", "Finite fields", "Foundations")


@dataclass(frozen=True)
class Entry:
    """One bibliography row. Attribute names match Plan 09's entry_view() reads."""
    key: str            # the registry key -- matches table.references / A.citations()
    bibtex_key: str
    formatted: str      # human citation string built from references.bib
    doi: str | None
    arxiv: str | None
    topic: str          # mapped from kind
    annotation: str


def _bib_fields(bibtex_key: str) -> dict:
    """Pull author/title/journal/year/volume/pages/doi/note from one .bib entry."""
    text = references_bib_path().read_text(encoding="utf-8")
    m = re.search(r"@\w+\{" + re.escape(bibtex_key) + r",(.*?\n)\}", text, re.S)
    body = m.group(1) if m else ""
    out = {}
    for fm in re.finditer(r"(\w+)\s*=\s*\{(.*?)\}(?=,\s*\n|\s*\n\})", body, re.S):
        out[fm.group(1).lower()] = re.sub(r"\s+", " ", fm.group(2)).strip()
    return out


def _clean(s: str) -> str:
    return (s.replace("{", "").replace("}", "").replace("\\", "")
             .replace("--", "-").replace("~", " ").strip())


def _format(bibtex_key: str, f: dict) -> str:
    authors = _clean(f.get("author", "")).replace(" and ", "; ")
    year = f.get("year", "")
    title = _clean(f.get("title", bibtex_key))
    venue = _clean(f.get("journal") or f.get("booktitle") or f.get("publisher") or "")
    vol = f.get("volume", "")
    pages = _clean(f.get("pages", ""))
    bits = [b for b in [f"{authors} ({year})." if authors else "", f"{title}.",
                        f"{venue} {vol}".strip() + (f", {pages}" if pages else "") + "."
                        if venue else ""] if b]
    return " ".join(bits)


def _arxiv_of(f: dict) -> str | None:
    m = re.search(r"arXiv:([\w./-]+)", f.get("note", "") + " " + f.get("eprint", ""))
    return m.group(1) if m else None


def _entry(key: str) -> Entry:
    r = REGISTRY[key]
    f = _bib_fields(r.bibtex_key)
    return Entry(key=key, bibtex_key=r.bibtex_key, formatted=_format(r.bibtex_key, f),
                 doi=f.get("doi"), arxiv=_arxiv_of(f),
                 topic=_KIND_TOPIC.get(r.kind, "Other"), annotation=r.annotation)


@dataclass(frozen=True)
class Bibliography:
    keys: tuple
    _entries: tuple

    def __iter__(self):
        return iter(self._entries)             # Plan 09 consumes this

    @property
    def groups(self) -> dict:
        g = {}
        for e in self._entries:
            g.setdefault(e.topic, []).append(e)
        return g

    def bibtex(self) -> str:
        out, seen = [], set()
        for k in self.keys:
            bk = REGISTRY[k].bibtex_key
            if bk not in seen:
                seen.add(bk)
                out.append(_bibtex(k))
        return "\n\n".join(out) + ("\n" if out else "")

    def to_dict(self) -> dict:
        return {"keys": list(self.keys),
                "groups": {t: [e.__dict__ for e in es] for t, es in self.groups.items()}}

    def __str__(self) -> str:
        lines = []
        for t in _TOPIC_ORDER:
            es = self.groups.get(t)
            if not es:
                continue
            lines.append(f"## {t}")
            for e in es:
                lines.append(f"  [{e.key}] {e.formatted}")
                lines.append(f"      {e.annotation}")
            lines.append("")
        return "\n".join(lines).rstrip() + "\n"


def bibliography(keys=None, grouped=True, annotated=True) -> Bibliography:
    if keys is None:
        keys = list(REGISTRY)
    seen, ordered = set(), []
    for k in keys:
        reference(k)                      # loud on unknown
        if k not in seen:
            seen.add(k)
            ordered.append(k)
    return Bibliography(tuple(ordered), tuple(_entry(k) for k in ordered))
