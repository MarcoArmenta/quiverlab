"""paper.md conforms to JOSS: required sections present, word count in band, and every
cited @key resolves in the single PACKAGED bib (Plan 08 Task 9), read via the library's
own accessor (never a docs/ path)."""
import pathlib
import re

from quiverlab.citations import references_bib_path

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
PAPER = (ROOT / "paper" / "paper.md").read_text()
BIB = references_bib_path().read_text(encoding="utf-8")

_REQUIRED = ["# Summary", "# Statement of need", "# State of the field",
             "# Software design", "# Research impact statement",
             "# AI usage disclosure", "# Acknowledgements", "# References"]


def _body(text):
    # drop the YAML front matter (it contains an @-email that is not a citation)
    return text.split("---", 2)[-1] if text.lstrip().startswith("---") else text


def test_required_joss_sections_present_in_order():
    idx = [PAPER.find(h) for h in _REQUIRED]
    assert all(i >= 0 for i in idx), [h for h, i in zip(_REQUIRED, idx) if i < 0]
    assert idx == sorted(idx), "JOSS sections out of order"


def test_word_count_in_joss_band():
    body = re.sub(r"```.*?```", "", _body(PAPER), flags=re.S)
    words = len(re.findall(r"\S+", body))
    assert 750 <= words <= 1750, f"paper body is {words} words (JOSS wants 750-1750)"


def test_every_citation_key_resolves_in_packaged_bib():
    keys = set(re.findall(r"@([A-Za-z][A-Za-z0-9_]+)", _body(PAPER)))   # body only, not the email
    bib_keys = set(re.findall(r"@\w+\{([^,]+),", BIB))
    missing = sorted(keys - bib_keys)
    assert not missing, f"paper cites keys absent from the packaged references.bib: {missing}"


def test_paper_bib_is_not_committed_separately():
    assert not (ROOT / "paper" / "paper.bib").exists() or \
        "paper/paper.bib" in (ROOT / ".gitignore").read_text()
