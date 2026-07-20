"""The citations registry core (spec §3.9). Real verified BibTeX; loud on unknown keys."""
import re

import pytest

from quiverlab import citations
from quiverlab.errors import CitationError


def test_every_registry_key_resolves_to_a_bibtex_entry():
    text = citations.references_bib_path().read_text(encoding="utf-8")
    bib_ids = set(re.findall(r"@\w+\{([^,]+),", text))
    for key in citations.all_keys():
        ref = citations.reference(key)
        assert ref.bibtex_key in bib_ids, f"{key} -> {ref.bibtex_key} missing from .bib"
        assert ref.annotation.strip(), f"{key} has no annotation"


def test_core_keys_present_with_right_kind():
    assert citations.reference("bardzell").kind == "algorithm"
    assert citations.reference("nakayama").kind == "family"
    assert citations.reference("conway").kind == "field"
    assert citations.reference("han_conjecture").kind == "foundation"


def test_bibtex_returns_the_entry_text():
    entry = citations.bibtex("chouhy_solotar")
    assert entry.startswith("@article{ChouhySolotar2015")
    assert "1406.2300" in entry


def test_unknown_key_fails_loudly():
    with pytest.raises(CitationError) as e:
        citations.reference("bardzel")     # typo
    assert "bardzell" in str(e.value)      # suggests the nearest known key


def test_bib_covers_registry():
    """Every registry-referenced bibtex_key resolves in references.bib, and the
    entry count is at least the current floor. This is a FLOOR, never an equality:
    Plan 08 appends software-citation entries (qpa, gap4, sagemath, quiverlab) and
    raises the floor; growth must not break this test."""
    text = citations.references_bib_path().read_text(encoding="utf-8")
    bib_ids = set(re.findall(r"@\w+\{([^,]+),", text))
    for key in citations.all_keys():
        assert citations.reference(key).bibtex_key in bib_ids
    assert len(re.findall(r"@\w+\{", text)) >= 20       # floor; Plan 08 appended qpa/gap4/sagemath/quiverlab (16 -> 20)
