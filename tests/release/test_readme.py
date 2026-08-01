"""README is the release front door: badges, real install, quickstart, links."""
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
README = (ROOT / "README.md").read_text(encoding="utf-8")


def test_badges_present():
    for b in ("actions/workflows/ci.yml/badge.svg", "img.shields.io/pypi/v/quiverlab",
              "License-MIT"):
        assert b in README


def test_tests_badge_matches_verification_page():
    """The README tests badge may never go stale: its count must equal the
    audited count stated on docs/verification.md (both are updated together by
    the verification-transparency standing rule), and it must link there."""
    badge = re.search(r"img\.shields\.io/badge/tests-(\d+)_oracle--pinned", README)
    assert badge, "README is missing the tests badge"
    page = (ROOT / "docs" / "verification.md").read_text(encoding="utf-8")
    stated = re.search(r"The suite is \*\*(\d+) tests\*\*", page)
    assert stated, "docs/verification.md no longer states the audited suite count"
    assert badge.group(1) == stated.group(1), (
        f"README badge says {badge.group(1)} tests but docs/verification.md "
        f"states {stated.group(1)} — update them together"
    )
    assert "marcoarmenta.github.io/quiverlab/verification/" in README


def test_readme_prose_count_matches_verification_page():
    """The prose "the suite is N tests" claim may never drift from the badge:
    like the badge, it must equal the audited count stated on docs/verification.md."""
    prose = re.search(r"the suite is (\d+) tests", README)
    assert prose, 'README no longer states "the suite is N tests" in prose'
    page = (ROOT / "docs" / "verification.md").read_text(encoding="utf-8")
    stated = re.search(r"The suite is \*\*(\d+) tests\*\*", page)
    assert stated, "docs/verification.md no longer states the audited suite count"
    assert prose.group(1) == stated.group(1), (
        f"README prose says {prose.group(1)} tests but docs/verification.md "
        f"states {stated.group(1)} — update them together"
    )


def test_install_and_quickstart():
    assert "pip install quiverlab" in README
    assert 'pip install "quiverlab[qpa]"' in README
    assert "hochschild_cohomology(3)" in README


def test_links_to_docs_and_tutorials_and_citation():
    assert "https://marcoarmenta.github.io/quiverlab/" in README
    assert "docs/tutorials" in README and "docs/internals" in README
    assert "CITATION.cff" in README
