"""README is the release front door: badges, real install, quickstart, links."""
import pathlib

README = (pathlib.Path(__file__).resolve().parent.parent.parent / "README.md").read_text()


def test_badges_present():
    for b in ("actions/workflows/ci.yml/badge.svg", "img.shields.io/pypi/v/quiverlab",
              "License-MIT"):
        assert b in README


def test_install_and_quickstart():
    assert "pip install quiverlab" in README
    assert 'pip install "quiverlab[qpa]"' in README
    assert "hochschild_cohomology(3)" in README


def test_links_to_docs_and_tutorials_and_citation():
    assert "https://marcoarmenta.github.io/quiverlab/" in README
    assert "docs/tutorials" in README and "docs/internals" in README
    assert "CITATION.cff" in README
