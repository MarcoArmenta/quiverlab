"""pyproject is PEP 639-modern (no deprecated license table / classifier), the
version is single-sourced, and the extras/urls are present (Plan 08 Task 3).

Also guards the load-bearing packaging that the Task 3 reconciliation flagged
critical: `matplotlib>=3.7` (Plan 07 viz floor) must stay in `dependencies`, and
the `[tool.setuptools.package-data]` block must keep shipping references.bib and
zoo_catalog.json (Plan 06 bibliography()/zoo() load them by __file__-relative
path; without the block a built wheel omits them and those APIs raise
FileNotFoundError in a real install)."""
import pathlib

import pytest

try:
    import tomllib
except ModuleNotFoundError:        # Python 3.10: tomllib is stdlib only from 3.11
    tomllib = pytest.importorskip("tomli")     # backport, or skip the module on 3.10

import quiverlab

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
PP = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))


def test_license_is_spdx_expression_not_table():
    lic = PP["project"]["license"]
    assert lic == "MIT", "use the PEP 639 SPDX string, not {text=...}"
    assert PP["project"]["license-files"] == ["LICENSE"]


def test_no_deprecated_license_classifier():
    bad = [c for c in PP["project"]["classifiers"] if c.startswith("License ::")]
    assert bad == [], f"drop the deprecated License classifier(s): {bad}"


def test_version_single_source():
    assert PP["project"]["version"] == quiverlab.__version__


def test_extras_and_urls_present():
    extras = PP["project"]["optional-dependencies"]
    assert {"fast", "qpa", "docs", "dev"} <= set(extras)
    assert any("passagemath-gap" in d for d in extras["qpa"])
    urls = PP["project"]["urls"]
    assert urls["Documentation"] == "https://marcoarmenta.github.io/quiverlab/"
    assert urls["Repository"] == "https://github.com/MarcoArmenta/quiverlab"


def test_matplotlib_floor_not_regressed():
    """matplotlib>=3.7 is a Plan-07 hard dep; never drop or downgrade it."""
    deps = PP["project"]["dependencies"]
    assert "matplotlib>=3.7" in deps, (
        "matplotlib>=3.7 must stay in dependencies (Plan 07 viz floor); "
        f"got dependencies={deps}"
    )


def test_package_data_ships_bib_and_zoo():
    """The load-bearing package-data block must keep shipping the data files
    that bibliography()/bibtex()/zoo() load by __file__-relative path."""
    pkg_data = PP["tool"]["setuptools"]["package-data"]
    assert pkg_data["quiverlab.citations"] == ["references.bib"]
    assert pkg_data["quiverlab.families"] == ["zoo_catalog.json"]


def test_numba_floor_not_regressed():
    """numba>=0.64 in the fast extra; kept identical to merged main."""
    assert "numba>=0.64" in PP["project"]["optional-dependencies"]["fast"]
