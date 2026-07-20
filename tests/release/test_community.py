"""Community/citation files exist and the canonical docs URL is consistent across
mkdocs.yml, README, and CITATION.cff (Plan 08 Task 12)."""
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
URL = "https://marcoarmenta.github.io/quiverlab/"


def test_community_files_present():
    for f in ("CONTRIBUTING.md", "CODE_OF_CONDUCT.md", "CITATION.cff"):
        assert (ROOT / f).is_file(), f"missing {f} (JOSS review checklist)"


def test_citation_cff_shape():
    cff = (ROOT / "CITATION.cff").read_text()
    assert "cff-version: 1.2.0" in cff and "Armenta" in cff
    assert 'repository-code: "https://github.com/MarcoArmenta/quiverlab"' in cff


def test_citation_version_matches_pyproject():
    import re
    import tomllib
    ver = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
    core = re.match(r"\d+\.\d+\.\d+", ver).group(0)        # 0.1.0.dev0 -> 0.1.0
    cff = (ROOT / "CITATION.cff").read_text()
    assert f"version: {core}" in cff, f"CITATION.cff version must be the release core {core}"


def test_docs_url_consistent_across_files():
    mk = (ROOT / "mkdocs.yml").read_text()
    rd = (ROOT / "README.md").read_text()
    cff = (ROOT / "CITATION.cff").read_text()
    assert f"site_url: {URL}" in mk
    assert URL in rd
    assert f'url: "{URL}"' in cff


def test_development_page_documents_qpa_math():
    dev = (ROOT / "docs" / "development" / "release.md").read_text()
    assert "Ext^n_{A^e}(A, A)" in dev and "EnvelopingAlgebra" in dev
