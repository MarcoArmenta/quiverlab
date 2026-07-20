"""Plan 08 acceptance: the whole release surface is present and coherent."""
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent


def test_all_workflows_present():
    wf = ROOT / ".github" / "workflows"
    for name in ("ci.yml", "qpa.yml", "docs.yml", "paper.yml", "release.yml"):
        assert (wf / name).is_file(), f"missing workflow {name}"


def test_release_artifacts_present():
    for rel in ("mkdocs.yml", "paper/paper.md", "CITATION.cff", "CONTRIBUTING.md",
                "CODE_OF_CONDUCT.md", "CHANGELOG.md", "scripts/gen_ref_pages.py",
                "scripts/release_freshness.py"):
        assert (ROOT / rel).is_file(), f"missing {rel}"


def test_pyproject_release_ready():
    import tomllib
    pp = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert pp["project"]["license"] == "MIT"                 # SPDX, not table
    extras = pp["project"]["optional-dependencies"]
    assert {"fast", "qpa", "docs", "dev"} <= set(extras)


def test_qpa_backend_optional_not_imported_by_core():
    import quiverlab                                          # must import w/o GAP
    assert hasattr(quiverlab.Algebra, "crosscheck")
