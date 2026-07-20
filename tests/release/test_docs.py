"""The docs site config is present and canonical; a real `mkdocs build` succeeds
when the [docs] extra is installed (else the build assertion skips)."""
import importlib.util
import pathlib
import shutil
import subprocess

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
DOCS_URL = "https://marcoarmenta.github.io/quiverlab/"


def test_mkdocs_config_is_canonical():
    y = (ROOT / "mkdocs.yml").read_text()
    assert f"site_url: {DOCS_URL}" in y
    assert "gen_ref_pages.py" in y and "mkdocstrings" in y and "mkdocs-jupyter" in y
    assert "exclude_docs" in y and "plans/" in y            # plan/spec md kept out
    for ch in ("internals/01-exact-fields.md", "tutorials/03-hochschild.ipynb"):
        assert ch in y


def test_gen_ref_script_present():
    assert (ROOT / "scripts" / "gen_ref_pages.py").is_file()


def test_gen_bibliography_reads_packaged_bib():
    src = (ROOT / "scripts" / "gen_bibliography.py").read_text()
    assert "references_bib_path" in src and "bibliography.md" in src
    # the packaged bib is the source of truth, never a docs/-tree copy
    assert "docs/references.bib" not in (ROOT / "mkdocs.yml").read_text()


def test_nav_covers_internals_and_tutorials():
    """Every internals chapter and tutorial notebook must appear in the mkdocs nav,
    so --strict (omitted_files -> error) never trips when Plans 04-07 add chapters."""
    y = (ROOT / "mkdocs.yml").read_text()
    chapters = sorted((ROOT / "docs" / "internals").glob("[0-9][0-9]-*.md"))
    notebooks = sorted((ROOT / "docs" / "tutorials").glob("*.ipynb"))
    missing = [f"internals/{p.name}" for p in chapters if f"internals/{p.name}" not in y]
    missing += [f"tutorials/{p.name}" for p in notebooks if f"tutorials/{p.name}" not in y]
    assert not missing, f"add these to mkdocs.yml nav (else --strict fails): {missing}"


def test_mkdocs_builds_strict_when_available():
    if importlib.util.find_spec("mkdocs") is None or shutil.which("mkdocs") is None:
        import pytest
        pytest.skip("[docs] extra not installed; the acceptance task runs the real build")
    out = subprocess.run(["mkdocs", "build", "--strict", "-d", "/tmp/quiverlab_site"],
                         cwd=ROOT, capture_output=True, text=True,
                         env={"NUMBA_NUM_THREADS": "2", "OMP_NUM_THREADS": "2",
                              "MPLBACKEND": "Agg",  # headless notebook execution
                              "PATH": __import__("os").environ["PATH"]})
    assert out.returncode == 0, out.stderr[-3000:]
