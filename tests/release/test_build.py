"""The wheel/sdist build cleanly and pass twine check (Plan 08 Task 10). Skips if
build/twine are not installed; the acceptance task installs [dev] and runs it."""
import importlib.util
import pathlib
import subprocess
import sys
import tarfile
import zipfile

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent


def _have(mod):
    return importlib.util.find_spec(mod) is not None


def test_build_and_twine_check(tmp_path):
    if not (_have("build") and _have("twine")):
        import pytest
        pytest.skip("build/twine not installed; acceptance task runs the real build")
    out = tmp_path / "dist"
    r = subprocess.run([sys.executable, "-m", "build", "--outdir", str(out)],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    dists = list(out.glob("quiverlab-*.whl")) + list(out.glob("quiverlab-*.tar.gz"))
    assert len(dists) == 2, [p.name for p in dists]
    r2 = subprocess.run([sys.executable, "-m", "twine", "check", *map(str, out.glob("*"))],
                        capture_output=True, text=True)
    assert r2.returncode == 0, r2.stdout + r2.stderr
    assert "PASSED" in r2.stdout


# Data-file (package-data) files loaded at runtime by __file__-relative path.
_DATA_FILES = ("quiverlab/citations/references.bib",
               "quiverlab/families/zoo_catalog.json")


def test_wheel_and_sdist_contain_package_data(tmp_path):
    """sdist/wheel data-file packaging check (open item #4).

    The BUILT wheel and sdist must CONTAIN the packaged data files that
    citations/registry.py (references.bib) and families/zoo.py (zoo_catalog.json)
    load by __file__-relative path -- otherwise bibliography()/zoo() break in an
    installed environment even though `twine check` still PASSES (twine checks
    metadata/readme rendering, never archive contents). This test FAILS if the
    [tool.setuptools.package-data] block is removed from pyproject.toml."""
    if not _have("build"):
        import pytest
        pytest.skip("build not installed; acceptance task runs the real build")
    out = tmp_path / "dist"
    # Robustness: a stale src/*.egg-info/SOURCES.txt freezes setuptools' file list and
    # can mask a dropped package-data block (false PASS on a dev machine). Clear it so
    # the build is honest -- CI checkouts are clean, but local runs may carry one.
    import shutil
    for egg in ROOT.glob("src/*.egg-info"):
        shutil.rmtree(egg, ignore_errors=True)
    r = subprocess.run([sys.executable, "-m", "build", "--outdir", str(out)],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]

    # wheel is a zip; package-data lands at quiverlab/<subpkg>/<file> (no src/ prefix).
    wheels = list(out.glob("quiverlab-*.whl"))
    assert wheels, [p.name for p in out.glob("*")]
    with zipfile.ZipFile(wheels[0]) as zf:
        wheel_names = set(zf.namelist())
    for rel in _DATA_FILES:
        assert rel in wheel_names, f"{rel} missing from the wheel (package-data dropped?)"

    # sdist is a tar.gz; entries are under <name-version>/src/quiverlab/...
    sdists = list(out.glob("quiverlab-*.tar.gz"))
    assert sdists, [p.name for p in out.glob("*")]
    with tarfile.open(sdists[0]) as tf:
        sdist_names = tf.getnames()
    for rel in _DATA_FILES:
        assert any(m.endswith(rel) for m in sdist_names), \
            f"{rel} missing from the sdist (package-data dropped?)"
