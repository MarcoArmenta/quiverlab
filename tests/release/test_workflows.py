"""Structural assertions on the committed workflows (Plan 08 Tasks 6, 7, 8, 9, 10).
Text-based so it needs no YAML parser and dodges the `on:`->True 1.1 gotcha."""
import pathlib

WF = pathlib.Path(__file__).resolve().parent.parent.parent / ".github" / "workflows"


def _read(name):
    return (WF / name).read_text()


def test_ci_matrix_covers_os_and_python():
    ci = _read("ci.yml")
    for os_ in ("ubuntu-latest", "macos-latest", "windows-latest"):
        assert os_ in ci
    for py in ("3.10", "3.11", "3.12", "3.13"):
        assert f'"{py}"' in ci
    assert "-m fast" in ci                       # fast leg
    assert "QUIVERLAB_NO_NUMBA" in ci            # pure engine leg
    assert 'NUMBA_NUM_THREADS: "2"' in ci        # thread throttle
    assert "test_no_floats.py" in ci             # float gate in lint


def test_qpa_workflow_is_linux_and_mandatory():
    qpa = _read("qpa.yml")
    assert "ubuntu-latest" in qpa                       # Linux only (GAP wheels)
    assert "macos-latest" not in qpa                    # no macOS/Windows legs
    assert "windows-latest" not in qpa
    assert 'QUIVERLAB_REQUIRE_QPA: "1"' in qpa          # absent QPA = HARD failure
    assert "[dev,qpa]" in qpa                            # installs the [qpa] extra
    assert "-m qpa" in qpa                               # runs the qpa-marked suite
    assert "gap_available()" in qpa                      # explicit second hard gate
    assert 'NUMBA_NUM_THREADS: "2"' in qpa               # thread throttle
    # separate from ci.yml: push-to-main + manual + weekly, NOT on pull_request
    assert "workflow_dispatch:" in qpa
    assert "schedule:" in qpa
    assert "pull_request:" not in qpa
    assert "@master" not in qpa                          # pinned action versions


def test_docs_workflow_deploys_pages():
    d = _read("docs.yml")
    assert "actions/upload-pages-artifact@v5" in d
    assert "actions/deploy-pages@v5" in d
    assert "mkdocs build --strict" in d
    assert "pages: write" in d and "id-token: write" in d
