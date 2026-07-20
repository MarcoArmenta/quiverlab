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
