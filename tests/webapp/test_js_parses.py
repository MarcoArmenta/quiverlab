"""Every shipped JS file must parse -- the merge-train's line-dedup ate
closing braces more than once (devil's-advocate fold, 2026-08-05)."""
import pathlib
import shutil
import subprocess

import pytest

_NODE = shutil.which("node")
_ROOT = pathlib.Path(__file__).resolve().parents[2]
_JS = sorted(
    list((_ROOT / "webapp" / "static").rglob("*.js"))
    + list((_ROOT / "docs" / "gui").glob("*.js")))


@pytest.mark.skipif(_NODE is None, reason="node not installed")
@pytest.mark.parametrize("path", _JS, ids=[str(p.relative_to(_ROOT)) for p in _JS])
def test_js_file_parses(path):
    proc = subprocess.run([_NODE, "--check", str(path)],
                          capture_output=True, encoding="utf-8")
    assert proc.returncode == 0, proc.stderr
