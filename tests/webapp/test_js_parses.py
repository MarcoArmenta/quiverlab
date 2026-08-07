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


@pytest.mark.parametrize("path", _JS, ids=[str(p.relative_to(_ROOT)) for p in _JS])
def test_js_file_is_clean_text(path):
    """Shipped JS is TEXT: valid UTF-8 with no raw control bytes. V8 tolerates
    an embedded NUL, but grep/diff/patch tooling flips to binary mode on it --
    a raw 0x00 inside a string literal must be written as the \\x00 escape
    (2026-08-06: two raw NULs shipped inside the commutativity preset's
    path-class separator and turned gui.js binary for every text tool)."""
    raw = path.read_bytes()
    raw.decode("utf-8")  # raises on invalid UTF-8
    bad = sorted({b for b in raw if b < 0x20 and b not in (0x09, 0x0A, 0x0D)})
    assert not bad, f"raw control bytes {bad} in {path.name}"
