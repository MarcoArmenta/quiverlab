"""Macaulay2 subprocess session: availability probe + script runner.

Unlike the QPA bridge (in-process libgap), M2 is driven as a subprocess:
each call writes the script to a temp file and runs ``M2 --script file``.
Importing this module never probes the binary (lazy, like qpa.session).
"""
from __future__ import annotations

import functools
import os
import shutil
import subprocess
import tempfile

from quiverlab.errors import M2UnavailableError

_HINT = ("Macaulay2 not found. Install it: macOS `brew tap Macaulay2/tap && "
         "brew install M2`; Ubuntu `sudo add-apt-repository ppa:macaulay2/macaulay2 "
         "&& sudo apt install macaulay2`. The m2 test bucket skips without it "
         "unless QUIVERLAB_REQUIRE_M2=1.")


def _which_m2():
    return shutil.which("M2")


@functools.lru_cache(maxsize=1)
def m2_available() -> bool:
    """True iff the M2 binary is on PATH and answers --version. Cached."""
    exe = _which_m2()
    if exe is None:
        return False
    try:
        out = subprocess.run([exe, "--version"], capture_output=True,
                             encoding="utf-8", timeout=30)
    except (OSError, subprocess.SubprocessError):
        return False
    return out.returncode == 0


def require_m2() -> None:
    """Raise M2UnavailableError (with a fix-it hint) unless M2 is live."""
    if not m2_available():
        raise M2UnavailableError(_HINT)


def m2_version() -> str:
    """The M2 version string (requires M2)."""
    require_m2()
    out = subprocess.run([_which_m2(), "--version"], capture_output=True,
                         encoding="utf-8", timeout=30)
    return out.stdout.strip() or out.stderr.strip()


def run_script(script: str, timeout: int = 120) -> str:
    """Run an M2 script headless; return raw stdout. Loud on failure."""
    require_m2()
    with tempfile.NamedTemporaryFile("w", suffix=".m2", delete=False,
                                     encoding="utf-8") as f:
        f.write(script)
        path = f.name
    try:
        out = subprocess.run([_which_m2(), "--script", path],
                             capture_output=True, encoding="utf-8",
                             timeout=timeout)
    finally:
        os.unlink(path)
    if out.returncode != 0:
        raise RuntimeError(f"M2 script failed (exit {out.returncode}):\n"
                           f"{out.stderr}")
    return out.stdout


def should_skip_m2() -> bool:
    """Skip predicate for tests/m2 (collection-time, like should_skip_qpa):
    skip when M2 is absent -- unless QUIVERLAB_REQUIRE_M2=1 (CI), where
    absence must fail loudly, not skip."""
    return not m2_available() and os.environ.get("QUIVERLAB_REQUIRE_M2") != "1"
