"""libgap session harness for the optional [qpa] backend (spec §5 component 12).

GAP is NOT a hard dependency: `pip install quiverlab[qpa]` pulls passagemath-gap
(prebuilt GAP + QPA + GBNP; macOS/Linux, Python 3.11-3.14). Everything here is
guarded by gap_available(); qpa-marked tests SKIP when GAP is absent (local dev)
and are MANDATORY in CI (QUIVERLAB_REQUIRE_QPA=1 -> require_gap() raises).

Windows: passagemath-gap ships no native wheel; gap_available() is False and
crosscheck() raises QpaUnavailableError pointing at WSL2 / conda-forge (loud,
graceful -- spec §5 c.12). No floats anywhere: GAP scripts are strings, results
are read back as exact integers."""
from __future__ import annotations

import functools
import os
import sys

from quiverlab.errors import QpaUnavailableError

_WINDOWS_MSG = (
    "quiverlab[qpa] needs GAP + QPA, for which passagemath-gap ships no native "
    "Windows wheel. Use WSL2 (Linux wheels work) or conda-forge gap-core + "
    "gap-pkg-qpa. The pure-Python quiverlab core is fully functional without "
    "[qpa]; crosscheck() is a validation convenience only."
)
_INSTALL_MSG = "install the backend with:  pip install 'quiverlab[qpa]'"


@functools.lru_cache(maxsize=1)
def _import_libgap():
    """Return (libgap, None) if importable, else (None, reason)."""
    if sys.platform == "win32":
        return None, _WINDOWS_MSG
    try:
        from sage.libs.gap.libgap import libgap
        return libgap, None
    except Exception as e_primary:  # noqa: BLE001
        try:
            from passagemath_gap import libgap
            return libgap, None
        except Exception:  # noqa: BLE001
            return None, f"passagemath-gap not importable ({e_primary!r}); {_INSTALL_MSG}"


@functools.lru_cache(maxsize=1)
def _qpa_loaded():
    lg, reason = _import_libgap()
    if lg is None:
        return False, reason
    try:
        ok = lg.LoadPackage("qpa")
        if ok == lg.eval("fail"):
            return False, 'GAP is present but LoadPackage("qpa") returned fail'
    except Exception as e:  # noqa: BLE001
        return False, f'LoadPackage("qpa") raised {e!r}'
    return True, None


def gap_available() -> bool:
    """True iff libgap imports AND QPA loads. Cached after the first call."""
    return _qpa_loaded()[0]


def require_gap() -> None:
    """Raise QpaUnavailableError (message + fix-it hint) unless QPA is live."""
    ok, reason = _qpa_loaded()
    if not ok:
        raise QpaUnavailableError("QPA backend unavailable", hint=reason)


def libgap_handle():
    """The live libgap handle (QPA loaded). Raises QpaUnavailableError otherwise."""
    require_gap()
    return _import_libgap()[0]


def run(script: str):
    """Eval a GAP script string in the QPA-loaded session; return the libgap value."""
    return libgap_handle().eval(script)


def should_skip_qpa() -> bool:
    """The skip predicate for qpa-marked tests. True (skip) only when GAP is absent
    AND we are NOT in the mandatory CI job. Under QUIVERLAB_REQUIRE_QPA=1 this returns
    False, so the tests RUN and fail naturally if GAP is missing/broken -- never a
    silent green skip.

    Why this and not a setup_module() escalation: a `@pytest.mark.skipif(...)` is
    evaluated at COLLECTION, before any setup_module() body runs, so a setup-time
    'enforce' can never turn a skipped test into a failure. Folding the environment
    check INTO the skip predicate is the only correct place."""
    return not gap_available() and os.environ.get("QUIVERLAB_REQUIRE_QPA") != "1"
