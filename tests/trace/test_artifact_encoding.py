"""Every artifact is written as UTF-8, EXPLICITLY.

``Path.write_text`` / ``open()`` default to the LOCALE codec. On Windows that is
cp1252, so the report's em dashes were written as single cp1252 bytes and every
reader -- which correctly asks for utf-8 -- blew up with ``UnicodeDecodeError``.
That failed the whole Windows CI matrix (runs 30414159886 and 30455304569,
2026-07-29), and it is invisible on macOS/Linux where the locale codec IS utf-8.

Two gates: a SOURCE scan (no text write in the shipping tree may omit the
encoding) and a live round-trip in a subprocess whose locale codec is NOT utf-8
(US-ASCII on Unix via LC_ALL=C, cp1252 on Windows). Both fail without the fix.
"""
import ast
import os
import pathlib
import subprocess
import sys
import textwrap

import pytest

pytestmark = [pytest.mark.oracle_selfcert]

_ROOT = pathlib.Path(__file__).resolve().parents[2]
_TREES = ("src/quiverlab", "webapp", "docs/gui")
# Reading from Linux /proc + /sys, which is ASCII by construction and never
# round-trips through a report.
_EXEMPT = {"src/quiverlab/hpc/resources.py"}


def _writes_without_encoding(path):
    """Every ``x.write_text(...)`` / ``x.read_text(...)`` / ``open(...)`` call in the
    file that does not pass ``encoding=``, as (line, call) pairs. Binary modes are
    fine -- they carry no codec."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    bad = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        name = (fn.attr if isinstance(fn, ast.Attribute)
                else fn.id if isinstance(fn, ast.Name) else None)
        if name not in ("write_text", "read_text", "open"):
            continue
        if any(k.arg == "encoding" for k in node.keywords):
            continue
        if name == "open":                      # binary modes need no encoding
            mode = next((a for a in node.args[1:2]), None)
            if isinstance(mode, ast.Constant) and "b" in str(mode.value):
                continue
            if len(node.args) < 2:              # open(path) -- text, needs encoding
                bad.append((node.lineno, name))
                continue
            if not isinstance(mode, ast.Constant):
                continue                        # computed mode: cannot judge, skip
        bad.append((node.lineno, name))
    return bad


def test_no_shipping_text_io_relies_on_the_locale_codec():
    offenders = []
    for tree in _TREES:
        for path in sorted((_ROOT / tree).rglob("*.py")):
            rel = path.relative_to(_ROOT).as_posix()
            if rel in _EXEMPT or "__pycache__" in rel:
                continue
            for line, call in _writes_without_encoding(path):
                offenders.append("%s:%d %s()" % (rel, line, call))
    assert not offenders, (
        "text I/O without an explicit encoding= (the locale codec is cp1252 on "
        "Windows and mangles the report): " + "; ".join(offenders))


def test_report_round_trips_under_a_non_utf8_locale(tmp_path):
    """The live reproduction: run the writer in a subprocess whose LOCALE CODEC is
    not UTF-8 -- ``LC_ALL=C`` with PEP-538 coercion and PEP-540 UTF-8 mode both off
    gives US-ASCII on Unix, and Windows uses its cp1252 ANSI code page -- and check
    the report is still written and still reads back as UTF-8.

    This test is only worth having if it FAILS without the fix -- it does: dropping
    ``encoding="utf-8"`` from ``writer.write_trace`` makes the child raise
    ``UnicodeEncodeError`` on the report's em dash. (An in-process monkeypatch of
    ``locale.getpreferredencoding`` does NOT reproduce it: CPython's text layer does
    not consult that function, so such a test passes with the bug present and is
    worse than no test at all.)"""
    script = textwrap.dedent(
        """
        import json, pathlib, sys
        sys.path.insert(0, %r)
        from quiverlab import truncated_polynomial, CC
        from quiverlab.trace.recorder import Trace
        from quiverlab.trace.writer import write_trace
        A = truncated_polynomial(2, field=CC)
        tr = Trace()
        table = A.hochschild_cohomology(1, trace=tr)
        produced = pathlib.Path(write_trace(list(tr), table, algebra=A, kind="HH^",
                                            top=1, out_dir=%r))
        html = produced.read_text(encoding="utf-8")
        assert chr(0x2014) in html, "em dash lost in the round trip"
        assert produced.with_suffix(".json").read_text(encoding="utf-8")
        print("ROUNDTRIP-OK")
        """
    ) % (str(_ROOT / "src"), str(tmp_path))
    env = dict(os.environ,
               LC_ALL="C", LANG="C",          # locale codec -> US-ASCII ...
               PYTHONUTF8="0",                # ... PEP 540 UTF-8 mode off ...
               PYTHONCOERCECLOCALE="0",       # ... and PEP 538 coercion off
               PYTHONIOENCODING="utf-8")      # (only so the child can print)
    # The child script travels through argv, which the child DECODES with its own
    # filesystem encoding (ascii under LC_ALL=C on Linux) -- so it must stay ASCII;
    # non-ASCII is spelled chr(...) inside it.
    assert script.isascii()
    proc = subprocess.run([sys.executable, "-c", script], capture_output=True,
                          text=True, env=env, timeout=180)
    assert proc.returncode == 0, proc.stderr
    assert "ROUNDTRIP-OK" in proc.stdout
