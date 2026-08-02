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


# --------------------------------------------------------------------------- #
# The symmetric hazard on the READING side: a TEST that reads a tracked repo
# file (gui.js / .html / .css / .json / .md) with a bare ``read_text()`` blows up
# the same way under Windows' cp1252 locale -- the em dashes in the report and the
# non-ASCII in the vendored gui.js decode-fail. This bit two p35 render tests
# (CI runs 30414159886-class, 2026-07-29/31). Guard the whole test tree so it
# cannot recur: every read of a repo-file literal must ask for utf-8.
_REPO_EXTS = (".js", ".html", ".css", ".json", ".md")


def _literal_has_repo_ext(subtree):
    """True if the path expression contains a string literal naming a tracked
    repo file (by extension) -- i.e. this read targets a repo file, not tmp_path."""
    for n in ast.walk(subtree):
        if (isinstance(n, ast.Constant) and isinstance(n.value, str)
                and n.value.endswith(_REPO_EXTS)):
            return True
    return False


def _repo_file_reads_without_encoding(path):
    """Each ``read_text(...)`` / read-mode ``open(...)`` of a repo-file literal in
    the file that omits ``encoding=``, as (line, call) pairs."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    bad = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        name = (fn.attr if isinstance(fn, ast.Attribute)
                else fn.id if isinstance(fn, ast.Name) else None)
        if name not in ("read_text", "open"):
            continue
        if any(k.arg == "encoding" for k in node.keywords):
            continue
        if name == "open":
            mode = next((a for a in node.args[1:2]), None)
            if isinstance(mode, ast.Constant) and any(
                    c in str(mode.value) for c in ("b", "w", "a")):
                continue                          # binary / write / append: not this bug
            recv = node.args[0] if node.args else None
        else:
            recv = fn.value if isinstance(fn, ast.Attribute) else None
        if recv is not None and _literal_has_repo_ext(recv):
            bad.append((node.lineno, name))
    return bad


def test_no_test_reads_a_repo_file_without_encoding():
    offenders = []
    for path in sorted((_ROOT / "tests").rglob("test_*.py")):
        if "__pycache__" in path.as_posix():
            continue
        for line, call in _repo_file_reads_without_encoding(path):
            offenders.append("%s:%d %s()"
                             % (path.relative_to(_ROOT).as_posix(), line, call))
    assert not offenders, (
        "test reads a tracked repo file without encoding='utf-8' (Windows' cp1252 "
        "locale decode-fails the report/gui.js): " + "; ".join(offenders))


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
                          text=True, encoding="utf-8", env=env, timeout=180)
    assert proc.returncode == 0, proc.stderr
    assert "ROUNDTRIP-OK" in proc.stdout


# --------------------------------------------------------------------------- #
# The subprocess READING side: a TEST that runs a child in TEXT mode
# (``text=True`` / ``universal_newlines=True``) and reads its stdout/stderr
# DECODES with the parent's LOCALE codec unless ``encoding=`` is passed. On
# Windows that codec is cp1252, so a node/python child that writes UTF-8 -- the
# report's em dash, non-ASCII in gui.js output -- decode-mangles into U+FFFD.
# That is exactly the ``assert '—' == '�'`` that failed every Windows fast
# cell on ``tests/gui/test_cayley_gui.py`` (the THIRD instance of this bug, after
# the write side and the repo-file read side above). Guard the whole test tree so
# it cannot recur: every text-mode ``subprocess.run``/``check_output`` must ask
# for utf-8.
# --------------------------------------------------------------------------- #
_SUBPROCESS_VERBS = ("run", "check_output")


def _is_true(node):
    return isinstance(node, ast.Constant) and node.value is True


def _scan_subprocess_text(tree):
    """Each ``subprocess.run(...)`` / ``subprocess.check_output(...)`` call in the
    AST that runs in TEXT mode (``text=True`` or ``universal_newlines=True``) but
    passes no ``encoding=``, as (line, verb) pairs. Byte-mode calls carry no codec
    and are fine; ``capture_output=True`` alone is bytes and is not flagged."""
    bad = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        if isinstance(fn, ast.Attribute):
            verb = fn.attr
        elif isinstance(fn, ast.Name):
            verb = fn.id
        else:
            continue
        if verb not in _SUBPROCESS_VERBS:
            continue
        kw = {k.arg: k.value for k in node.keywords if k.arg}
        text_mode = _is_true(kw.get("text")) or _is_true(kw.get("universal_newlines"))
        if not text_mode or "encoding" in kw:
            continue
        bad.append((node.lineno, verb))
    return bad


def _subprocess_text_without_encoding(path):
    return _scan_subprocess_text(ast.parse(path.read_text(encoding="utf-8")))


def test_no_test_runs_a_subprocess_in_text_mode_without_encoding():
    offenders = []
    for path in sorted((_ROOT / "tests").rglob("test_*.py")):
        if "__pycache__" in path.as_posix():
            continue
        for line, verb in _subprocess_text_without_encoding(path):
            offenders.append("%s:%d subprocess.%s()"
                             % (path.relative_to(_ROOT).as_posix(), line, verb))
    assert not offenders, (
        "text-mode subprocess.run/check_output without encoding='utf-8' (the child's "
        "UTF-8 stdout decodes with Windows' cp1252 locale and mangles em dashes / "
        "non-ASCII program output): " + "; ".join(offenders))


def test_the_subprocess_gate_catches_the_pre_fix_pattern():
    """The detector must FLAG the pre-fix pattern and CLEAR the fixed / byte-mode
    forms -- otherwise the gate above is a rubber stamp."""
    pre = "subprocess.run(cmd, capture_output=True, text=True)"
    uni = "subprocess.check_output(cmd, universal_newlines=True)"
    post = "subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8')"
    byte = "subprocess.run(cmd, capture_output=True)"
    assert _scan_subprocess_text(ast.parse(pre)), "gate missed the pre-fix text=True pattern"
    assert _scan_subprocess_text(ast.parse(uni)), "gate missed universal_newlines=True"
    assert not _scan_subprocess_text(ast.parse(post)), "gate false-flags the encoding= fix"
    assert not _scan_subprocess_text(ast.parse(byte)), "gate flags a byte-mode call"
