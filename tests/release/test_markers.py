"""The CI buckets fast/deep/qpa are a PARTITION of the collected suite: pairwise
disjoint and exhaustive; slow implies deep; no unknown-marker warnings (Task 2).

Marked `deep` so it runs once on the deep leg, not on all 12 fast matrix cells
(it shells out four collections)."""
import pathlib
import subprocess
import sys

import pytest

pytestmark = pytest.mark.deep
# Portable -- NEVER hardcode a laptop path here: CI runners check out under a
# different root and would FileNotFoundError on `cwd=ROOT`, turning the deep leg
# and lint job red on the first push. parents[2] = repo root (tests/release/<file>).
ROOT = str(pathlib.Path(__file__).resolve().parents[2])
VENV = sys.executable


def _ids(expr):
    """Collected node ids under a -m expression ('' == the whole suite)."""
    cmd = [VENV, "-m", "pytest", "-q", "--collect-only", "-p", "no:cacheprovider"]
    if expr:
        cmd += ["-m", expr]
    out = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    assert "PytestUnknownMarkWarning" not in (out.stdout + out.stderr), out.stdout + out.stderr
    return {ln.strip() for ln in out.stdout.splitlines() if "::" in ln}


def test_buckets_partition_the_suite():
    everything = _ids("")
    fast, deep, qpa = _ids("fast"), _ids("deep"), _ids("qpa")
    # pairwise disjoint
    assert not (fast & deep), sorted(fast & deep)[:5]
    assert not (fast & qpa), sorted(fast & qpa)[:5]
    assert not (deep & qpa), sorted(deep & qpa)[:5]
    # exhaustive
    assert fast | deep | qpa == everything, sorted(everything - (fast | deep | qpa))[:5]


def test_slow_is_a_subset_of_deep():
    slow, deep = _ids("slow"), _ids("deep")
    assert slow <= deep, "every `slow` test must ride the deep leg"


def test_known_anchors():
    assert "tests/test_no_floats.py::test_no_float_literals_or_calls_in_src" in _ids("fast")
