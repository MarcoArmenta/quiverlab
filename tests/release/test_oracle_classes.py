"""The oracle-class marker counts on docs/verification.md are AUDITED against the
live pytest collection -- the badge==page doctrine (cf. test_markers.py for the
fast/deep/qpa buckets, test_readme.py for the tests badge), now applied to the
four-way oracle taxonomy (Plan 32).

For each oracle-class marker expression the page states a count; this test shells
out `pytest --collect-only -m <expr>` and requires the live number to match. A
stale table (a plan that adds a battery but forgets to bump the count) fails here.

Mechanism = OUT-OF-PROCESS collection, exactly like test_markers.py: it is proven
robust in this repo, inherits the installed extras, and stays well under 30 s. The
oracle markers live only on the pure-library deep/fast suites (engine,
resolutions_cs, modules, invariants, families, batch, hochschild), which collect
with the base deps alone -- so the audited counts do NOT depend on the [web]/[hpc]
extras being present (those tiers are UNMARKED). The `qpa` count is the QPA oracle
class; qpa tests collect without GAP (module-level skipif), as test_markers already
relies on.

Marked `deep` so it shells out once on the deep leg, not on all 12 fast cells."""
import pathlib
import re
import subprocess
import sys

import pytest

pytestmark = pytest.mark.deep

# parents[2] = repo root (tests/release/<file>). NEVER hardcode a laptop path:
# CI checks out under a different root and would FileNotFoundError on cwd=ROOT.
ROOT = pathlib.Path(__file__).resolve().parents[2]
VENV = sys.executable
PAGE = ROOT / "docs" / "verification.md"

# The canonical marker expressions the page MUST document, one row each.
CANON = (
    "oracle_literature",
    "oracle_crossengine",
    "oracle_selfcert",
    "qpa",
    "oracle_literature or oracle_crossengine or oracle_selfcert or qpa",
)


def _ids(expr):
    """Collected node ids under a -m expression."""
    cmd = [VENV, "-m", "pytest", "-q", "--collect-only", "-p", "no:cacheprovider",
           "-m", expr]
    out = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                         cwd=str(ROOT))
    blob = out.stdout + out.stderr
    assert "PytestUnknownMarkWarning" not in blob, blob
    # returncode 0 == a clean, complete collection with matches. A collection
    # error (e.g. a missing extra) is 2/3; "no tests" is 5 -- either would make
    # the count meaningless, so fail loudly rather than compare a partial number.
    assert out.returncode == 0, (out.returncode, blob)
    return {ln.strip() for ln in out.stdout.splitlines() if "::" in ln}


def _page_counts():
    """Parse {marker-expression: stated-count} out of the audited table.

    A row looks like:  | <label> | `-m oracle_literature` | 342 | <meaning> |
    keyed by the backticked `-m <expr>` cell; the count is the next numeric cell.
    """
    counts = {}
    for line in PAGE.read_text(encoding="utf-8").splitlines():
        if not line.lstrip().startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        expr = None
        n = None
        for c in cells:
            m = re.fullmatch(r"`-m\s+\"?(.+?)\"?`", c)
            if m:
                expr = m.group(1).strip()
            elif re.fullmatch(r"\d+", c):
                n = int(c)
        if expr is not None and n is not None:
            counts[expr] = n
    return counts


def test_page_documents_every_class():
    page = _page_counts()
    assert set(page) == set(CANON), (
        "docs/verification.md oracle-class table must have exactly one row per "
        "canonical marker expression.\n"
        f"  missing: {sorted(set(CANON) - set(page))}\n"
        f"  extra:   {sorted(set(page) - set(CANON))}"
    )


def test_page_counts_match_live_collection():
    page = _page_counts()
    for expr in CANON:
        stated = page[expr]
        live = len(_ids(expr))
        assert live == stated, (
            f"oracle-class count drift for `-m {expr}`: docs/verification.md "
            f"states {stated}, live collection has {live}. Update the table "
            f"(and re-run) -- the page is the audited record."
        )


def test_union_is_bounded_by_the_parts():
    """Sanity: the union count is at least each part and at most their sum
    (classes overlap on purpose, so it is usually strictly less than the sum)."""
    page = _page_counts()
    parts = [page[e] for e in CANON[:4]]
    union = page[CANON[4]]
    assert union >= max(parts), (union, parts)
    assert union <= sum(parts), (union, parts)


def test_no_oracle_markers_in_extras_gated_dirs():
    """Oracle-class markers may live ONLY in always-collectible directories.
    tests/{webapp,gui,hpc} collect only when their extras are installed, so a
    marker there makes the audited class counts environment-dependent -- the
    float-gate CI job ([dev]-only) would see different numbers than a full
    checkout (the 2026-07-27 shakeout; same lesson as the bank-battery
    collection guard). Those directories are contract & infrastructure by the
    Plan-32 ruling; their cross-checks stay unmarked."""
    import re
    root = pathlib.Path(__file__).resolve().parent.parent
    offenders = []
    for d in ("webapp", "gui", "hpc"):
        for f in (root / d).rglob("test_*.py"):
            src = f.read_text(encoding="utf-8")
            if re.search(r"pytest\.mark\.oracle_(literature|crossengine|selfcert)", src):
                offenders.append(str(f.relative_to(root)))
    assert not offenders, (
        "oracle-class markers inside extras-gated test dirs (breaks the audited "
        "count invariant): %s" % offenders)
