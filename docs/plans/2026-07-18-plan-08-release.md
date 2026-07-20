# quiverlab Plan 08 — QPA Extra, CI, Docs Site, PyPI Release, JOSS Paper (release)

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax. Do every cycle in order: write the failing check, run it, see the stated FAIL, write the artifact, run it, see the stated PASS, commit. This plan builds **release infrastructure** (CI, packaging, docs site, JOSS paper, an optional GAP cross-check backend), not new mathematics — but the same discipline holds: every workflow/config is validated by an executable check before it is committed, and the full test suite stays green at every commit.

**Goal:** turn the finished quiverlab library (Plans 01–07 delivered) into a **release candidate**: `pip install quiverlab` works from a built wheel in a fresh venv on macOS/Linux/Windows; a GitHub Actions matrix proves the suite green across OS × Python 3.10–3.13 on both the numba and pure-Python engine paths; a `mkdocs-material` documentation site (auto API reference, CI-executed tutorial notebooks, the "Under the hood" internals chapters) builds `--strict` and deploys to GitHub Pages; an optional `pip install quiverlab[qpa]` GAP backend provides `A.crosscheck(...)` (an independent QPA recomputation of module Ext and Hochschild dimensions via the enveloping-algebra route) wired into a **CI-only Linux cross-check job**; a JOSS `paper.md`/`paper.bib` draft compiles to PDF in CI; and the front-door `README.md` carries badges and a three-line quickstart. **No package is ever published to real PyPI during this plan** — every step goes up to (and validates) the tag; pushing the tag and clicking publish is Marco's button (§9). The `pyproject.toml` `license = {text = "MIT"}` deprecation deferred from Plan 01 is fixed here (PEP 639 SPDX expression).

**Architecture:** all new artifacts are release-infrastructure, side-by-side with the untouched library:
- `.github/workflows/` — **greenfield** (no `.github/` exists today). Five workflows: `ci.yml` (test matrix + engine-path legs + float-gate/lint), `qpa.yml` (GAP+QPA cross-check, Linux), `docs.yml` (build `--strict` + deploy to Pages), `paper.yml` (JOSS draft PDF), `release.yml` (build + `twine check` + OIDC trusted publishing on tag).
- `src/quiverlab/qpa/` — the `[qpa]` extra: a `libgap` session harness with a **skip-if-absent guard** (`gap_available()`), GAP script builders for the enveloping-algebra HH route and module Ext, and `Algebra.crosscheck(...)` (lazy-imported so the pure core never depends on GAP).
- `mkdocs.yml` + `docs/` site sources: `docs/index.md`, `docs/gen_ref_pages.py` (mkdocstrings auto API reference), tutorials rendered from `docs/tutorials/*.ipynb`, the internals chapters, a Development/Release page. Deployed to `https://marcoarmenta.github.io/quiverlab/`.
- `paper/paper.md` (JOSS) — bibliography reuses the single **packaged** canonical `src/quiverlab/citations/references.bib` (Plan 06, accessed via `quiverlab.citations.references_bib_path()`); no second bib is committed (§ "Bibliography single-source").
- Root: modernized `pyproject.toml`, overhauled `README.md`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `CITATION.cff`, `tests/conftest.py` (marker auto-assignment), `tests/release/` (freshness gate, packaging/paper/docs validation).

**Tech Stack:** Python ≥ 3.10 (venv is 3.12). Build/dev tooling added to `[project.optional-dependencies]`: `build`, `twine` (packaging), `mkdocs-material` + `mkdocstrings[python]` + `mkdocs-jupyter` + `mkdocs-gen-files` + `mkdocs-literate-nav` + `mkdocs-section-index` + `mike` (docs), `pytest` (already present). The `[qpa]` extra depends on `passagemath-gap` (macOS/Linux; loud graceful message on Windows — see Task 4, **VERIFY AT EXECUTION** the exact pin/wheels). No new **hard** runtime dependency beyond what Plans 01–07 already require (`numpy`, `sympy`, `matplotlib`); `numba` stays the optional `[fast]` extra; `passagemath-gap` is the optional `[qpa]` extra. Exact arithmetic only; the float-ban AST gate stays green.

---

## Global Constraints

- **Repo root:** `/Users/marco/Desktop/HomologicalNetworks/quiverlab`. All paths below are relative to it. **GitHub owner/repo:** `MarcoArmenta/quiverlab` (from `git remote`). **Docs (Pages) URL:** `https://marcoarmenta.github.io/quiverlab/` — this exact string is the single source for `mkdocs.yml` `site_url`, the README docs badge/link, the JOSS repo link, and the Plan-09 web config `QLWEB_DOCS_URL` (§ "Cross-plan contracts"). If Marco publishes under a different account/handle, change it in exactly these places (a consistency test in Task 8 pins that they agree).
- **Interpreter:** use the project venv **`/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python`** (Python 3.12). The system `python`/`python3` is older and MUST NOT be used.
- **Thread throttle:** prefix **every** test command with `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2` (Marco's machine has crashed under agent fleets; keep thread/memory pressure low). CI workflows set the same env on every test step.
- **Suite discipline (binding).** Local iteration runs a **focused foreground** command (the single file/marker relevant to the task). Before each commit run the **full suite in the foreground**. Once, at the acceptance task, run **ONE tracked background full suite**, await it, then commit + report. All test commands run from the repo root:
  `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q`
- **Green at every commit.** The full suite (`pytest -q`, i.e. `-m "fast or deep"` after Task 2, excluding `qpa` which skips without GAP) passes before every commit. The float-ban AST gate `tests/test_no_floats.py` MUST stay green at every commit — this plan adds **no** `src/` float/complex literals (the `qpa` package passes GAP scripts as strings and reads exact integers back; no floats).
- **AST gate green.** `tests/test_no_floats.py` is part of every task's suite.
- **Banks are READ-ONLY.** `/Users/marco/Desktop/HomologicalNetworks/HomologicalAlgebra/` and every other bank are never read or written in this plan. Plan 08 touches only `src/quiverlab/`, tests, docs, and root release files already in the repo.
- **No real PyPI upload during this plan.** Everything is built, `twine check`ed, and dry-run installed into a fresh venv locally; the `release.yml` workflow is written, YAML-validated, and committed, but it triggers only on a pushed `v*` tag and publishes through PyPI **trusted publishing** — the tag push and the PyPI "create pending publisher" step are **Marco's** actions, done after this plan lands. No API tokens are committed anywhere.
- **GAP is not installed locally.** Every GAP/QPA interaction is fully scripted (GAP code as strings) and guarded by `gap_available()`: `qpa`-marked tests **skip** locally and are **mandatory** in the `qpa.yml` CI job (which sets `QUIVERLAB_REQUIRE_QPA=1`, turning the skip into a hard failure). The one worked GAP fixture (Task 5, HH of `kA_2` over `GF(2)`) has its expected output written out and flagged **VERIFY AT EXECUTION** (it cannot run on this machine).
- **Commits:** conventional prefixes (`feat:`, `test:`, `chore:`, `docs:`, `ci:`, `build:`); every commit message ends with the **two-line** trailer the repo uses (verified in `git log`):
  ```
  Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd
  ```
  The `Claude-Session` line carries the **executing** session's `claude.ai/code` URL (the harness substitutes it per its Git convention; the URL shown is the one in scope now). The per-task commit blocks below inline both lines.
- **Prerequisite plans.** Plan 08 releases the *whole* library; it assumes Plans 03–07 delivered (general `kQ/I`, Chouhy–Solotar, modules/invariants, families/`bibliography()`/`references.bib`, viz/trace). **Task 1 is an executable freshness gate that STOPS the plan if any prerequisite surface is missing** — this plan is not runnable against a half-built library, and the gate says exactly which plan is incomplete.

---

## Cross-plan contracts (state once, honored throughout)

- **Docs URL.** `mkdocs.yml:site_url = https://marcoarmenta.github.io/quiverlab/`. The Plan-09 web app reads this as `QLWEB_DOCS_URL` (web-spec §3/§16, the `/literature` page and docs links). Plan 08 owns the value; Plan 09 consumes it. Task 8's consistency test asserts the README docs link, the `CITATION.cff` url, and `mkdocs.yml:site_url` are byte-identical, so there is one string to change.
- **Bibliography single-source.** Canonical file: the **packaged** **`src/quiverlab/citations/references.bib`** (a Plan-06 deliverable, resolved at runtime by `quiverlab.citations.references_bib_path()`; also surfaced by `quiverlab.bibliography()`). It ships **inside the installed package** (Plan-06 `pyproject` package-data), so every consumer reads it by import, never by repo path. The JOSS paper does **not** commit a second `.bib`: `paper/paper.md` front matter sets `bibliography: paper.bib`, and the `paper.yml` job **copies** `src/quiverlab/citations/references.bib → paper/paper.bib` at build time (`paper/paper.bib` is `.gitignore`d). A freshness test (Task 9) asserts every `@key` cited in `paper.md` resolves in the packaged bib (read via `references_bib_path()`). The **docs site** renders a `References` page generated at build time from the same packaged file (Task 8, a `mkdocs-gen-files` hook — no repo copy). Plan 09's `/literature` page renders it via `bibliography()`. One packaged bib, three consumers (paper, docs, web), no duplication and no `docs/`-tree copy.
- **Marker vocabulary.** `fast` / `deep` / `qpa` / `slow` (Task 2). Every CI cell's `pytest -m` expression is listed in Task 6; the web-app CI (Plan 09) reuses the same conftest auto-marking for its `tests/webapp/` (marked `fast`).

---

## Task overview (13 tasks)

1. Executable **freshness gate** — STOP on drift (prerequisite surfaces, pyproject state, docs sources, marker state).
2. **pytest markers + `conftest.py`** — `fast`/`deep`/`qpa`/`slow` auto-assignment (the CI-split foundation).
3. **`pyproject.toml` modernization** — SPDX license fix (the Plan-01 deferral), classifiers, urls, `[qpa]`/`[docs]` extras, build-system bump, version single-source test.
4. **`[qpa]` extra** — `quiverlab.qpa` (libgap session + `gap_available()` guard), GAP script builders, `Algebra.crosscheck(...)`.
5. **QPA cross-check tests** — `tests/qpa/` (marked `qpa`), the worked `kA_2`/`GF(2)` enveloping-algebra HH GAP fixture + module Ext cross-check.
6. **GitHub Actions CI** — `ci.yml`: fast matrix (OS × py3.10–3.13), deep numba + pure legs, float-gate/lint.
7. **QPA cross-check CI job** — `qpa.yml`: install GAP+QPA on Linux, mandatory `-m qpa`.
8. **Docs site** — `mkdocs.yml` + mkdocstrings auto API reference + CI-executed tutorials + internals chapters + `docs.yml` Pages deploy; versioned-docs decision.
9. **JOSS paper** — `paper/paper.md` full draft + bib single-source + `paper.yml` draft-PDF workflow.
10. **PyPI packaging + release** — `release.yml` (trusted publishing on tag), name-availability check, fresh-venv install dry-run (no real upload).
11. **README overhaul** — badges + three-line quickstart + links (docs, web GUI, tutorials, internals, JOSS).
12. **Community + citation files** — `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `CITATION.cff`; docs Development/Release page; internals README honest-coverage update; the internals-chapter decision (skip `12-release-infrastructure.md`, justified).
13. **Release-candidate acceptance** — the full checklist; one tracked background full suite awaited; docs build live; paper compiles; fresh-venv install; commit + report.

---

### Task 1: Executable freshness gate — STOP on drift

**Files:**
- Create: `tests/release/__init__.py` (empty), `tests/release/test_freshness.py`, `scripts/release_freshness.py`

**Interfaces:**
- Consumes: the public library surface (`quiverlab`), the repo tree.
- Produces: a hard STOP if a prerequisite is missing. This task writes **no** release infrastructure — it verifies the ground is solid before Tasks 2–13 build on it. **If the gate fails, STOP the plan and report which prerequisite plan is incomplete.** (Reconciled to merged reality: Plans 03–07 are all merged to `main`, so the gate now exits **0** — it prints an informational NOTE that the license is not yet SPDX-fixed, which is a Task-3 TODO, not a STOP. The task is therefore near-vacuous on current `main`; its load-bearing content is the **exit-code contract**: STOP only on genuine prerequisite drift, never on the deferred-license line.)

**Exit-code reconciliation (binding).** The gate separates two concerns: (1) genuine **prerequisite** drift (Plans 03–07 surfaces + docs sources) — the only thing that returns exit 1; and (2) the deprecated `license = { text = "MIT" }` table, which is a Task-3 fix-it **TODO** surfaced by `license_todo()` as an **informational** note that does NOT force exit 1. Without this split the gate would exit 1 on current `main` (license not SPDX-fixed until Task 3), contradicting Step 2's "exit=0" and Task 13's "release_freshness.py exits 0" checklist item. With it, Step 2's exit=0 holds today and Task 3 still fixes the license.

- [ ] **Step 1: Write the gate (test + script share one checker)**

`scripts/release_freshness.py`:

```python
"""Plan 08 freshness gate: refuse to build release infrastructure against a
half-built library. Prints a report; exits nonzero (STOP) on any drift.

Run:  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 \
      /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python scripts/release_freshness.py
"""
from __future__ import annotations

import importlib
import importlib.util
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent


def check() -> list[str]:
    """Return a list of PREREQUISITE drift messages (Plans 03-07 library surfaces +
    docs sources); empty list == fresh. The deprecated-license line is NOT included
    here -- it is an informational Task-3 TODO returned by license_todo() -- so a
    not-yet-SPDX license never forces a STOP."""
    drift: list[str] = []

    # --- prerequisite library surfaces (Plans 03-07) -----------------------
    import quiverlab
    from quiverlab import Quiver, CC

    # Plan 03: general (non-monomial) kQ/I lowers, groebner.system present.
    try:
        importlib.import_module("quiverlab.groebner.system")
        A = Quiver([1, 2, 3, 4],
                   {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
                   ).algebra(relations=["a*b - c*d"], field=CC)
        assert A.dim == 9
    except Exception as e:  # noqa: BLE001
        drift.append(f"Plan 03 (general kQ/I / groebner.system) incomplete: {e!r}")

    # Plan 04: Chouhy-Solotar resolution package.
    if importlib.util.find_spec("quiverlab.resolutions_cs") is None:
        drift.append("Plan 04 (quiverlab.resolutions_cs) not present")

    # Plan 05: modules + Ext + gl.dim on the public Algebra.
    for name in ("ext", "simple", "projective", "global_dimension"):
        if not hasattr(quiverlab.Algebra, name):
            drift.append(f"Plan 05 surface Algebra.{name} missing")

    # Plan 06: families catalog + bibliography() + the PACKAGED citations registry.
    for name in ("zoo", "families", "bibliography"):
        if not hasattr(quiverlab, name):
            drift.append(f"Plan 06 surface quiverlab.{name}() missing")
    try:
        from quiverlab import citations
        if not citations.references_bib_path().is_file():
            drift.append("Plan 06 packaged citations.references_bib_path() does not "
                         "resolve to a file (JOSS paper.bib + docs + web /literature depend on it)")
    except Exception as e:  # noqa: BLE001
        drift.append(f"Plan 06 quiverlab.citations registry missing: {e!r}")

    # Marker policy: heavy test dirs must be in the deep bucket. Any NEW top-level
    # tests/ dir not covered by tests/conftest.py buckets is triaged loudly by
    # tests/release/test_markers.py::test_buckets_partition_the_suite; this note is
    # the standing reminder to add Plan 05/06/07 heavy dirs to _DEEP_DIRS.

    # Plan 07: viz + trace on the public Algebra.
    for name in ("draw", "tikz"):
        if not hasattr(quiverlab.Algebra, name):
            drift.append(f"Plan 07 surface Algebra.{name} missing")

    # --- docs sources this plan will publish -------------------------------
    tut = sorted((ROOT / "docs" / "tutorials").glob("*.ipynb"))
    if len(tut) < 3:
        drift.append(f"expected >=3 tutorial notebooks, found {len(tut)}")
    internals = sorted((ROOT / "docs" / "internals").glob("[0-9][0-9]-*.md"))
    if len(internals) < 7:
        drift.append(f"expected >=7 internals chapters, found {len(internals)}")

    # The deprecated license={text=...} table is NOT prerequisite drift -- it is a
    # Task-3 fix-it TODO, reported separately by license_todo() as an INFORMATIONAL
    # note and never a STOP. Keeping it out of `drift` is what lets this gate exit 0
    # on current `main` (where the license is not SPDX-fixed until Task 3 runs), so
    # Step 2's "exit=0" and Task 13's "release_freshness.py exits 0" both hold.
    return drift


def license_todo() -> list[str]:
    """Informational only (never a STOP): the deprecated PEP-639 license={text=...}
    table that Task 3 replaces with the SPDX `license = "MIT"` string. main() prints
    these notes but they do NOT affect the exit code."""
    notes: list[str] = []
    pp = (ROOT / "pyproject.toml").read_text()
    if 'license = { text = "MIT" }' in pp or 'license = {text = "MIT"}' in pp:
        notes.append("pyproject still uses the deprecated license={text=...} table "
                     "(Task 3 replaces it with the PEP 639 SPDX string)")
    return notes


def main() -> int:
    drift = check()
    for note in license_todo():
        print(f"PLAN 08 FRESHNESS GATE: NOTE (informational, not a STOP) -- {note}")
    if drift:
        print("PLAN 08 FRESHNESS GATE: STOP -- prerequisites drifted:\n")
        for d in drift:
            print(f"  - {d}")
        print("\nComplete the named prerequisite plan(s) before running Plan 08.")
        return 1
    print("PLAN 08 FRESHNESS GATE: OK -- library surfaces, docs sources, and "
          "pyproject baseline are as expected.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

`tests/release/test_freshness.py`:

```python
"""Plan 08 freshness gate as a test. `check()` returns only PREREQUISITE drift
(Plans 03-07 surfaces + docs sources); the deprecated-license line is an
informational Task-3 TODO handled by `license_todo()` and never a STOP, so this
test asserts prerequisite freshness directly and is agnostic to whether the license
has been SPDX-fixed yet."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent.parent / "scripts"))
import release_freshness as rf  # noqa: E402


def test_prerequisite_surfaces_present():
    drift = rf.check()
    assert drift == [], "Plan 08 prerequisites drifted:\n" + "\n".join(drift)
```

- [ ] **Step 2: Run the gate**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python scripts/release_freshness.py; echo "exit=$?"`
Expected (Plans 03–07 merged to `main`): `PLAN 08 FRESHNESS GATE: OK`, `exit=0`. On current `main` it also prints one informational line — `NOTE (informational, not a STOP) -- pyproject still uses the deprecated license={text=...} table` — which Task 3 clears; this NOTE does **not** change the exit code (still 0). **If it prints STOP (exit 1), halt the plan** and report the listed missing prerequisites to Marco — Plan 08 cannot proceed until they land.

- [ ] **Step 3: (no implementation — this is the gate)**

- [ ] **Step 4: Run the focused suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/release/test_freshness.py tests/test_no_floats.py -q`
Expected: pass (gate green) — only if prerequisites are present. Otherwise STOP as in Step 2.

- [ ] **Step 5: Commit**

```bash
git add scripts/release_freshness.py tests/release/
git commit -m "test(release): Plan 08 freshness gate (STOP on prerequisite drift)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
```

---

### Task 2: pytest markers + `conftest.py` — the CI-split foundation

**Files:**
- Create: `tests/release/test_markers.py`
- Modify: `tests/conftest.py` (**MERGE** the bucket hook into the existing Plan-07 conftest — do NOT overwrite it), `pyproject.toml` (`[tool.pytest.ini_options].markers`)

> **Load-bearing (do NOT regress `main`):** Plan 07 already shipped `tests/conftest.py` with a `pytest_configure` that registers the `verbose_default` marker **and** an autouse `_quiet_traces` fixture that forces `quiverlab.verbose = False` for the whole suite. **The Plan-07 `_quiet_traces` autouse fixture + `verbose_default` marker registration MUST be preserved** — overwriting the file would delete them, breaking 5+ trace tests, un-quieting the suite (stray `./quiverlab_traces/` writes), and re-introducing a `PytestUnknownMarkWarning` that fails this task's own `test_markers.py`. This task therefore **modifies** the existing conftest: it ADDS `pytest_collection_modifyitems` (and its `_bucket`/`_DEEP_*` helpers) while KEEPING the existing `pytest_configure` and `_quiet_traces` verbatim. The CI-split markers (`fast`/`deep`/`qpa`/`slow`) are registered in `pyproject.toml` (below), so `pytest_configure` only needs to keep registering `verbose_default`; if a future revision registers markers via `pytest_configure` instead, it must APPEND to the existing body (register `verbose_default` **and** the new markers), never replace it.

**Interfaces:**
- Consumes: the collected test tree.
- Produces: three **mutually-exclusive, exhaustive** CI buckets — `fast` (every matrix cell), `deep` (heavy suites; one Linux leg), `qpa` (needs `[qpa]`; the CI QPA job only) — plus one **orthogonal sub-tag** `slow`. Policy (binding): **every test is in exactly one bucket** (enforced by `test_buckets_partition_the_suite`), and **`slow` implies `deep`** (a slow test always rides the deep leg — it is never in the fast matrix). Explicit `@pytest.mark.fast/deep/qpa` on a test wins over the path default; a bare `@pytest.mark.slow` is auto-promoted to `deep`.

**Marker plan (exact).** The suite has ~615 tests today (Plan 03 state) and grows through Plans 04–07. Heavy time lives in `tests/engine/` (32 files), the Chouhy–Solotar suite (Plan 04), deep groebner completion, module resolutions (Plan 05), and the families/batch scans (Plan 06). `conftest.py` assigns buckets by top-level test directory / filename, with per-test overrides:

- `tests/qpa/**` → `qpa`
- `tests/{engine,resolutions_cs,modules,families,batch}/**` and named heavy files (`test_complete.py`, `test_deepen.py`, `test_properties.py`, `test_acceptance.py`, `test_cs_*`, `test_bardzell*`, `test_minimal*`) → `deep`
- everything else (incl. `tests/viz/**` — small drawings, worth cross-platform coverage) → `fast`
- any test explicitly marked `slow` (with no bucket) → promoted to `deep`

**Future-heavy-dir policy:** when Plans 05–07 add heavy top-level dirs, add them to `_DEEP_DIRS`. If one is forgotten, `test_buckets_partition_the_suite` still passes (it only checks disjoint+exhaustive, not wall-time) but its docstring names the triage; the freshness-gate note (Task 1) is the standing reminder. A genuinely slow test that slips into `fast` is caught by the fast-matrix wall-time, and can be pinned with `@pytest.mark.slow` (→ deep).

`tests/conftest.py` — **MERGE** the following into the EXISTING Plan-07 conftest. Keep its module docstring, `import quiverlab`, the `pytest_configure` that registers `verbose_default`, and the autouse `_quiet_traces` fixture EXACTLY as shipped; only ADD the bucket helpers + `pytest_collection_modifyitems`. The complete merged file (Plan-07 parts preserved verbatim, Plan-08 additions marked) is:

```python
"""Global test fixtures + CI-bucket auto-assignment.

Plan 07: the worked-steps trace subsystem is ON by default (quiverlab.verbose =
True, spec D9); force it OFF for the whole suite so that unrelated computations do
not write ./quiverlab_traces/ files. Trace tests that need it opt back in explicitly
(set quiverlab.verbose = True or pass verbose=True / trace=[...]).

Plan 08 Task 2: auto-assign CI buckets by test location so the suite can be split
across the GitHub Actions matrix. Explicit bucket markers on a test win.
Buckets (exactly one per test; disjoint + exhaustive, enforced by the partition test):
  qpa  -- needs the [qpa] extra (passagemath-gap + QPA); CI QPA job only.
  deep -- heavy engine / resolution / module / families / batch suites; one Linux leg.
  fast -- everything else; runs on every OS x Python matrix cell.
Orthogonal sub-tag:
  slow -- an individually long test; IMPLIES deep (never runs in the fast matrix).
"""
import pytest

import quiverlab


# --- Plan 07 (PRESERVE VERBATIM) -------------------------------------------
def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "verbose_default: observe the shipped quiverlab.verbose default "
        "(opt out of the quiet-suite fixture)",
    )


@pytest.fixture(autouse=True)
def _quiet_traces(request):
    # Always capture + restore, so even a `verbose_default`-marked test that mutates
    # quiverlab.verbose cannot leak state into the next test (review hardening).
    prev = getattr(quiverlab, "verbose", True)
    if request.node.get_closest_marker("verbose_default") is None:
        # Force the trace subsystem OFF so unrelated tests don't write ./quiverlab_traces/.
        quiverlab.verbose = False
    # else: leave the shipped default (quiverlab.verbose = True) observable.
    try:
        yield
    finally:
        quiverlab.verbose = prev


# --- Plan 08 Task 2 (ADD) --------------------------------------------------
# Top-level dirs (relative to tests/) whose tests are heavy -> the deep leg.
# (resolutions_cs is the Chouhy-Solotar suite; there is NO tests/chouhy_solotar dir.)
_DEEP_DIRS = ("engine", "resolutions_cs", "modules", "families", "batch")
# Individually heavy files that may live outside the deep dirs.
_DEEP_FILES = ("test_complete.py", "test_deepen.py", "test_properties.py",
               "test_acceptance.py", "test_cs_", "test_bardzell", "test_minimal")
_BUCKETS = {"fast", "deep", "qpa"}


def _bucket(nodeid: str) -> str:
    # nodeid looks like "tests/engine/test_foo.py::test_bar"
    parts = nodeid.replace("\\", "/").split("/")
    top = parts[1] if len(parts) > 1 else ""
    fname = parts[-1].split("::")[0]
    if top == "qpa":
        return "qpa"
    if top in _DEEP_DIRS or any(fname.startswith(f) or fname == f for f in _DEEP_FILES):
        return "deep"
    return "fast"


def pytest_collection_modifyitems(config, items):
    for item in items:
        names = {m.name for m in item.iter_markers()}
        # `slow` always rides the deep leg (unless it is a qpa test).
        if "slow" in names and not (names & {"deep", "qpa"}):
            item.add_marker(pytest.mark.deep)
            names.add("deep")
        if names & _BUCKETS:          # explicit bucket (fast/deep/qpa) wins -> no auto-mark
            continue
        item.add_marker(getattr(pytest.mark, _bucket(item.nodeid)))
```

Register the markers in `pyproject.toml` (replace the existing `[tool.pytest.ini_options]` block):

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
markers = [
    "fast: quick tests; run on every CI matrix cell (OS x Python)",
    "deep: heavy engine/resolution/Chouhy-Solotar suites; one Linux cell (~19 min)",
    "qpa: requires the [qpa] extra (passagemath-gap + QPA); the CI QPA job only",
    "slow: individually long tests; opt-in via -m slow",
]
```

- [ ] **Step 1: Write the failing test**

`tests/release/test_markers.py`:

```python
"""The CI buckets fast/deep/qpa are a PARTITION of the collected suite: pairwise
disjoint and exhaustive; slow implies deep; no unknown-marker warnings (Task 2).

Marked `deep` so it runs once on the deep leg, not on all 12 fast matrix cells
(it shells out four collections)."""
import subprocess

import pytest

pytestmark = pytest.mark.deep
ROOT = "/Users/marco/Desktop/HomologicalNetworks/quiverlab"
VENV = f"{ROOT}/.venv/bin/python"


def _ids(expr):
    """Collected node ids under a -m expression ('' == the whole suite)."""
    cmd = [VENV, "-m", "pytest", "-qq", "--collect-only", "-p", "no:cacheprovider"]
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
```

- [ ] **Step 2: Run to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/release/test_markers.py -q`
Expected: FAIL — markers not registered yet (`PytestUnknownMarkWarning`) / `conftest.py` absent.

- [ ] **Step 3: MERGE the bucket hook into the existing `tests/conftest.py` and register markers in `pyproject.toml`** (both above). Preserve Plan 07's `pytest_configure` (`verbose_default`) and the `_quiet_traces` autouse fixture verbatim; only add the bucket helpers + `pytest_collection_modifyitems`. Do NOT overwrite the file.

- [ ] **Step 4: Verify the split by hand + run the suite**

Run:
```bash
cd /Users/marco/Desktop/HomologicalNetworks/quiverlab
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q --collect-only -m fast | tail -3
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q --collect-only -m deep | tail -3
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m fast
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/release/test_markers.py tests/test_no_floats.py -q
```
Expected: `-m fast` and `-m deep` partition the non-qpa suite (their collected counts sum to the total); `-m fast` runs green quickly; the marker test and float gate pass. No `PytestUnknownMarkWarning`.

- [ ] **Step 5: Commit**

```bash
git add tests/conftest.py tests/release/test_markers.py pyproject.toml
git commit -m "test(ci): fast/deep/qpa/slow markers auto-assigned by path (CI-split foundation)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
```

---

### Task 3: `pyproject.toml` modernization — SPDX license fix, classifiers, urls, extras

**Files:**
- Modify: `pyproject.toml`
- Create: `tests/release/test_packaging_metadata.py`

**Interfaces:**
- Consumes: `pyproject.toml`, `quiverlab.__version__`.
- Produces: a PEP 639-compliant, non-deprecated build configuration; the `[qpa]`/`[docs]` extras; project URLs; and a version single-source consistency test. **Fixes the Plan-01-deferred `license = {text = "MIT"}` deprecation.**

**The deprecation.** setuptools ≥ 77 supports the PEP 639 SPDX **license expression** `license = "MIT"` (a string, not a table) plus `license-files`, and **warns** on the legacy `license = {text = "MIT"}` table and on the `License :: OSI Approved :: MIT License` **classifier** when a SPDX expression is present. Fix: switch to the string form, add `license-files = ["LICENSE"]`, and **drop** the `License ::` classifier. (**VERIFY AT EXECUTION**: confirm the installed setuptools is ≥ 77; if the build pins an older setuptools, bump `build-system.requires`.)

Full modernized `pyproject.toml` (complete file):

```toml
[build-system]
requires = ["setuptools>=77"]
build-backend = "setuptools.build_meta"

[project]
name = "quiverlab"
version = "0.1.0.dev0"
description = "Quivers with relations and Hochschild theory, exactly, for algebraists"
readme = "README.md"
requires-python = ">=3.10"
license = "MIT"
license-files = ["LICENSE"]
authors = [{ name = "Marco Armenta", email = "drmarcoarmenta@gmail.com" }]
keywords = [
    "quiver", "path algebra", "Hochschild cohomology", "Gerstenhaber bracket",
    "representation theory", "homological algebra", "exact arithmetic",
]
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Mathematics",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Programming Language :: Python :: 3.13",
    "Operating System :: OS Independent",
    "Typing :: Typed",
]
dependencies = ["numpy>=1.21", "sympy>=1.12", "matplotlib>=3.7"]

[project.optional-dependencies]
fast = ["numba>=0.64"]
# [qpa]: passagemath-gap 10.8.x ships prebuilt GAP + (via its own [qpa] sub-extra)
# QPA v1.37 + GBNP, in manylinux/musllinux + macOS wheels (Python 3.11-3.14). No
# Windows wheel -> the `sys_platform != 'win32'` marker makes the extra installable-
# but-empty on Windows, and quiverlab.qpa prints the loud graceful message at runtime
# (Task 4). VERIFY the exact pin/extra name at execution (research: 10.8.6, 2026-07-16).
qpa = ["passagemath-gap[qpa]>=10.8; sys_platform != 'win32'"]
docs = [
    "mkdocs-material>=9.5",
    "mkdocstrings[python]>=0.26",
    "mkdocs-jupyter>=0.25",
    "mkdocs-gen-files>=0.5",
    "mkdocs-literate-nav>=0.6",
    "mkdocs-section-index>=0.3",
    "mike>=2.1",
]
dev = ["pytest>=8", "build>=1.2", "twine>=5.1"]

[project.urls]
Homepage = "https://marcoarmenta.github.io/quiverlab/"
Documentation = "https://marcoarmenta.github.io/quiverlab/"
Repository = "https://github.com/MarcoArmenta/quiverlab"
Issues = "https://github.com/MarcoArmenta/quiverlab/issues"
Changelog = "https://github.com/MarcoArmenta/quiverlab/blob/main/CHANGELOG.md"
# "Web GUI" and the JOSS DOI are added at web-deploy / JOSS-acceptance time.

[tool.setuptools.packages.find]
where = ["src"]

# LOAD-BEARING -- do NOT drop. citations/registry.py loads references.bib and
# families/zoo.py loads zoo_catalog.json by __file__-relative path; without this
# block both files are absent from the built wheel/sdist and bibliography()/zoo()
# break in an installed environment (Plan 06). See the sdist/wheel data-file
# packaging check in Task 10 (test_build.py, open item #4).
[tool.setuptools.package-data]
"quiverlab.citations" = ["references.bib"]
"quiverlab.families" = ["zoo_catalog.json"]

[tool.pytest.ini_options]
testpaths = ["tests"]
markers = [
    "fast: quick tests; run on every CI matrix cell (OS x Python)",
    "deep: heavy engine/resolution/Chouhy-Solotar suites; one Linux cell (~19 min)",
    "qpa: requires the [qpa] extra (passagemath-gap + QPA); the CI QPA job only",
    "slow: individually long tests; opt-in via -m slow",
]
```

Notes:
- **`package-data` + `matplotlib>=3.7` are load-bearing (Plan 06 bib/zoo, Plan 07 viz) — never drop or downgrade them.** This modernization PRESERVES the merged-`main` `dependencies` (incl. `matplotlib>=3.7`) and the entire `[tool.setuptools.package-data]` block verbatim; it only ADDS the SPDX license, classifiers, urls, and extras. Dropping the package-data block regresses `main` (built wheels would ship without `references.bib`/`zoo_catalog.json`); downgrading matplotlib below `3.7` regresses the Plan-07 viz floor.
- **`matplotlib>=3.7`** is a **hard** dependency by Plan 07 (viz, spec §9); kept identical to merged `main`.
- **`numba>=0.64`** — kept identical to merged `main` (the venv runs numba 0.66, which satisfies `>=0.64`); this modernization deliberately does NOT change the numba floor.
- `Development Status :: 4 - Beta` matches semver 0.x battle-testing (§9); bumped to `5 - Production/Stable` at 1.0 / JOSS acceptance.

`tests/release/test_packaging_metadata.py`:

```python
"""pyproject is PEP 639-modern (no deprecated license table / classifier), the
version is single-sourced, and the extras/urls are present (Plan 08 Task 3)."""
import pathlib
import tomllib

import quiverlab

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
PP = tomllib.loads((ROOT / "pyproject.toml").read_text())


def test_license_is_spdx_expression_not_table():
    lic = PP["project"]["license"]
    assert lic == "MIT", "use the PEP 639 SPDX string, not {text=...}"
    assert PP["project"]["license-files"] == ["LICENSE"]


def test_no_deprecated_license_classifier():
    bad = [c for c in PP["project"]["classifiers"] if c.startswith("License ::")]
    assert bad == [], f"drop the deprecated License classifier(s): {bad}"


def test_version_single_source():
    assert PP["project"]["version"] == quiverlab.__version__


def test_extras_and_urls_present():
    extras = PP["project"]["optional-dependencies"]
    assert {"fast", "qpa", "docs", "dev"} <= set(extras)
    assert any("passagemath-gap" in d for d in extras["qpa"])
    urls = PP["project"]["urls"]
    assert urls["Documentation"] == "https://marcoarmenta.github.io/quiverlab/"
    assert urls["Repository"] == "https://github.com/MarcoArmenta/quiverlab"
```

- [ ] **Step 1: Write the failing test** (above).
- [ ] **Step 2: Run to verify it fails**
  Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/release/test_packaging_metadata.py -q`
  Expected: FAIL — current `license = { text = "MIT" }`, no extras/urls/classifiers.
- [ ] **Step 3: Rewrite `pyproject.toml`** (complete file above).
- [ ] **Step 4: Verify build config + suite**
  Run:
  ```bash
  cd /Users/marco/Desktop/HomologicalNetworks/quiverlab
  .venv/bin/python -m pip install -q build >/dev/null 2>&1 || true
  .venv/bin/python -m build --sdist --wheel 2>&1 | tail -5   # expect no SetuptoolsDeprecationWarning about license
  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/release/test_packaging_metadata.py tests/test_no_floats.py -q
  ```
  Expected: the build emits **no** license/classifier deprecation warning; the metadata test and float gate pass. (If `build`/`setuptools>=77` is not yet installed, `pip install build` first; the release env installs `.[dev]`.)
  **Env note:** the `python -m build` verify needs `build` in the venv — `/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pip install build` (it is in the `[dev]` extra this task adds; install it into `.venv` if not already present).
- [ ] **Step 5: Commit**
  ```bash
  git add pyproject.toml tests/release/test_packaging_metadata.py
  git commit -m "build: modernize pyproject (PEP 639 SPDX license, classifiers, urls, qpa/docs extras)

Fixes the Plan-01-deferred license={text} deprecation; adds [qpa]/[docs] extras
and project URLs.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
  ```

---

### Task 4: `[qpa]` extra — `quiverlab.qpa` (libgap session + guard) + `A.crosscheck(...)`

**Files:**
- Create: `src/quiverlab/qpa/__init__.py`, `src/quiverlab/qpa/session.py`, `src/quiverlab/qpa/scripts.py`, `src/quiverlab/qpa/crosscheck.py`
- Modify: `src/quiverlab/errors.py` (add `QpaUnavailableError`), `src/quiverlab/__init__.py` (export it), `src/quiverlab/core/algebra.py` (add `Algebra.crosscheck`, lazy-imported)
- Create: `tests/qpa/__init__.py` (empty), `tests/qpa/test_session_guard.py`

**Interfaces:**
- Consumes: `Algebra`, `Quiver`, `Domain` (via the algebra), `libgap` from `passagemath-gap` (optional).
- Produces: a GAP-session harness (`gap_available`, `require_gap`, `run`), GAP script builders (path algebra `kQ/I` in QPA; the enveloping-algebra HH route; module Ext), and `Algebra.crosscheck(what, ...)` — an **independent** QPA recomputation used for validation. **The pure core never imports GAP**; `Algebra.crosscheck` lazy-imports `quiverlab.qpa` and raises `QpaUnavailableError` (with a fix-it) if the extra is absent.

**Research ground truth (2026, sourced; some flagged VERIFY AT EXECUTION):**
- PyPI distribution `passagemath-gap` (latest **10.8.6**, 2026-07-16); Python **3.11–3.14**; wheels for **Linux** (manylinux+musllinux x86_64/aarch64) and **macOS** (arm64 + x86_64, macOS 13+); **no native Windows wheel** (WSL2 or conda-forge `gap-core`/`gap-pkg-qpa`).
- QPA (Quivers and Path Algebras) v1.37 (2026-05-13); pulled by the `passagemath-gap[qpa]` sub-extra (fallback `[full]`), which also resolves its dependency GBNP. Load with `libgap.LoadPackage("qpa")`.
- Import: `from sage.libs.gap.libgap import libgap` (robust) or `from passagemath_gap import libgap`.
- QPA has **no packaged Hochschild cohomology** (spec §2): assemble `HH^n(A) = Ext^n_{A^e}(A, A)`, `A^e = EnvelopingAlgebra(A)`. Confirmed QPA symbols: `Quiver`, `PathAlgebra`, `QuotientOfPathAlgebra` / `kQ/rels`, `EnvelopingAlgebra`, `TensorProductOfAlgebras`, `RightModuleOverPathAlgebra`, `DimensionVector`, `Dimension`, `ProjectiveResolutionOfPathAlgebraModule(M, n)`, `ExtOverAlgebra(M, N)`, `ExtAlgebraGenerators(M, n)`, `1stSyzygy`, `NthSyzygy`. **VERIFY AT EXECUTION:** the exact QPA operation that presents `A` as a right module over `EnvelopingAlgebra(A)` (the object fed to the resolution) — QPA has bimodule/enveloping-module machinery; pin the constructor name from QPA manual ch. 6/8 before relying on it.

`src/quiverlab/errors.py` — append:

```python
class QpaUnavailableError(QuiverlabError):
    """The optional [qpa] GAP backend is not available (not installed, wrong
    platform, or QPA failed to load). The pure-Python core does not need it."""
```

`src/quiverlab/qpa/session.py`:

```python
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
```

`src/quiverlab/qpa/scripts.py`:

```python
"""GAP/QPA script builders. Translate a quiverlab Algebra's presentation (its
quiver + relations, over QQ or GF(p)) into QPA constructor calls, and assemble
the enveloping-algebra Hochschild route (HH^n = Ext^n_{A^e}(A,A)) since QPA ships
no HH function (spec §2). Scripts are GAP source strings; no floats.

Scope: cross-check runs on algebras presented over QQ or a prime field GF(p) --
the fields QPA supports exactly. Number-field CC entries and GF(p^n) are out of
the cross-check scope (raise QpaUnavailableError with that reason)."""
from __future__ import annotations


def _gap_field(domain) -> str:
    """QPA base field literal for a quiverlab Domain (QQ or GF(p))."""
    name = getattr(domain, "name", "")
    char = domain.characteristic()
    if char == 0 and name in ("QQ", "Rationals"):
        return "Rationals"
    if char > 0 and getattr(domain, "degree", 1) == 1:      # prime field GF(p)
        return f"GF({char})"
    raise ValueError(
        f"QPA cross-check supports QQ or prime GF(p) only; got domain {name!r} "
        f"(characteristic {char}). Number-field CC and GF(p^n) are out of scope."
    )


def quiver_and_algebra_script(algebra) -> str:
    """Emit GAP source binding `A := kQ/rels` (or `A := kQ` when no relations),
    reconstructing the quiver from algebra._quiver and its relations. Vertices are
    numbered 1..n in quiver order; arrows carry their quiverlab names. VERIFY the
    relation-string translation against QPA's element grammar at execution."""
    Q = algebra._quiver
    verts = list(Q.vertices)
    idx = {v: i + 1 for i, v in enumerate(verts)}               # QPA is 1-based
    arrows = [[idx[Q.source(a)], idx[Q.target(a)], a] for a in Q.arrows]
    arrow_gap = ", ".join(f'[{s}, {t}, "{name}"]' for s, t, name in arrows)
    field = _gap_field(algebra.domain)
    lines = [
        f"Q := Quiver({len(verts)}, [{arrow_gap}]);;",
        f"kQ := PathAlgebra({field}, Q);;",
    ]
    rels = getattr(algebra, "_relations", None)
    if rels:
        # Each relation is a linear combo of parallel paths; render as a QPA element
        # over the generators kQ.<arrow>. VERIFY the exact element syntax at run.
        terms = _relations_to_gap(rels, "kQ")
        lines.append(f"rels := [{terms}];;")
        lines.append("A := kQ/rels;;")
    else:
        lines.append("A := kQ;;")
    return "\n".join(lines)


def _relations_to_gap(relations, kq: str) -> str:
    """Render quiverlab relations (tuples of (coeff, word)) as QPA algebra
    elements `sum coeff * kQ.a1*kQ.a2*...`. Coefficients are exact integers/
    fractions -> GAP integer/rational literals (no floats). VERIFY grammar."""
    out = []
    for rel in relations:
        parts = []
        for coeff, word in rel.terms:              # Relation.terms: ((Fraction, (arrow,...)),...)
            num, den = coeff.numerator, coeff.denominator
            path = "*".join(f"{kq}.{a}" for a in word)
            scal = f"{num}" if den == 1 else f"({num}/{den})"
            parts.append(f"{scal}*{path}")
        out.append(" + ".join(parts))
    return ", ".join(out)


def hochschild_dims_script(algebra, top: int) -> str:
    """Append the enveloping-algebra HH route to the algebra script, binding a GAP
    list `hh := [dim HH^0, ..., dim HH^top]`. HH^n(A) = Ext^n_{A^e}(A,A).

    Dim read: QPA's `ExtAlgebraGenerators(M, n)` returns a list whose FIRST component
    is the list of `dim Ext^i(M, M)` for `i = 0..n` (the standard QPA idiom for reading
    an Ext/HH dimension series). Here `M = AA` is `A` as an `A^e`-module, so
    `ExtAlgebraGenerators(AA, top)[1]` is exactly `[dim HH^0, ..., dim HH^top]`.

    VERIFY AT EXECUTION: (a) `AlgebraAsModuleOverEnvelopingAlgebra` -- the QPA op that
    presents A as a right A^e-module (QPA manual ch.6/8); if the name differs, build the
    bimodule explicitly from `EnvelopingAlgebra` + the regular representation. (b) that
    `ExtAlgebraGenerators(M, n)[1]` is the degreewise-dimension component (GAP is
    1-indexed: `[1]` is the first return value); the `[1,0,0]` fixture is the oracle."""
    base = quiver_and_algebra_script(algebra)
    return base + "\n" + "\n".join([
        "Ae := EnvelopingAlgebra(A);;",
        "AA := AlgebraAsModuleOverEnvelopingAlgebra(A);;    # VERIFY constructor name",
        f"info := ExtAlgebraGenerators(AA, {top});;",
        f"hh := info[1];;                                    # dims Ext^0..Ext^{top}  [VERIFY component]",
        "Print(hh);",
    ])


def module_self_ext_dims_script(algebra, dimvec_M, top: int) -> str:
    """Bind `ext := [dim Ext^0(M,M), ..., dim Ext^top(M,M)]` (self-Ext of one module
    given by its dimension vector) via the SAME idiom `ExtAlgebraGenerators(M, top)[1]`.

    Self-Ext keeps the cross-check on the one confirmed QPA idiom. Distinct-module
    Ext(M,N) (M != N) needs `ExtOverAlgebra` + iterated `NthSyzygy` instead and is a
    flagged post-v1 extension. VERIFY the `RightModuleOverPathAlgebra` args + the
    `[1]` component read at execution."""
    base = quiver_and_algebra_script(algebra)
    return base + "\n" + "\n".join([
        f"M := RightModuleOverPathAlgebra(A, {list(dimvec_M)}, []);;   # VERIFY args",
        f"info := ExtAlgebraGenerators(M, {top});;",
        f"ext := info[1];;                                   # dims Ext^0..Ext^{top}(M,M)  [VERIFY]",
        "Print(ext);",
    ])
```

`src/quiverlab/qpa/crosscheck.py`:

```python
"""A.crosscheck(...): independent QPA recomputation of Hochschild dims / module Ext
for validation workflows (spec §5 c.12, §8 ring 3). Returns a CrosscheckReport;
never silently disagrees -- .assert_agree() raises on mismatch."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QpaUnavailableError
from quiverlab.qpa import scripts, session


@dataclass
class CrosscheckReport:
    what: str                 # "hochschild" | "module_ext"
    ours: list                # quiverlab dims
    qpa: list                 # QPA dims
    agree: bool

    def assert_agree(self):
        if not self.agree:
            raise AssertionError(
                f"QPA cross-check DISAGREES on {self.what}: quiverlab {self.ours} "
                f"vs QPA {self.qpa}")
        return self


def _read_int_list(gap_value) -> list:
    """Convert a GAP list of integers into a Python list[int] (exact; no floats)."""
    return [int(x) for x in gap_value]


def crosscheck_hochschild(algebra, top: int) -> CrosscheckReport:
    session.require_gap()
    ours = algebra.hochschild_cohomology(top).dims
    gap = session.run(scripts.hochschild_dims_script(algebra, top) + " hh;")
    qpa = _read_int_list(gap)
    return CrosscheckReport("hochschild", list(ours), qpa, list(ours) == qpa)


def crosscheck_module_ext(algebra, M, top: int) -> CrosscheckReport:
    """Self-Ext Ext^*(M, M) vs QPA (via ExtAlgebraGenerators). Distinct-module
    Ext(M, N) is a flagged post-v1 extension (needs ExtOverAlgebra + syzygies)."""
    session.require_gap()
    ours = [algebra.ext(M, M, n).dimension() for n in range(top + 1)]  # Plan 05 surface
    gap = session.run(
        scripts.module_self_ext_dims_script(algebra, M.dimension_vector(), top) + " ext;")
    qpa = _read_int_list(gap)
    return CrosscheckReport("module_ext", list(ours), qpa, list(ours) == qpa)


def crosscheck(algebra, what: str, *args, **kwargs) -> CrosscheckReport:
    """Dispatch. what="hochschild" -> crosscheck_hochschild(algebra, top);
    what="module_ext" -> crosscheck_module_ext(algebra, M, top) (self-Ext)."""
    if what == "hochschild":
        return crosscheck_hochschild(algebra, *args, **kwargs)
    if what == "module_ext":
        return crosscheck_module_ext(algebra, *args, **kwargs)
    raise QpaUnavailableError(f"unknown cross-check {what!r}",
                              hint='use "hochschild" or "module_ext"')
```

`src/quiverlab/qpa/__init__.py`:

```python
"""quiverlab.qpa: optional GAP/QPA cross-check backend (spec §5 component 12).
Imported lazily by Algebra.crosscheck; importing this module does NOT import GAP
(that happens on first gap_available()/crosscheck call)."""
from quiverlab.qpa.session import gap_available, require_gap  # noqa: F401
from quiverlab.qpa.crosscheck import crosscheck, CrosscheckReport  # noqa: F401

__all__ = ["gap_available", "require_gap", "crosscheck", "CrosscheckReport"]
```

`src/quiverlab/core/algebra.py` — add the method (lazy import keeps GAP out of the core import path):

```python
    def crosscheck(self, what="hochschild", *args, **kwargs):
        """Independently recompute an invariant via the optional QPA backend and
        compare (spec §5 c.12). Requires `pip install quiverlab[qpa]`; raises
        QpaUnavailableError otherwise. Examples:
            A.crosscheck("hochschild", 3)          # HH^0..HH^3 vs QPA enveloping route
            A.crosscheck("module_ext", M, 4)       # Ext^0..Ext^4(M,M) vs QPA (self-Ext)
        Returns a CrosscheckReport; call .assert_agree() to fail loudly on mismatch."""
        from quiverlab.qpa.crosscheck import crosscheck as _cc
        return _cc(self, what, *args, **kwargs)
```

Export `QpaUnavailableError` from `src/quiverlab/__init__.py` (add to the `from quiverlab.errors import (...)` list and to `__all__`).

`tests/qpa/test_session_guard.py` (marked `qpa` by path — but the guard tests themselves must run everywhere; put the guard-behaviour tests here and mark them `fast` explicitly, keeping the actual GAP tests in Task 5):

```python
"""The [qpa] guard behaves without GAP installed (runs on every cell).

These are marked `fast` explicitly so they run in the normal matrix even though
they live under tests/qpa/ (the path default would be `qpa`)."""
import pytest

from quiverlab import Quiver, GF, QpaUnavailableError
from quiverlab.qpa import gap_available

pytestmark = pytest.mark.fast


def test_gap_available_is_boolean_and_cheap():
    assert isinstance(gap_available(), bool)      # no exception when GAP absent


def test_crosscheck_raises_cleanly_without_backend():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(2))
    if gap_available():
        pytest.skip("GAP present; guard-absence path not exercised here")
    with pytest.raises(QpaUnavailableError) as e:
        A.crosscheck("hochschild", 2)
    assert "quiverlab[qpa]" in str(e.value) or "WSL" in str(e.value)
```

- [ ] **Step 1: Write the failing test** (`tests/qpa/test_session_guard.py`, above) + `tests/qpa/__init__.py`.
- [ ] **Step 2: Run to verify it fails**
  Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/qpa/test_session_guard.py -q`
  Expected: FAIL — `quiverlab.qpa` / `QpaUnavailableError` / `Algebra.crosscheck` do not exist.
- [ ] **Step 3: Write the package** (`errors.py` addition, `qpa/` modules, `Algebra.crosscheck`, `__init__` export) above.
- [ ] **Step 4: Run focused suite**
  Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/qpa/test_session_guard.py tests/test_no_floats.py -q`
  Expected: pass (guard returns False cleanly; `crosscheck` raises `QpaUnavailableError`; float gate green — the `qpa` package holds no float literals).
- [ ] **Step 5: Commit**
  ```bash
  git add src/quiverlab/qpa src/quiverlab/errors.py src/quiverlab/__init__.py src/quiverlab/core/algebra.py tests/qpa/
  git commit -m "feat(qpa): optional GAP/QPA backend — session guard, script builders, A.crosscheck

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
  ```

---

### Task 5: QPA cross-check test suite + the `kA_2`/`GF(2)` enveloping-algebra HH fixture

**Files:**
- Create: `tests/qpa/test_crosscheck.py`, `tests/qpa/test_gap_fixture.py`

**Interfaces:**
- Consumes: `Algebra.crosscheck`, `quiverlab.qpa.session` (`gap_available`, `should_skip_qpa`), the GAP scripts.
- Produces: the `qpa`-marked cross-check battery (spec §8 ring 3). Locally these **skip** (no GAP); in the `qpa.yml` CI job (`QUIVERLAB_REQUIRE_QPA=1`) they are **mandatory**. One fixture writes the full GAP script and its expected output, flagged **VERIFY AT EXECUTION** (cannot run on this machine).

**The worked fixture — HH of `kA_2` over `GF(2)`.** `kA_2` = path algebra of `1 --a--> 2`, no relations (hereditary, directed; `= T_2`, upper-triangular 2×2). Its Hochschild cohomology is `HH^0 = k` (center = scalars), `HH^{≥1} = 0` (Happel: hereditary, tree quiver, no oriented cycles). Over `GF(2)`, identical (dimensions characteristic-independent here). quiverlab returns `[1, 0, 0]`; the QPA enveloping-algebra route must return `[1, 0, 0]`.

The exact GAP script the builder emits (Task 4) for this algebra, annotated:

```gap
LoadPackage("qpa");;
Q := Quiver(2, [[1, 2, "a"]]);;      # vertices 1,2 ; arrow a: 1 -> 2
kQ := PathAlgebra(GF(2), Q);;
A := kQ;;                            # no relations => A = kQ (hereditary)
Ae := EnvelopingAlgebra(A);;         # A^e = A^op (x) A
AA := AlgebraAsModuleOverEnvelopingAlgebra(A);;   # A as right A^e-module  [VERIFY name]
res := ProjectiveResolutionOfPathAlgebraModule(AA, 3);;
hh := List([0..2], n -> <dim of H^n Hom_{A^e}(res, AA)>);;   # [VERIFY the dim read]
Print(hh);                          # EXPECTED OUTPUT:  [ 1, 0, 0 ]   [VERIFY AT EXECUTION]
```

**Expected GAP stdout:** `[ 1, 0, 0 ]` — flagged **VERIFY AT EXECUTION** (GAP is not installed locally; confirm both the `AlgebraAsModuleOverEnvelopingAlgebra` constructor name and the `Ext^n` dimension read against the QPA manual when the `qpa.yml` job first runs; the *mathematics* `[1,0,0]` is certain — Happel — but the QPA call surface must be pinned).

`tests/qpa/test_gap_fixture.py`:

```python
"""The kA_2/GF(2) enveloping-algebra HH fixture: quiverlab vs QPA both give
[1,0,0] (spec §8 ring 3). Skips locally (no GAP); mandatory in the qpa.yml CI job."""
import pytest

from quiverlab import Quiver, GF
from quiverlab.qpa import session, scripts


def _ka2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(2))


@pytest.mark.skipif(session.should_skip_qpa(), reason="[qpa] backend not installed")
def test_ka2_gf2_hochschild_matches_qpa():
    A = _ka2()
    assert A.hochschild_cohomology(2).dims == [1, 0, 0]     # Happel; char-independent
    report = A.crosscheck("hochschild", 2)
    report.assert_agree()
    assert report.qpa == [1, 0, 0]                          # VERIFY AT EXECUTION


def test_ka2_gf2_script_is_wellformed_even_without_gap():
    """The emitted GAP source is stable and mentions the enveloping route -- this
    part runs everywhere (no GAP needed), pinning the script the CI job will run."""
    src = scripts.hochschild_dims_script(_ka2(), 2)
    assert 'Quiver(2, [[1, 2, "a"]])' in src
    assert "PathAlgebra(GF(2), Q)" in src
    assert "EnvelopingAlgebra(A)" in src
```

`tests/qpa/test_crosscheck.py`:

```python
"""QPA cross-check battery over a small zoo (spec §8 ring 3). qpa-marked: skips
locally, mandatory under QUIVERLAB_REQUIRE_QPA=1 in CI."""
import pytest

from quiverlab import Quiver, GF, QQ
from quiverlab.qpa import session

# skip locally (no GAP); under QUIVERLAB_REQUIRE_QPA=1 the predicate is False, so the
# tests RUN and fail naturally if GAP is missing/broken (no silent green skip).
pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


@pytest.mark.parametrize("field", [GF(2), GF(3), QQ])
def test_hochschild_crosscheck_small_zoo(field):
    # commutative square kQ/(a*b - c*d) = kA_2 (x) kA_2 : HH = [1,0,0] (Kunneth)
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    A = Q.algebra(relations=["a*b - c*d"], field=field)
    A.crosscheck("hochschild", 2).assert_agree()


def test_module_self_ext_crosscheck():
    # self-Ext of the simple S1 in kA_2 over GF(2): Ext^0=Hom(S1,S1)=1, Ext^{>=1}=0
    # (no loop at vertex 1). ExtAlgebraGenerators(M, n)[1] gives the dim series.
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(2))
    A.crosscheck("module_ext", A.simple(1), 2).assert_agree()   # -> [1, 0, 0]
```

- [ ] **Step 1: Write the tests** (above).
- [ ] **Step 2: Run locally (expect SKIP)**
  Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/qpa/ -q`
  Expected: the GAP-dependent tests **skip** ("[qpa] backend not installed"); the script-shape test and the guard tests pass. No failures, no errors. (`should_skip_qpa()` returns True only because GAP is absent AND `QUIVERLAB_REQUIRE_QPA` is unset; under `=1` it returns False, so the tests run and would fail on absent GAP.)
- [ ] **Step 3: (no src changes — Task 4 supplies the machinery)**
  If a script-shape assertion fails, fix the Task-4 builder (do not weaken the assertion).
- [ ] **Step 4: Confirm the marker + float gate**
  Run:
  ```bash
  cd /Users/marco/Desktop/HomologicalNetworks/quiverlab
  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q --collect-only -m qpa | tail -3   # tests/qpa/* collected under qpa
  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/qpa tests/test_no_floats.py -q
  ```
  Expected: `-m qpa` collects the GAP tests (the guard/script tests are explicitly `fast`); the local run skips cleanly; float gate green.
- [ ] **Step 5: Commit**
  ```bash
  git add tests/qpa/test_crosscheck.py tests/qpa/test_gap_fixture.py
  git commit -m "test(qpa): enveloping-algebra HH cross-check + kA_2/GF(2) GAP fixture (skips w/o GAP)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
  ```

---

### Task 6: GitHub Actions CI — test matrix + engine-path legs + float-gate/lint

**Files:**
- Create: `.github/workflows/ci.yml`, `tests/release/test_workflows.py`

**Interfaces:**
- Consumes: the marker vocabulary (Task 2), the `[dev]`/`[fast]` extras (Task 3).
- Produces: the concrete CI split. `.github/` is greenfield.

**The CI-split design (concrete).** The full suite is ~19 min heavy (deep engine/CS/resolution tests). Split:
- **`fast` matrix — 12 cells** (`{ubuntu, macos, windows} × {3.10, 3.11, 3.12, 3.13}`): install `.[dev]` (pure Python — no numba, whose wheels lag on the newest Python), run `pytest -m fast`. Broad cross-platform coverage, each cell a few minutes. `qpa`-marked tests are excluded (the guard tests are `fast` and pass with GAP absent, incl. Windows).
- **`deep` — 2 Linux cells** (`py3.12`): the ~19-min full suite, once on the **numba** engine path (`.[dev,fast]`) and once on the **pure** path (`.[dev]` + `QUIVERLAB_NO_NUMBA=1`) — this is the spec's "both paths equality-tested" gate. `qpa` tests skip here (no GAP).
- **`lint` — 1 Linux cell**: the float-ban AST gate + the release-metadata/marker/workflow tests.
- The **QPA** cross-check is a **separate workflow** (Task 7, `qpa.yml`) so its heavy GAP install never blocks the matrix.

`.github/workflows/ci.yml`:

```yaml
name: CI

on:
  push:
    branches: [main]
  pull_request:

concurrency:
  group: ci-${{ github.ref }}
  cancel-in-progress: true

env:
  MPLBACKEND: "Agg"          # headless matplotlib on every job (matplotlib is a hard dep)

jobs:
  fast:
    name: fast · ${{ matrix.os }} · py${{ matrix.python }}
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]
        python: ["3.10", "3.11", "3.12", "3.13"]
    env:
      NUMBA_NUM_THREADS: "2"
      OMP_NUM_THREADS: "2"
    steps:
      - uses: actions/checkout@v7
      - uses: actions/setup-python@v6
        with:
          python-version: ${{ matrix.python }}
      - run: python -m pip install --upgrade pip
      - run: pip install -e ".[dev]"          # pure path; numba is a separate deep leg
      - run: python -m pytest -m fast -q

  deep:
    name: deep · linux · py3.12 · ${{ matrix.leg }}
    runs-on: ubuntu-latest
    strategy:
      fail-fast: false
      matrix:
        include:
          - { leg: numba, extras: ".[dev,fast]", no_numba: "0" }
          - { leg: pure,  extras: ".[dev]",      no_numba: "1" }
    env:
      NUMBA_NUM_THREADS: "2"
      OMP_NUM_THREADS: "2"
      QUIVERLAB_NO_NUMBA: ${{ matrix.no_numba }}
    steps:
      - uses: actions/checkout@v7
      - uses: actions/setup-python@v6
        with:
          python-version: "3.12"
      - run: python -m pip install --upgrade pip
      - run: pip install -e "${{ matrix.extras }}"
      - run: python -m pytest -q              # full suite (fast + deep); qpa skips (no GAP)

  lint:
    name: float-gate + metadata
    runs-on: ubuntu-latest
    env:
      NUMBA_NUM_THREADS: "2"
      OMP_NUM_THREADS: "2"
    steps:
      - uses: actions/checkout@v7
      - uses: actions/setup-python@v6
        with:
          python-version: "3.12"
      - run: python -m pip install --upgrade pip
      - run: pip install -e ".[dev]"
      - run: python -m pytest tests/test_no_floats.py tests/release -q
```

`tests/release/test_workflows.py` (dependency-free structural check — no PyYAML needed; GitHub validates the YAML itself on push, and the `on:`→`True` YAML-1.1 key gotcha makes `safe_load` structural asserts brittle):

```python
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
```

This file grows one assertion per workflow: `test_ci_matrix_covers_os_and_python`
now (Task 6); `test_qpa_workflow_is_linux_and_mandatory` in Task 7;
`test_docs_workflow_deploys_pages` in Task 8; `test_paper_workflow_builds_pdf` in
Task 9; `test_release_workflow_trusted_publishing` in Task 10 — each added together
with the workflow it checks, so the suite is green at every commit (no assertion
ever references a not-yet-written file).

- [ ] **Step 1: Write `tests/release/test_workflows.py`** with **only** `test_ci_matrix_covers_os_and_python` (above) — it references `ci.yml`, which Step 3 creates.
- [ ] **Step 2: Run to verify it fails**
  Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/release/test_workflows.py -q`
  Expected: FAIL — `.github/workflows/ci.yml` does not exist yet.
- [ ] **Step 3: Write `.github/workflows/ci.yml`** (above).
- [ ] **Step 4: Validate + run**
  Run:
  ```bash
  cd /Users/marco/Desktop/HomologicalNetworks/quiverlab
  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/release/test_workflows.py::test_ci_matrix_covers_os_and_python tests/test_no_floats.py -q
  ```
  Expected: pass. (Optional: run `actionlint` if available to lint the YAML; GitHub validates it on push regardless.)
- [ ] **Step 5: Commit**
  ```bash
  git add .github/workflows/ci.yml tests/release/test_workflows.py
  git commit -m "ci: test matrix (OS x py3.10-3.13), deep numba+pure legs, float-gate/metadata lint

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
  ```

---

### Task 7: QPA cross-check CI job — Linux, GAP+QPA installed, mandatory `-m qpa`

**Files:**
- Create: `.github/workflows/qpa.yml`
- Modify: `tests/release/test_workflows.py` (add `test_qpa_workflow_is_linux_and_mandatory`)

**Interfaces:**
- Consumes: the `[qpa]` extra (Task 3), the `qpa` marker + `QUIVERLAB_REQUIRE_QPA` guard (Tasks 4–5).
- Produces: the CI-only Linux cross-check job. It is **separate** from `ci.yml` (heavy GAP install; runs on push-to-main, manual dispatch, and weekly — not on every PR).

`.github/workflows/qpa.yml`:

```yaml
name: QPA cross-check

on:
  push:
    branches: [main]
  workflow_dispatch:
  schedule:
    - cron: "0 6 * * 1"          # weekly (Mon 06:00 UTC); GAP install is heavy

jobs:
  qpa:
    name: QPA cross-check · linux · py3.12
    runs-on: ubuntu-latest
    env:
      NUMBA_NUM_THREADS: "2"
      OMP_NUM_THREADS: "2"
      MPLBACKEND: "Agg"
      QUIVERLAB_REQUIRE_QPA: "1"   # absent/broken QPA is a HARD failure in this job
    steps:
      - uses: actions/checkout@v7
      - uses: actions/setup-python@v6
        with:
          python-version: "3.12"   # passagemath-gap wheels: py3.11-3.14, linux
          cache: pip               # cache the heavy passagemath-gap wheel between runs
          cache-dependency-path: pyproject.toml
      - name: Install GAP + QPA (passagemath-gap) and quiverlab[dev,qpa]
        run: |
          python -m pip install --upgrade pip
          # passagemath-gap[qpa] pulls prebuilt GAP + QPA v1.37 + GBNP (manylinux wheel).
          # Fallback if the [qpa] sub-extra name drifts: pip install "passagemath-gap[full]"
          pip install -e ".[dev,qpa]" \
            || { pip install "passagemath-gap[full]" && pip install -e ".[dev]"; }
      - name: Sanity — libgap imports and QPA loads
        run: python -c "from quiverlab.qpa import gap_available; assert gap_available(), 'QPA failed to load'"
      - name: Run the QPA cross-check suite (mandatory)
        run: python -m pytest -m qpa -q
```

Notes:
- **VERIFY AT EXECUTION** (Task 4/5 research): the `passagemath-gap[qpa]` extra name (fallback `[full]`), and the two QPA call surfaces in `scripts.py` (`AlgebraAsModuleOverEnvelopingAlgebra`, the `Ext^n` dim read). The first time this job runs, confirm both against the QPA manual and pin them; the `[1,0,0]` fixture output is the oracle.
- `QUIVERLAB_REQUIRE_QPA=1` is folded **into the skip predicate itself** (`session.should_skip_qpa()` returns False under it), so the qpa tests **run** and fail naturally if GAP is missing/broken — no silent green-skip. (A `setup_module()`-based escalation would be dead code: `skipif` is evaluated at collection, before any `setup_module()` runs.) The "Sanity — QPA loads" step is a second, explicit hard gate.

- [ ] **Step 1: Add `test_qpa_workflow_is_linux_and_mandatory`** to `tests/release/test_workflows.py` (as in Task 6's file). Expected FAIL (no `qpa.yml`).
- [ ] **Step 2: Run to verify it fails** — `... -m pytest tests/release/test_workflows.py -q`.
- [ ] **Step 3: Write `.github/workflows/qpa.yml`** (above).
- [ ] **Step 4: Validate + run**
  Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/release/test_workflows.py -q`
  Expected: the CI + QPA workflow assertions pass. (The GAP job itself only exercises on GitHub's Linux runners; locally we only assert its shape.)
- [ ] **Step 5: Commit**
  ```bash
  git add .github/workflows/qpa.yml tests/release/test_workflows.py
  git commit -m "ci(qpa): Linux QPA cross-check job (passagemath-gap+QPA, mandatory -m qpa)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
  ```

---

### Task 8: Docs site — mkdocs-material, auto API reference, executed tutorials, internals, Pages deploy

**Files:**
- Create: `mkdocs.yml`, `docs/index.md`, `scripts/gen_ref_pages.py`, `scripts/gen_bibliography.py`, `docs/javascripts/mathjax.js`, `.github/workflows/docs.yml`
- Create: `tests/release/test_docs.py`
- Modify: `tests/release/test_workflows.py` (add `test_docs_workflow_deploys_pages`)

**Interfaces:**
- Consumes: the `[docs]` extra (Task 3), `docs/tutorials/*.ipynb`, `docs/internals/*.md`, the packaged `quiverlab.citations.references_bib_path()` bib, the installed `quiverlab` package (mkdocstrings introspects it).
- Produces: a `mkdocs-material` site that builds `--strict`, with auto API reference (mkdocstrings + gen-files/literate-nav/section-index recipe), CI-executed tutorial notebooks, the "Under the hood" internals chapters, and GitHub Pages deployment via the native artifact flow. `site_url = https://marcoarmenta.github.io/quiverlab/`.

**Tooling (research, 2026):** `mkdocs-material` 9.7.x (note: Material entered maintenance mode early 2026 — fine for a 0.x release), `mkdocstrings` 1.0.x + `mkdocstrings-python` 2.0.x (the `[python]` extra pulls the matching handler), `mkdocs-jupyter` 0.26.x, `mkdocs-gen-files` 0.6.x, `mkdocs-literate-nav` 0.6.x, `mkdocs-section-index` 0.3.x, `mike` 2.2.x.

**Versioned-docs decision (`mike`): DEFER for the 0.x line.** Justification: a single "latest" site is simpler while the API churns pre-1.0; a version selector adds noise and every old build would have to be kept green. `mike` (already in the `[docs]` extra) is adopted at **1.0 / JOSS acceptance**, when old-version docs must be pinned. Enabling it later is a two-line change (add `extra.version.provider: mike` to `mkdocs.yml` and switch `docs.yml` to `mike deploy --push --update-aliases <ver> latest`); recorded here, not built now.

`scripts/gen_ref_pages.py` (the official mkdocstrings recipe; `src/` layout — `Path(__file__).parent.parent` == repo root, package at `src/quiverlab/`):

```python
"""Generate one API-reference page per module + a literate-nav SUMMARY (mkdocstrings
recipe). Run by the gen-files plugin at build time; writes virtual files only."""
from pathlib import Path

import mkdocs_gen_files

nav = mkdocs_gen_files.Nav()
root = Path(__file__).parent.parent
src = root / "src"

for path in sorted(src.rglob("*.py")):
    module_path = path.relative_to(src).with_suffix("")
    doc_path = path.relative_to(src).with_suffix(".md")
    full_doc_path = Path("reference", doc_path)
    parts = tuple(module_path.parts)
    if parts[-1] == "__init__":
        parts = parts[:-1]
        doc_path = doc_path.with_name("index.md")
        full_doc_path = full_doc_path.with_name("index.md")
    elif parts[-1] == "__main__":
        continue
    if not parts or parts[-1].startswith("_"):
        continue                                  # skip private modules (_kernels, etc.)
    nav[parts] = doc_path.as_posix()
    with mkdocs_gen_files.open(full_doc_path, "w") as fd:
        fd.write(f"::: {'.'.join(parts)}")
    mkdocs_gen_files.set_edit_path(full_doc_path, path.relative_to(root))

with mkdocs_gen_files.open("reference/SUMMARY.md", "w") as nav_file:
    nav_file.writelines(nav.build_literate_nav())
```

`mkdocs.yml` (complete):

```yaml
site_name: quiverlab
site_description: Quivers with relations and Hochschild theory, exactly, for algebraists
site_url: https://marcoarmenta.github.io/quiverlab/
repo_url: https://github.com/MarcoArmenta/quiverlab
repo_name: MarcoArmenta/quiverlab
copyright: "MIT © 2026 Marco Armenta"

theme:
  name: material
  features:
    - navigation.sections
    - navigation.indexes         # required by section-index
    - navigation.top
    - content.code.copy
    - search.suggest
  palette:
    - media: "(prefers-color-scheme: light)"
      scheme: default
      toggle: { icon: material/brightness-7, name: Switch to dark mode }
    - media: "(prefers-color-scheme: dark)"
      scheme: slate
      toggle: { icon: material/brightness-4, name: Switch to light mode }

plugins:
  - search
  - gen-files:
      scripts:
        - scripts/gen_ref_pages.py       # API reference pages from docstrings
        - scripts/gen_bibliography.py     # References page from the PACKAGED bib
  - literate-nav:
      nav_file: SUMMARY.md
  - section-index
  - mkdocstrings:
      handlers:
        python:
          options:
            docstring_style: google
            show_source: true
            show_root_heading: true
            merge_init_into_class: true
            separate_signature: true
            show_signature_annotations: true
  - mkdocs-jupyter:
      execute: true               # CI-executes the tutorials (spec §10)
      include_source: true        # "download notebook" link
      allow_errors: false
      theme: auto

nav:
  - Home: index.md
  - Tutorials:
      - tutorials/01-exact-fields.ipynb
      - tutorials/02-quivers-and-algebras.ipynb
      - tutorials/03-hochschild.ipynb
  - Under the hood:
      # LIST EVERY docs/internals/NN-*.md that exists at write time (Plans 04-07 add
      # chapters 08-11 = CS, modules, families, viz/trace BEFORE Plan 08 runs). The
      # nav-coverage test (test_docs.py::test_nav_covers_internals_and_tutorials) FAILS
      # loudly if any internals chapter or tutorial notebook is missing from this nav,
      # and `validation.omitted_files: warn` -> error under --strict is the CI backstop.
      - internals/README.md
      - internals/01-exact-fields.md
      - internals/02-quivers-relations.md
      - internals/03-algebra.md
      - internals/04-hochschild-bar.md
      - internals/05-resolutions.md
      - internals/06-invariants.md
      - internals/07-dispatch.md
      # - internals/08-chouhy-solotar.md   (Plan 04)  <- add when present
      # - internals/09-modules.md          (Plan 05)
      # - internals/10-...                  (Plan 05/06/07)
  - Development:
      - development/release.md     # Task 12 (CI/PyPI/JOSS overview; the QPA route math)
  - References: bibliography.md    # generated by scripts/gen_bibliography.py from the packaged bib
  - API Reference: reference/      # literate-nav fills this from the generated SUMMARY.md

markdown_extensions:
  - pymdownx.arithmatex: { generic: true }
  - pymdownx.highlight
  - pymdownx.superfences
  - admonition
  - pymdownx.details
  - toc: { permalink: true }

extra_javascript:
  - javascripts/mathjax.js
  # PINNED third-party CDN (never a floating @3 tag) from a reputable source, same
  # supply-chain hygiene as stripping polyfill.io from the trace HTML. Fully-offline
  # option (Marco's call): vendor this asset under docs/javascripts/ and reference the
  # local path instead, so the published docs load no third-party JS at all.
  - https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-mml-chtml.js

# Keep the plan/spec markdown out of the site (they live under docs/). The
# bibliography is NOT in docs/ (it is packaged in src/quiverlab/citations/); the
# References page is generated from it by gen_bibliography.py, so nothing to exclude.
exclude_docs: |
  plans/
  specs/
  tutorials/README.md

validation:
  omitted_files: warn             # -> error under --strict: any internals chapter not in nav fails the build
  absolute_links: ignore          # internals chapters cross-link to plans (excluded)
  unrecognized_links: ignore
  nav:
    not_found: warn
```

`scripts/gen_bibliography.py` (the docs "References" page — generated at build time
from the **packaged** bib; no copy lands in the repo tree):

```python
"""Generate the docs References page from the single packaged bibliography
(src/quiverlab/citations/references.bib) via the library's own accessor. Run by
the gen-files plugin; writes a virtual bibliography.md only."""
import mkdocs_gen_files

from quiverlab import bibliography              # Plan 06: grouped/annotated text
from quiverlab.citations import references_bib_path

_raw = references_bib_path().read_text(encoding="utf-8")
with mkdocs_gen_files.open("bibliography.md", "w") as fd:
    fd.write("# References\n\n")
    fd.write("The curated, verified bibliography that quiverlab cites, rendered from "
             "the single packaged source `src/quiverlab/citations/references.bib`.\n\n")
    fd.write(bibliography())                    # annotated, grouped by algorithm/family/field
    fd.write("\n\n## Raw BibTeX\n\n```bibtex\n" + _raw + "\n```\n")
```

`docs/javascripts/mathjax.js` (Material's documented MathJax hook — plays well with mkdocs-jupyter's rendered math):

```javascript
window.MathJax = {
  tex: { inlineMath: [["\\(", "\\)"]], displayMath: [["\\[", "\\]"]],
         processEscapes: true, processEnvironments: true },
  options: { ignoreHtmlClass: ".*|", processHtmlClass: "arithmatex" }
};
document$.subscribe(() => { MathJax.startup.output.clearCache(); MathJax.typesetClear();
                            MathJax.texReset(); MathJax.typesetPromise(); });
```

`docs/index.md` (landing page — mirrors the README quickstart, links out):

```markdown
# quiverlab

**Quivers with relations and Hochschild theory, exactly, for algebraists.**

Finite-dimensional algebras `kQ/I` over the complex numbers (exactly — no floating
point, ever) and all finite fields: certified finite-dimensionality, Hochschild
(co)homology with cup products and Gerstenhaber brackets, the first full
Chouhy–Solotar resolution, module Ext, Cartan/Coxeter invariants, drawings, and
worked-steps traces.

```python
from quiverlab import Quiver, CC

Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})
A = Q.algebra(relations=["a*b"], field=CC)
print(A.hochschild_cohomology(3))
```

- **Tutorials** — start here (executable notebooks).
- **Under the hood** — how each object is represented and each number produced.
- **API Reference** — every public function and class.
- **Web GUI** — compute without installing anything (Plan 09).
- **Cite** — see the JOSS paper and `CITATION.cff`.
```

`.github/workflows/docs.yml` (native Pages artifact flow; build verifies `--strict` on every PR, deploy only on `main`/tags):

```yaml
name: Docs

on:
  push:
    branches: [main]
    tags: ["v*"]
  pull_request:

permissions:
  contents: read
  pages: write
  id-token: write

concurrency:
  group: pages
  cancel-in-progress: false

jobs:
  build:
    runs-on: ubuntu-latest
    env:
      NUMBA_NUM_THREADS: "2"
      OMP_NUM_THREADS: "2"
      MPLBACKEND: "Agg"                       # headless matplotlib for tutorial execution
    steps:
      - uses: actions/checkout@v7
      - uses: actions/setup-python@v6
        with:
          python-version: "3.12"
      - run: python -m pip install --upgrade pip
      - run: pip install -e ".[dev,docs]"     # installs quiverlab (mkdocstrings + notebook exec)
      - run: mkdocs build --strict            # executes tutorials, builds API ref, fails on warnings
      - uses: actions/configure-pages@v5
      - uses: actions/upload-pages-artifact@v5
        with:
          path: site

  deploy:
    needs: build
    if: github.ref == 'refs/heads/main' || startsWith(github.ref, 'refs/tags/')
    runs-on: ubuntu-latest
    environment:
      name: github-pages
      url: ${{ steps.deployment.outputs.page_url }}
    steps:
      - id: deployment
        uses: actions/deploy-pages@v5
```

(One-time repo setting Marco does: **Settings → Pages → Source = GitHub Actions.**)

**Action pins (bumped past the Node-20 runner cutoff, hard-breaks 2026-09-16):**
`actions/checkout@v7`, `actions/setup-python@v6`, `actions/upload-artifact@v7`,
`actions/download-artifact@v8`, `actions/configure-pages@v5`,
`actions/upload-pages-artifact@v5`, `actions/deploy-pages@v5`,
`pypa/gh-action-pypi-publish@release/v1` (moving tag by design),
`openjournals/openjournals-draft-action@master` (**moving ref — VERIFY/pin a commit
SHA at execution if reproducibility matters**). Re-confirm the latest majors at
execution; these are chosen to run on the post-Node-20 runners.

`tests/release/test_docs.py`:

```python
"""The docs site config is present and canonical; a real `mkdocs build` succeeds
when the [docs] extra is installed (else the build assertion skips)."""
import importlib.util
import pathlib
import shutil
import subprocess

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
DOCS_URL = "https://marcoarmenta.github.io/quiverlab/"


def test_mkdocs_config_is_canonical():
    y = (ROOT / "mkdocs.yml").read_text()
    assert f"site_url: {DOCS_URL}" in y
    assert "gen_ref_pages.py" in y and "mkdocstrings" in y and "mkdocs-jupyter" in y
    assert "exclude_docs" in y and "plans/" in y            # plan/spec md kept out
    for ch in ("internals/01-exact-fields.md", "tutorials/03-hochschild.ipynb"):
        assert ch in y


def test_gen_ref_script_present():
    assert (ROOT / "scripts" / "gen_ref_pages.py").is_file()


def test_gen_bibliography_reads_packaged_bib():
    src = (ROOT / "scripts" / "gen_bibliography.py").read_text()
    assert "references_bib_path" in src and "bibliography.md" in src
    # the packaged bib is the source of truth, never a docs/-tree copy
    assert "docs/references.bib" not in (ROOT / "mkdocs.yml").read_text()


def test_nav_covers_internals_and_tutorials():
    """Every internals chapter and tutorial notebook must appear in the mkdocs nav,
    so --strict (omitted_files -> error) never trips when Plans 04-07 add chapters."""
    y = (ROOT / "mkdocs.yml").read_text()
    chapters = sorted((ROOT / "docs" / "internals").glob("[0-9][0-9]-*.md"))
    notebooks = sorted((ROOT / "docs" / "tutorials").glob("*.ipynb"))
    missing = [f"internals/{p.name}" for p in chapters if f"internals/{p.name}" not in y]
    missing += [f"tutorials/{p.name}" for p in notebooks if f"tutorials/{p.name}" not in y]
    assert not missing, f"add these to mkdocs.yml nav (else --strict fails): {missing}"


def test_mkdocs_builds_strict_when_available():
    if importlib.util.find_spec("mkdocs") is None or shutil.which("mkdocs") is None:
        import pytest
        pytest.skip("[docs] extra not installed; the acceptance task runs the real build")
    out = subprocess.run(["mkdocs", "build", "--strict", "-d", "/tmp/quiverlab_site"],
                         cwd=ROOT, capture_output=True, text=True,
                         env={"NUMBA_NUM_THREADS": "2", "OMP_NUM_THREADS": "2",
                              "MPLBACKEND": "Agg",  # headless notebook execution
                              "PATH": __import__("os").environ["PATH"]})
    assert out.returncode == 0, out.stderr[-3000:]
```

Add `test_docs_workflow_deploys_pages` to `tests/release/test_workflows.py`:

```python
def test_docs_workflow_deploys_pages():
    d = _read("docs.yml")
    assert "actions/upload-pages-artifact@v5" in d
    assert "actions/deploy-pages@v5" in d
    assert "mkdocs build --strict" in d
    assert "pages: write" in d and "id-token: write" in d
```

- [ ] **Step 1: Write `tests/release/test_docs.py` + the workflow assertion** (above). Expected FAIL (no mkdocs.yml).
- [ ] **Step 2: Run to verify it fails** — `... -m pytest tests/release/test_docs.py -q`.
- [ ] **Step 3: Write `mkdocs.yml`, `scripts/gen_ref_pages.py`, `scripts/gen_bibliography.py`, `docs/index.md`, `docs/javascripts/mathjax.js`, `.github/workflows/docs.yml`** (above). List **all** existing `docs/internals/NN-*.md` chapters (08-11 if Plans 04-07 delivered them) in the nav — the nav-coverage test enforces it.
- [ ] **Step 4: Install `[docs]` and build for real**
  Run:
  ```bash
  cd /Users/marco/Desktop/HomologicalNetworks/quiverlab
  .venv/bin/python -m pip install -e ".[docs]"
  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/mkdocs build --strict -d /tmp/quiverlab_site 2>&1 | tail -20
  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/release/test_docs.py tests/release/test_workflows.py tests/test_no_floats.py -q
  ```
  Expected: `mkdocs build --strict` exits 0 — the three tutorials execute, the API reference generates, the internals chapters render, no strict warnings. Docs + workflow tests pass. (If a tutorial fails to execute under the current library, fix the notebook or the library — do not set `allow_errors: true`.)
  **Env note:** the docs build needs the `[docs]` toolchain (mkdocs-material, mkdocstrings[python], mkdocs-jupyter, gen-files, literate-nav, section-index, mike) in the venv — `/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pip install -e ".[docs]"` before running `mkdocs build`.
- [ ] **Step 5: Commit**
  ```bash
  git add mkdocs.yml scripts/gen_ref_pages.py scripts/gen_bibliography.py docs/index.md docs/javascripts .github/workflows/docs.yml tests/release/test_docs.py tests/release/test_workflows.py
  git commit -m "docs: mkdocs-material site (auto API ref, executed tutorials, internals, packaged-bib References) + Pages deploy

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
  ```

---

### Task 9: JOSS paper — `paper.md` full draft + single-source bib + draft-PDF workflow

**Files:**
- Create: `paper/paper.md`, `.github/workflows/paper.yml`, `tests/release/test_paper.py`
- Modify: `.gitignore` (add `paper/paper.bib`), `tests/release/test_workflows.py` (add `test_paper_workflow_builds_pdf`)
- Modify: `src/quiverlab/citations/references.bib` (append the four paper-only **software** keys — `qpa`, `gap4`, `sagemath`, the `quiverlab` self-cite — to the **single packaged** canonical bib; do NOT create a second `.bib`)
- Modify: Plan 06's `tests/citations/test_registry.py` bib-count assertion (see "Coverage-test edit" below) — Plan 08 **owns** bumping it when it appends entries.

**Interfaces:**
- Consumes: the packaged canonical `src/quiverlab/citations/references.bib` via `quiverlab.citations.references_bib_path()` / `bibtex()` / `all_keys()` (Plan 06), the honest CS scope (Plan 04), the §2 gap table (spec).
- Produces: a JOSS `paper.md` conforming to the **current (2026) required section set**, citing Plan-06's **real** BibTeX ids, and a workflow that compiles it to PDF. **Bibliography single-source:** `paper.md` sets `bibliography: paper.bib`; the workflow copies `src/quiverlab/citations/references.bib → paper/paper.bib` at build time; `paper/paper.bib` is git-ignored; a test asserts every cited key resolves in the packaged bib (read via `references_bib_path()`).

**JOSS 2026 required sections (research-confirmed), in order:** Summary · Statement of need · State of the field · Software design · Research impact statement · AI usage disclosure · Acknowledgements · References. Length **750–1750 words**. Citations pandoc-style (`[@key]`, `@key`).

`paper/paper.md` (full draft):

```markdown
---
title: "quiverlab: exact Hochschild theory for quivers with relations in Python"
tags:
  - Python
  - representation theory
  - quivers with relations
  - Hochschild cohomology
  - Gerstenhaber bracket
  - homological algebra
authors:
  - name: Marco Armenta
    orcid: 0000-0000-0000-0000        # VERIFY: insert Marco's real ORCID before submission
    corresponding: true
    affiliation: 1
affiliations:
  - index: 1
    name: Independent researcher      # VERIFY: set Marco's affiliation before submission
date: 18 July 2026
bibliography: paper.bib
---

# Summary

`quiverlab` is a pure-Python library for computing with finite-dimensional
associative algebras presented as **quivers with relations**, `A = kQ/I`, over
exact fields. Given a quiver and a list of relation strings, it certifies that the
algebra is finite-dimensional, builds an exact multiplication table, and computes
**Hochschild cohomology and homology** together with their algebraic operations —
the cup product, the Gerstenhaber bracket, and the cap action — as well as cyclic
homology, module Ext, and Cartan/Coxeter invariants. Every number is exact: the
library works over the rationals, over exact subfields of the complex numbers
(algebraic number fields `Q(α)`), and over every finite field `GF(p^n)`, and it
fails loudly on any floating-point input rather than returning an approximation.
It is aimed at research algebraists who barely program: three lines take a user
from a presentation to a certified Hochschild table, and every computation can emit
a human-readable worked-steps document.

# Statement of need

Hochschild cohomology `HH^\bullet(A)`, with its Gerstenhaber algebra structure, is
a central invariant in representation theory and deformation theory, yet it is
remarkably hard to compute by hand beyond the smallest examples. Researchers who
study quivers with relations — the standard presentation of finite-dimensional
algebras — have no installable tool that computes these invariants exactly and
without programming overhead. The existing options each stop short: they either do
not compute Hochschild cohomology at all, do not implement its operations, are not
exact, or require a substantial computer-algebra installation and scripting effort.
`quiverlab` fills this gap. It is `pip install quiverlab` with no external system
dependencies, exposes a flat, discoverable API, and returns certified exact answers
across characteristics — so an algebraist can, for instance, watch the
characteristic-`p` pathology of `k[x]/(x^2)` appear by changing a single argument.

# State of the field

The strongest existing system is **QPA** [@qpa], a mature GAP package for quivers
and path algebras: it constructs `kQ/I` by admissible ideals, computes minimal
module resolutions, module Ext, and Auslander–Reiten theory. QPA ships **no
Hochschild cohomology** — it must be assembled by hand via the enveloping algebra
and module Ext — and **no cup product or Gerstenhaber bracket**; installation
requires GAP, non-trivial for non-programmers, and historically awkward on Windows.
**SageMath** [@sagemath] provides only the *free* path algebra, with no
quotient-by-relations object, and an unreduced-bar Hochschild complex usable only at
toy sizes. **Magma**, **Macaulay2/Singular**, and **QuiverTools** address adjacent
problems (Ext algebras, noncommutative Gröbner bases, moduli of representations) but
none computes finite-dimensional Hochschild theory with its operations. On PyPI there
is **nothing**: no quivers-with-relations, no Hochschild cohomology, no bracket. In
short, no system on earth ships Hochschild cohomology *with its operations* for
finite-dimensional algebras, none ships the Chouhy–Solotar resolution, and nothing is
pip-installable for non-coders. `quiverlab` is built to be that system rather than a
contribution to any of the above, because its exact-only, non-coder-first design and
its resolution engine differ fundamentally from each existing architecture.

# Software design

`quiverlab` is layered. A `Domain` protocol carries all coefficient arithmetic
exactly (rationals with fraction-free elimination, number fields via `sympy`, finite
fields via a fast integer kernel with a pure-Python fallback), so no engine has a
floating-point code path; a static analysis gate forbids float literals in the source
tree. On top of it, a quiver-with-relations front end runs an exact noncommutative
**Gröbner (Buchberger–Mora overlap) completion** with a degree bound and an
**admissibility certificate**: it either returns a certified finite-dimensional
algebra with an irreducible-path basis, or fails loudly with the offending cycle —
never a hang and never a guess. The resolution layer offers four interchangeable
backends (normalized bar, minimal `A^e` syzygy, Bardzell for monomial algebras, and
the general **Chouhy–Solotar** resolution [@ChouhySolotar2015], the first full
implementation in any system, specializing exactly to Bardzell [@Bardzell1997] in the
monomial case). This general resolution is certified for quadratic-tip and all
monomial presentations; a non-quadratic non-monomial presentation raises at an
explicit boundary rather than risk a wrong answer. Deep Hochschild *dimensions* are
certified for closed-form families to arbitrary degree and, for general admissible
`kQ/I`, are computed and gated per instance by three independent checks (`d∘d = 0`, an
order condition, and degreewise agreement with the bar/minimal engines within their
window); the Gerstenhaber *operations* [@NegronWitherspoon2016; @Volkov2019] are
transported to bar cochains and are certified inside the bar-buildable degree window.
Module Ext uses minimal projective resolutions in the Green–Solberg–Zacharia style
[@GSZ2001]. An optional `pip install quiverlab[qpa]` backend runs
an **independent** recomputation of module Ext and Hochschild dimensions in QPA (via
the enveloping-algebra route) as a validation oracle in continuous integration.

# Research impact statement

`quiverlab` makes exact Hochschild computations routine where they were previously
manual, enabling systematic study of Gerstenhaber structure across families of
algebras and across characteristics — for example reproducing published Hochschild
dimensions of quantum complete intersections [@BGMS2005] and the vanishing behaviour
of hereditary algebras [@Happel1989], then extending them along parameter sweeps that
would be impractical by hand. Because it is exact and pip-installable, it lowers the
barrier to reproducible experiments in representation theory and provides a common,
citable reference implementation of the Chouhy–Solotar resolution against which future
work can be checked. The library unifies and generalizes the author's prior research
software into a tool independent of the application that produced it.

# AI usage disclosure

The library's design and mathematics are the author's. AI coding assistants were used
under close human supervision for implementation scaffolding, test authoring, and
documentation drafting; every mathematical claim, algorithm, and certified value was
verified by the author against hand computations, published literature, and an
independent computer-algebra cross-check (QPA). No result reported by the software or
this paper rests on unverified AI output.

# Acknowledgements

We thank the developers of QPA [@qpa], GAP [@gap4], SymPy, and NumPy, on whose ideas
and tools this work builds. *(Financial support to be acknowledged before submission.)*

# References
```

**Word count:** the body above is ~760–810 words of prose (front matter, headings, and
fenced code excluded) — comfortably within JOSS's 750–1750 band; the word-count test
flags if editing pushes it out of range.

**Bib keys used** (all must resolve in the **packaged** `src/quiverlab/citations/references.bib`,
Plan 06's **real** BibTeX ids): `@ChouhySolotar2015`, `@Bardzell1997`,
`@NegronWitherspoon2016`, `@Volkov2019` (note the year — **not** 2016), `@GSZ2001`,
`@BGMS2005`, `@Happel1989` (Plan 06 seeds these 14 entries; check
`scratchpad/plan-06-families-draft.md` — it may also gain `CartanEilenberg1956` and a
Gerstenhaber–Schack/incidence key in its own fix round, which does not affect us since
we read the live bib). **This task appends four software keys** — `@qpa`, `@gap4`,
`@sagemath`, and the `@software{quiverlab}` self-citation — to the **same packaged**
bib. Entries to append verbatim:

```bibtex
@misc{qpa,
  author = {Green, Edward L. and Solberg, {\O}yvind},
  title  = {{QPA} -- {Quivers}, path algebras and representations, a {GAP} package, Version 1.37},
  year   = {2026}, howpublished = {\url{https://folk.ntnu.no/oyvinso/QPA/}}
}
@manual{gap4,
  key = {GAP}, organization = {The GAP Group},
  title = {{GAP -- Groups, Algorithms, and Programming}, Version 4.13}, year = {2026},
  note = {\url{https://www.gap-system.org}}
}
@misc{sagemath,
  author = {{The Sage Developers}}, title = {{SageMath}, the {Sage} {Mathematics} {Software} {System}},
  year = {2026}, note = {\url{https://www.sagemath.org}}
}
@software{quiverlab,
  author = {Armenta, Marco}, title = {quiverlab: exact Hochschild theory for quivers with relations},
  year = {2026}, version = {0.1.0}, url = {https://github.com/MarcoArmenta/quiverlab}
}
```

**Coverage-test edit (Plan 08 owns this).** Plan 06's `tests/citations/test_registry.py`
carries the bib-count assertion — being converted (coordinated with the P06 author) from
`count == 14` to a **coverage** assertion: *every registry key resolves to a `@`-entry in
`references_bib_path()`* AND *the total entry count `>= N`*. Appending these four entries
takes the file from **16 → 20** entries (Plan 06's final draft ships 16: its fix round added Hochschild1945 and CartanEilenberg1956), all four being **registry-less software keys**.
So the Plan-08 edit is: (1) bump the floor `N` from `16` to `20`; (2) if that test also
forbids orphan entries (every `@`-id must be a registry `bibtex_key`), add
`{"qpa", "gap4", "sagemath", "quiverlab"}` to its allowed-software-keys set. The
"all registry keys resolve" half stays green untouched (the 14 registry `bibtex_key`s are
not modified). State this in the commit body.

`.github/workflows/paper.yml`:

```yaml
name: JOSS draft PDF

on:
  push:
    paths:
      - "paper/**"
      - "src/quiverlab/citations/references.bib"     # the single packaged bib
      - ".github/workflows/paper.yml"
  workflow_dispatch:

jobs:
  paper:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v7
      - name: Assemble single-source bibliography
        # one bib, packaged; paper.bib is a build-time copy (gitignored)
        run: cp src/quiverlab/citations/references.bib paper/paper.bib
      - name: Build draft PDF
        uses: openjournals/openjournals-draft-action@master   # moving ref; pin a SHA if needed
        with:
          journal: joss
          paper-path: paper/paper.md
      - uses: actions/upload-artifact@v7
        with:
          name: paper
          path: paper/paper.pdf
```

Add to `.gitignore`: `paper/paper.bib`.

`tests/release/test_paper.py`:

```python
"""paper.md conforms to JOSS: required sections present, word count in band, and every
cited @key resolves in the single PACKAGED bib (Plan 08 Task 9), read via the library's
own accessor (never a docs/ path)."""
import pathlib
import re

from quiverlab.citations import references_bib_path

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
PAPER = (ROOT / "paper" / "paper.md").read_text()
BIB = references_bib_path().read_text(encoding="utf-8")

_REQUIRED = ["# Summary", "# Statement of need", "# State of the field",
             "# Software design", "# Research impact statement",
             "# AI usage disclosure", "# Acknowledgements", "# References"]


def _body(text):
    # drop the YAML front matter (it contains an @-email that is not a citation)
    return text.split("---", 2)[-1] if text.lstrip().startswith("---") else text


def test_required_joss_sections_present_in_order():
    idx = [PAPER.find(h) for h in _REQUIRED]
    assert all(i >= 0 for i in idx), [h for h, i in zip(_REQUIRED, idx) if i < 0]
    assert idx == sorted(idx), "JOSS sections out of order"


def test_word_count_in_joss_band():
    body = re.sub(r"```.*?```", "", _body(PAPER), flags=re.S)
    words = len(re.findall(r"\S+", body))
    assert 750 <= words <= 1750, f"paper body is {words} words (JOSS wants 750-1750)"


def test_every_citation_key_resolves_in_packaged_bib():
    keys = set(re.findall(r"@([A-Za-z][A-Za-z0-9_]+)", _body(PAPER)))   # body only, not the email
    bib_keys = set(re.findall(r"@\w+\{([^,]+),", BIB))
    missing = sorted(keys - bib_keys)
    assert not missing, f"paper cites keys absent from the packaged references.bib: {missing}"


def test_paper_bib_is_not_committed_separately():
    assert not (ROOT / "paper" / "paper.bib").exists() or \
        "paper/paper.bib" in (ROOT / ".gitignore").read_text()
```

Add to `tests/release/test_workflows.py`:

```python
def test_paper_workflow_builds_pdf():
    p = _read("paper.yml")
    assert "openjournals/openjournals-draft-action@master" in p
    assert "paper-path: paper/paper.md" in p
    assert "cp src/quiverlab/citations/references.bib paper/paper.bib" in p   # packaged single-source bib
```

- [ ] **Step 1: Write `tests/release/test_paper.py` + the workflow assertion** (above). Expected FAIL (no paper).
- [ ] **Step 2: Run to verify it fails** — `... -m pytest tests/release/test_paper.py -q`.
- [ ] **Step 3: Write `paper/paper.md`, `.github/workflows/paper.yml`, the `.gitignore` line; append the four software keys to `src/quiverlab/citations/references.bib`; and bump Plan 06's coverage-test floor (16 → 20)** (all above). Insert Marco's real ORCID/affiliation (flagged) — the test does not require them, but they must be set before submission.
- [ ] **Step 4: Run**
  Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/release/test_paper.py tests/citations/test_registry.py tests/release/test_workflows.py tests/test_no_floats.py -q`
  Expected: pass — sections present/ordered, 750–1750 words, every `@key` resolves in the packaged bib, `paper.bib` not committed, and Plan 06's citation registry/coverage tests stay green with the bumped floor. (Real PDF compilation happens in `paper.yml` on GitHub via the `inara`/pandoc container; locally we validate structure and citations.)
- [ ] **Step 5: Commit**
  ```bash
  git add paper/paper.md .github/workflows/paper.yml .gitignore src/quiverlab/citations/references.bib tests/citations/test_registry.py tests/release/test_paper.py tests/release/test_workflows.py
  git commit -m "docs(paper): JOSS paper.md draft (2026 sections) + draft-PDF workflow; append 4 software keys to packaged bib

Appends qpa/gap4/sagemath/quiverlab to src/quiverlab/citations/references.bib
(16 -> 20 entries) and bumps Plan 06's coverage floor; paper.bib is a build-time
copy of the single packaged bib.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
  ```

---

### Task 10: PyPI packaging + release workflow (trusted publishing) — no real upload

**Files:**
- Create: `.github/workflows/release.yml`, `CHANGELOG.md`, `scripts/check_pypi_name.py`, `tests/release/test_build.py`
- Modify: `tests/release/test_workflows.py` (add `test_release_workflow_trusted_publishing`)

**Interfaces:**
- Consumes: the modernized `pyproject.toml` (Task 3), `build`/`twine` (`[dev]` extra).
- Produces: a tag-gated OIDC **trusted-publishing** workflow (no API token ever committed), a name-availability check, a local build + `twine check`, and the **sdist/wheel data-file packaging check (open item #4)** — `test_build.py::test_wheel_and_sdist_contain_package_data` asserts the built wheel and sdist actually CONTAIN `quiverlab/citations/references.bib` and `quiverlab/families/zoo_catalog.json` (a `twine check` PASS does not verify archive contents, so a dropped `[tool.setuptools.package-data]` block would otherwise ship silently and break `bibliography()`/`zoo()`). **Nothing is uploaded to real PyPI during this plan** — the workflow triggers only on a pushed `v*` tag, which is Marco's action.

**Name availability (research, point-in-time 2026-07-18):** both **`quiverlab`** and **`quiver-lab`** return HTTP 404 on the PyPI JSON API → **both AVAILABLE** (they are distinct under PEP 503 normalization). D8: register **`quiverlab`** (primary); optionally also register `quiver-lab` to defend the namespace. Re-check immediately before registering — availability is point-in-time. Use the **JSON API** (`https://pypi.org/pypi/<name>/json`), not the project-page URL (Cloudflare-challenges `curl`, returning 200 for every name).

`scripts/check_pypi_name.py`:

```python
"""Advisory: is a name free on PyPI? 404 on the JSON API == available.
(The /project/<name>/ page is Cloudflare-challenged and unreliable via curl.)"""
import urllib.error
import urllib.request


def available(name: str):
    try:
        urllib.request.urlopen(f"https://pypi.org/pypi/{name}/json", timeout=10)
        return False                       # 200 -> taken
    except urllib.error.HTTPError as e:
        return True if e.code == 404 else None
    except Exception:                      # noqa: BLE001  (offline -> unknown)
        return None


if __name__ == "__main__":
    for n in ("quiverlab", "quiver-lab"):
        a = available(n)
        print(f"{n}: {'AVAILABLE' if a else 'TAKEN' if a is False else 'UNKNOWN (offline?)'}")
```

`.github/workflows/release.yml` (trusted publishing; tag-gated; version/tag consistency guard):

```yaml
name: Release

on:
  push:
    tags: ["v*"]        # ONLY on a version tag Marco pushes; never on branch pushes

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v7
      - uses: actions/setup-python@v6
        with:
          python-version: "3.12"
      - name: Tag matches pyproject version
        run: |
          TAG="${GITHUB_REF_NAME#v}"
          VER=$(python -c "import tomllib; print(tomllib.load(open('pyproject.toml','rb'))['project']['version'])")
          test "$TAG" = "$VER" || { echo "tag $TAG != pyproject version $VER"; exit 1; }
      - run: python -m pip install --upgrade pip build twine
      - run: python -m build                       # sdist + wheel into dist/
      - run: python -m twine check dist/*          # metadata/readme render check
      - uses: actions/upload-artifact@v7
        with:
          name: dist
          path: dist/

  publish:
    needs: build
    runs-on: ubuntu-latest
    environment:
      name: pypi
      url: https://pypi.org/p/quiverlab
    permissions:
      id-token: write                              # OIDC trusted publishing (no token)
    steps:
      - uses: actions/download-artifact@v8
        with:
          name: dist
          path: dist/
      - uses: pypa/gh-action-pypi-publish@release/v1
```

**PyPI-side setup Marco does once (after this plan lands, before the first tag):** on pypi.org → Account settings → Publishing → *Add a pending publisher (GitHub)*: PyPI project name `quiverlab`, owner `MarcoArmenta`, repository `quiverlab`, workflow filename `release.yml`, environment name `pypi`. (Optionally mirror on test.pypi.org with environment `testpypi` for a dry run.) After the first successful publish the pending publisher becomes permanent. **No secrets are added to the repo.**

`CHANGELOG.md` (seed):

```markdown
# Changelog

All notable changes to quiverlab are documented here. This project adheres to
[Semantic Versioning](https://semver.org) (0.x during battle-testing; 1.0 at JOSS
acceptance).

## [Unreleased]
### Added
- Optional `[qpa]` backend: `A.crosscheck(...)` (independent QPA recomputation).
- GitHub Actions CI (matrix + engine-path legs), docs site, JOSS paper draft.
- Modernized packaging (PEP 639 SPDX license), README, community files.
```

`tests/release/test_build.py`:

```python
"""The wheel/sdist build cleanly and pass twine check (Plan 08 Task 10). Skips if
build/twine are not installed; the acceptance task installs [dev] and runs it."""
import importlib.util
import pathlib
import subprocess
import sys
import tarfile
import zipfile

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent


def _have(mod):
    return importlib.util.find_spec(mod) is not None


def test_build_and_twine_check(tmp_path):
    if not (_have("build") and _have("twine")):
        import pytest
        pytest.skip("build/twine not installed; acceptance task runs the real build")
    out = tmp_path / "dist"
    r = subprocess.run([sys.executable, "-m", "build", "--outdir", str(out)],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    dists = list(out.glob("quiverlab-*.whl")) + list(out.glob("quiverlab-*.tar.gz"))
    assert len(dists) == 2, [p.name for p in dists]
    r2 = subprocess.run([sys.executable, "-m", "twine", "check", *map(str, out.glob("*"))],
                        capture_output=True, text=True)
    assert r2.returncode == 0, r2.stdout + r2.stderr
    assert "PASSED" in r2.stdout


# Data-file (package-data) files loaded at runtime by __file__-relative path.
_DATA_FILES = ("quiverlab/citations/references.bib",
               "quiverlab/families/zoo_catalog.json")


def test_wheel_and_sdist_contain_package_data(tmp_path):
    """sdist/wheel data-file packaging check (open item #4).

    The BUILT wheel and sdist must CONTAIN the packaged data files that
    citations/registry.py (references.bib) and families/zoo.py (zoo_catalog.json)
    load by __file__-relative path -- otherwise bibliography()/zoo() break in an
    installed environment even though `twine check` still PASSES (twine checks
    metadata/readme rendering, never archive contents). This test FAILS if the
    [tool.setuptools.package-data] block is removed from pyproject.toml."""
    if not _have("build"):
        import pytest
        pytest.skip("build not installed; acceptance task runs the real build")
    out = tmp_path / "dist"
    r = subprocess.run([sys.executable, "-m", "build", "--outdir", str(out)],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]

    # wheel is a zip; package-data lands at quiverlab/<subpkg>/<file> (no src/ prefix).
    wheels = list(out.glob("quiverlab-*.whl"))
    assert wheels, [p.name for p in out.glob("*")]
    with zipfile.ZipFile(wheels[0]) as zf:
        wheel_names = set(zf.namelist())
    for rel in _DATA_FILES:
        assert rel in wheel_names, f"{rel} missing from the wheel (package-data dropped?)"

    # sdist is a tar.gz; entries are under <name-version>/src/quiverlab/...
    sdists = list(out.glob("quiverlab-*.tar.gz"))
    assert sdists, [p.name for p in out.glob("*")]
    with tarfile.open(sdists[0]) as tf:
        sdist_names = tf.getnames()
    for rel in _DATA_FILES:
        assert any(m.endswith(rel) for m in sdist_names), \
            f"{rel} missing from the sdist (package-data dropped?)"
```

Add to `tests/release/test_workflows.py`:

```python
def test_release_workflow_trusted_publishing():
    r = _read("release.yml")
    assert 'tags: ["v*"]' in r                     # tag-gated only
    assert "pypa/gh-action-pypi-publish@release/v1" in r
    assert "id-token: write" in r                  # OIDC, no token
    assert "password:" not in r                    # no committed secret
    assert "twine check dist/*" in r
```

- [ ] **Step 1: Write `tests/release/test_build.py` + the workflow assertion** (above). Expected FAIL (no `release.yml`).
- [ ] **Step 2: Run to verify it fails** — `... -m pytest tests/release/test_build.py tests/release/test_workflows.py -q`.
- [ ] **Step 3: Write `release.yml`, `CHANGELOG.md`, `scripts/check_pypi_name.py`** (above).
- [ ] **Step 4: Build locally + check name + run**
  Run:
  ```bash
  cd /Users/marco/Desktop/HomologicalNetworks/quiverlab
  .venv/bin/python -m pip install -q build twine
  .venv/bin/python scripts/check_pypi_name.py            # advisory (needs network)
  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/release/test_build.py tests/release/test_workflows.py tests/test_no_floats.py -q
  ```
  Expected: `python -m build` produces `quiverlab-<ver>.tar.gz` + `.whl`; `twine check` reports `PASSED`; the package-data check confirms `references.bib` + `zoo_catalog.json` are inside both artifacts (open item #4); the workflow test confirms trusted publishing + tag gating + no token. **No upload happens.**
  **Env note:** this step needs `build` + `twine` in the venv — `/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pip install build twine` (both are in the `[dev]` extra).
- [ ] **Step 5: Commit**
  ```bash
  git add .github/workflows/release.yml CHANGELOG.md scripts/check_pypi_name.py tests/release/test_build.py tests/release/test_workflows.py
  git commit -m "build(release): tag-gated trusted-publishing workflow + build/twine check + name check

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
  ```

---

### Task 11: README overhaul — the front door

**Files:**
- Modify: `README.md`
- Create: `tests/release/test_readme.py`

**Interfaces:**
- Consumes: the docs URL, repo URL, marker/CI names.
- Produces: a release front door — badges, the spec §1 three-line quickstart, and links to the docs site, web GUI, tutorials, internals, and JOSS. Preserves the existing status/examples body (updated: real `pip install quiverlab`, not the dev-preview line).

Replace the top of `README.md` (through the install section) with:

```markdown
# quiverlab

[![CI](https://github.com/MarcoArmenta/quiverlab/actions/workflows/ci.yml/badge.svg)](https://github.com/MarcoArmenta/quiverlab/actions/workflows/ci.yml)
[![Docs](https://github.com/MarcoArmenta/quiverlab/actions/workflows/docs.yml/badge.svg)](https://marcoarmenta.github.io/quiverlab/)
[![PyPI](https://img.shields.io/pypi/v/quiverlab.svg)](https://pypi.org/project/quiverlab/)
[![Python](https://img.shields.io/pypi/pyversions/quiverlab.svg)](https://pypi.org/project/quiverlab/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
<!-- [![DOI](https://joss.theoj.org/papers/<id>/status.svg)](https://doi.org/<doi>) -- added at JOSS acceptance -->

**Quivers with relations and Hochschild theory, exactly, for algebraists.**

quiverlab computes with finite-dimensional algebras `kQ/I` over the complex numbers
(exactly — no floating point, ever) and over all finite fields: certified
finite-dimensionality, Hochschild (co)homology with cup products and Gerstenhaber
brackets, the first full Chouhy–Solotar resolution, module Ext, and Cartan/Coxeter
invariants. Floats fail loudly by design.

## Install

```bash
pip install quiverlab                 # pure-Python core, no external systems
pip install "quiverlab[fast]"         # + numba GF(p) acceleration (optional)
pip install "quiverlab[qpa]"          # + GAP/QPA cross-check backend (macOS/Linux)
```

## Three lines to a Hochschild table

```python
from quiverlab import Quiver, CC

Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})
print(Q.algebra(relations=["a*b"], field=CC).hochschild_cohomology(3))
```

## Learn more

- **Documentation:** <https://marcoarmenta.github.io/quiverlab/>
- **Tutorials:** [executable notebooks](docs/tutorials/) — start here.
- **Under the hood:** [internals chapters](docs/internals/) — how each number is produced.
- **Web GUI:** compute without installing (Plan 09; linked from the docs when live).
- **Cite:** see the JOSS paper (`paper/paper.md`) and [`CITATION.cff`](CITATION.cff).
```

Keep the existing "The classic characteristic pathology", "Status", and "General
quivers with relations (kQ/I)" sections below (they are current and good); just
ensure the install snippet there reads `pip install quiverlab`, not the dev preview.

`tests/release/test_readme.py`:

```python
"""README is the release front door: badges, real install, quickstart, links."""
import pathlib

README = (pathlib.Path(__file__).resolve().parent.parent.parent / "README.md").read_text()


def test_badges_present():
    for b in ("actions/workflows/ci.yml/badge.svg", "img.shields.io/pypi/v/quiverlab",
              "License-MIT"):
        assert b in README


def test_install_and_quickstart():
    assert "pip install quiverlab" in README
    assert 'pip install "quiverlab[qpa]"' in README
    assert "hochschild_cohomology(3)" in README


def test_links_to_docs_and_tutorials_and_citation():
    assert "https://marcoarmenta.github.io/quiverlab/" in README
    assert "docs/tutorials" in README and "docs/internals" in README
    assert "CITATION.cff" in README
```

- [ ] **Step 1: Write `tests/release/test_readme.py`** (above). Expected FAIL (no badges/links yet).
- [ ] **Step 2: Run to verify it fails** — `... -m pytest tests/release/test_readme.py -q`.
- [ ] **Step 3: Rewrite the README top** (above), preserving the good body sections.
- [ ] **Step 4: Run**
  Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/release/test_readme.py tests/test_quickstart.py tests/test_no_floats.py -q`
  Expected: pass (README front door + the existing quickstart example still runs).
- [ ] **Step 5: Commit**
  ```bash
  git add README.md tests/release/test_readme.py
  git commit -m "docs(readme): release front door — badges, quickstart, docs/web/JOSS links

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
  ```

---

### Task 12: Community + citation files; docs Development page; internals-chapter decision

**Files:**
- Create: `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `CITATION.cff`, `docs/development/release.md`
- Modify: `docs/internals/README.md` (honest-coverage update)
- Create: `tests/release/test_community.py`

**Interfaces:**
- Consumes: the docs URL, repo URL, the QPA cross-check route (Task 4).
- Produces: the JOSS-required community artifacts (contribution guidelines, code of conduct, machine-readable citation), the docs-site Development/Release page (where release infrastructure **and** the QPA cross-check *mathematics* are documented), and the URL-consistency lock.

**Internals-chapter decision (binding, justified): DO NOT add `docs/internals/12-release-infrastructure.md`.** The `docs/internals/` tree has a stated audience and purpose — "how quiverlab represents its objects and produces each number, for algebraists who do not program" (internals README). CI, packaging, PyPI, and the JOSS paper are release *infrastructure*, not core mathematical representation; a release-infra chapter would violate the tree's contract. The roadmap's standing constraint ("each plan adds internals chapters for what it introduces") is satisfied differently for Plan 08: the one genuinely *mathematical* thing Plan 08 introduces is the **QPA cross-check route** (`HH^n(A) = Ext^n_{A^e}(A,A)` via the enveloping algebra — QPA ships no HH function), and that is documented on the **docs Development/Release page**, which the internals README points to. Rationale recorded so a future reader does not "fix" the missing chapter.

`CITATION.cff`:

```yaml
cff-version: 1.2.0
message: "If you use quiverlab in your research, please cite it as below."
title: "quiverlab: exact Hochschild theory for quivers with relations"
type: software
authors:
  - family-names: Armenta
    given-names: Marco
    email: drmarcoarmenta@gmail.com
    orcid: "https://orcid.org/0000-0000-0000-0000"   # VERIFY: real ORCID before release
version: 0.1.0
date-released: "2026-07-18"
license: MIT
repository-code: "https://github.com/MarcoArmenta/quiverlab"
url: "https://marcoarmenta.github.io/quiverlab/"
keywords:
  - quivers with relations
  - Hochschild cohomology
  - Gerstenhaber bracket
  - representation theory
  - exact arithmetic
# At JOSS acceptance, add the preferred-citation (type: article) with the JOSS DOI:
# preferred-citation:
#   type: article
#   title: "quiverlab: exact Hochschild theory for quivers with relations"
#   authors: [{family-names: Armenta, given-names: Marco}]
#   journal: "Journal of Open Source Software"
#   doi: "10.21105/joss.XXXXX"
#   year: 2026
```

`CONTRIBUTING.md`:

```markdown
# Contributing to quiverlab

Thank you for your interest! quiverlab is a pure-Python library for exact
computation with finite-dimensional algebras.

## Development setup

```bash
git clone https://github.com/MarcoArmenta/quiverlab
cd quiverlab
python3.12 -m venv .venv && . .venv/bin/activate
pip install -e ".[dev,fast]"          # add ,docs for the docs site; ,qpa on macOS/Linux
```

## Running tests

The suite is split by marker (see `pyproject.toml`):

```bash
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 pytest -m fast     # quick, cross-platform
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 pytest             # full (adds -m deep)
NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 pytest -m qpa      # needs [qpa] + GAP (macOS/Linux)
QUIVERLAB_NO_NUMBA=1 pytest                              # the pure-Python engine path
```

Always throttle threads with `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2`.

## Non-negotiables

- **Exact arithmetic only.** No float or complex literals anywhere under `src/`;
  the AST gate `tests/test_no_floats.py` enforces this and must stay green.
- **Loud failure over silent approximation.** Errors carry a fix-it hint.
- **Tests first.** Add a failing test, make it pass, keep the full suite green.
- **Conventional commits** (`feat:`, `fix:`, `test:`, `docs:`, `ci:`, `build:`).

## Reporting problems / asking for help

Open a GitHub issue. For the web GUI there is also an in-app feedback form.
```

`CODE_OF_CONDUCT.md`:

```markdown
# Code of Conduct

This project adopts the [Contributor Covenant](https://www.contributor-covenant.org),
version 2.1. In short: be respectful, welcoming, and constructive. Harassment and
exclusionary behavior are not tolerated.

Report concerns to Marco Armenta at drmarcoarmenta@gmail.com. Reports are handled
confidentially. The full text is at
<https://www.contributor-covenant.org/version/2/1/code_of_conduct/>.
```

`docs/development/release.md` (docs-site page — CI/packaging/JOSS **and** the QPA cross-check mathematics):

```markdown
# Development, release, and the QPA cross-check

## Continuous integration

- **`ci.yml`** — a fast test suite on every OS × Python 3.10–3.13 cell, plus two
  deep Linux legs running the full suite on the numba and pure-Python engine paths
  (`QUIVERLAB_NO_NUMBA=1`), plus a float-ban lint. Tests are bucketed by the
  `fast` / `deep` / `qpa` markers.
- **`qpa.yml`** — a Linux-only job that installs GAP + QPA (`passagemath-gap[qpa]`)
  and runs the cross-check suite (`-m qpa`), mandatory there.
- **`docs.yml`** — builds this site `--strict` (executing the tutorials) and deploys
  to GitHub Pages.
- **`paper.yml`** — compiles the JOSS paper draft to PDF.
- **`release.yml`** — on a `v*` tag: build, `twine check`, and publish to PyPI via
  OIDC trusted publishing (no API token).

## Releasing (semver 0.x → 1.0 at JOSS acceptance)

1. Bump `version` in `pyproject.toml` (and `__version__`); update `CHANGELOG.md`.
2. Commit, then `git tag vX.Y.Z && git push --tags`.
3. `release.yml` builds and publishes; the tag must equal the pyproject version.

## The QPA cross-check (the mathematics)

QPA has no Hochschild cohomology function, so `quiverlab[qpa]` assembles it via the
**enveloping algebra**. For `A = kQ/I` with enveloping algebra `A^e = A^{op} ⊗ A`,
Hochschild cohomology is `A^e`-module Ext of `A` with itself:

$$ HH^n(A) \;=\; \mathrm{Ext}^n_{A^e}(A, A). $$

`A.crosscheck("hochschild", n)` scripts, in GAP: build the quiver and `PathAlgebra`
over the same field, form `EnvelopingAlgebra(A)`, present `A` as a right `A^e`-module,
take a minimal projective resolution, and read `dim Ext^k` for `k = 0..n`. It then
compares these to quiverlab's own `hochschild_cohomology(n).dims` and fails loudly on
any disagreement (both use QPA's `ExtAlgebraGenerators(-, n)[1]` dimension series).
`A.crosscheck("module_ext", M, n)` does the analogous check for module **self-Ext**
`Ext^*(M, M)`; distinct-module `Ext(M, N)` (needing `ExtOverAlgebra` + iterated
syzygies) is a flagged post-v1 extension. This is an independent-implementation oracle
(spec §8 ring 3), not a dependency of the core.
```

Append to `docs/internals/README.md` (honest-coverage), a short note:

```markdown
## Release infrastructure and the QPA cross-check

Plan 08 (release) adds **no internals chapter**: CI, packaging, the docs site, and
the JOSS paper are infrastructure, not core mathematics. The one mathematical addition
— the optional QPA cross-check that recomputes Hochschild dimensions via the
enveloping algebra `HH^n(A) = Ext^n_{A^e}(A,A)` — is documented on the docs
[Development page](../development/release.md), not here.
```

`tests/release/test_community.py`:

```python
"""Community/citation files exist and the canonical docs URL is consistent across
mkdocs.yml, README, and CITATION.cff (Plan 08 Task 12)."""
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
URL = "https://marcoarmenta.github.io/quiverlab/"


def test_community_files_present():
    for f in ("CONTRIBUTING.md", "CODE_OF_CONDUCT.md", "CITATION.cff"):
        assert (ROOT / f).is_file(), f"missing {f} (JOSS review checklist)"


def test_citation_cff_shape():
    cff = (ROOT / "CITATION.cff").read_text()
    assert "cff-version: 1.2.0" in cff and "Armenta" in cff
    assert 'repository-code: "https://github.com/MarcoArmenta/quiverlab"' in cff


def test_citation_version_matches_pyproject():
    import re
    import tomllib
    ver = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
    core = re.match(r"\d+\.\d+\.\d+", ver).group(0)        # 0.1.0.dev0 -> 0.1.0
    cff = (ROOT / "CITATION.cff").read_text()
    assert f"version: {core}" in cff, f"CITATION.cff version must be the release core {core}"


def test_docs_url_consistent_across_files():
    mk = (ROOT / "mkdocs.yml").read_text()
    rd = (ROOT / "README.md").read_text()
    cff = (ROOT / "CITATION.cff").read_text()
    assert f"site_url: {URL}" in mk
    assert URL in rd
    assert f'url: "{URL}"' in cff


def test_development_page_documents_qpa_math():
    dev = (ROOT / "docs" / "development" / "release.md").read_text()
    assert "Ext^n_{A^e}(A, A)" in dev and "EnvelopingAlgebra" in dev
```

- [ ] **Step 1: Write `tests/release/test_community.py`** (above). Expected FAIL.
- [ ] **Step 2: Run to verify it fails** — `... -m pytest tests/release/test_community.py -q`.
- [ ] **Step 3: Write the four files + the internals README note** (above).
- [ ] **Step 4: Run**
  Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/release/test_community.py tests/release/test_docs.py tests/test_no_floats.py -q`
  Expected: pass — files present, CITATION shape, URL consistent across mkdocs/README/CITATION, Development page carries the enveloping-algebra formula. (Re-run `mkdocs build --strict` if `[docs]` is installed to confirm the new `development/release.md` page resolves in the nav.)
- [ ] **Step 5: Commit**
  ```bash
  git add CONTRIBUTING.md CODE_OF_CONDUCT.md CITATION.cff docs/development/release.md docs/internals/README.md tests/release/test_community.py
  git commit -m "docs: community + citation files, docs Development page (QPA cross-check math); skip release internals chapter

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
  ```

---

### Task 13: Release-candidate acceptance — everything green, docs live, paper compiles, fresh-venv install

**Files:**
- Create: `tests/release/test_acceptance.py`

**Interfaces:**
- Consumes: the entire release surface (Tasks 1–12).
- Produces: the release-candidate gate. This is the only task that runs **one tracked background full suite**, awaits it, then commits and reports.

`tests/release/test_acceptance.py`:

```python
"""Plan 08 acceptance: the whole release surface is present and coherent."""
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent


def test_all_workflows_present():
    wf = ROOT / ".github" / "workflows"
    for name in ("ci.yml", "qpa.yml", "docs.yml", "paper.yml", "release.yml"):
        assert (wf / name).is_file(), f"missing workflow {name}"


def test_release_artifacts_present():
    for rel in ("mkdocs.yml", "paper/paper.md", "CITATION.cff", "CONTRIBUTING.md",
                "CODE_OF_CONDUCT.md", "CHANGELOG.md", "scripts/gen_ref_pages.py",
                "scripts/release_freshness.py"):
        assert (ROOT / rel).is_file(), f"missing {rel}"


def test_pyproject_release_ready():
    import tomllib
    pp = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert pp["project"]["license"] == "MIT"                 # SPDX, not table
    extras = pp["project"]["optional-dependencies"]
    assert {"fast", "qpa", "docs", "dev"} <= set(extras)


def test_qpa_backend_optional_not_imported_by_core():
    import quiverlab                                          # must import w/o GAP
    assert hasattr(quiverlab.Algebra, "crosscheck")
```

- [ ] **Step 1: Write `tests/release/test_acceptance.py`** (above).
- [ ] **Step 2: Run the acceptance test + the release runbook**
  **Env note:** the runbook needs `build`, `twine`, and the `[docs]` toolchain (mkdocs) in the venv — `/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pip install -e ".[dev,docs]"` installs all three (build/twine from `[dev]`, mkdocs stack from `[docs]`) before running steps (c)–(e).
  ```bash
  cd /Users/marco/Desktop/HomologicalNetworks/quiverlab
  # a) freshness gate is green (prerequisites intact)
  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python scripts/release_freshness.py

  # b) acceptance + all release tests
  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/release -q

  # c) docs build live (strict; executes tutorials; builds API ref)
  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/mkdocs build --strict -d /tmp/quiverlab_site

  # d) build the distributions + twine check
  .venv/bin/python -m build && .venv/bin/python -m twine check dist/*

  # e) fresh-venv install dry-run (clean env, only hard deps from PyPI; NO upload)
  rm -rf /tmp/ql_fresh && .venv/bin/python -m venv /tmp/ql_fresh
  /tmp/ql_fresh/bin/pip install --upgrade pip
  /tmp/ql_fresh/bin/pip install dist/quiverlab-*.whl
  /tmp/ql_fresh/bin/python -c "from quiverlab import Quiver, CC; \
    print(Quiver([1,2,3],{'a':(1,2),'b':(2,3),'c':(1,3)}).algebra(relations=['a*b'],field=CC).hochschild_cohomology(3))"

  # f) paper structure (real PDF compiles in paper.yml on GitHub)
  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/release/test_paper.py -q
  ```
  Expected: (a) OK; (b) all `tests/release` green; (c) `mkdocs build --strict` exits 0; (d) two dists + `twine check` PASSED; (e) the wheel installs into a clean venv and the README quickstart prints a Hochschild table (proves `pip install quiverlab` works with only `numpy`/`sympy`/`matplotlib`); (f) paper tests green.
- [ ] **Step 3: (no new source — assembly/acceptance)**
- [ ] **Step 4: ONE tracked background full suite, awaited**
  Start the full suite in the background, await it, and only then commit:
  ```bash
  cd /Users/marco/Desktop/HomologicalNetworks/quiverlab
  NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q   # full suite (fast+deep; qpa skips w/o GAP)
  # (run as the single tracked background job for this plan; await completion; expect all green)
  ```
  Expected: the entire suite (Plans 01–08) passes on this machine; the float gate is green; `tests/release/*` green; `tests/qpa/*` skip cleanly (no local GAP). The `qpa` job proves itself on GitHub's Linux runner.
- [ ] **Step 5: Commit**
  ```bash
  git add tests/release/test_acceptance.py
  git commit -m "test(release): Plan 08 acceptance — release-candidate surface complete and green

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R7bMM4JBnSWUHbUV1DFoMd"
  ```

---

## Release-candidate checklist (verbatim — the acceptance gate)

A quiverlab checkout is a **release candidate** when every line below is true:

- [ ] **Freshness:** `scripts/release_freshness.py` exits 0 (Plans 03–07 surfaces, docs sources, packaged `citations.references_bib_path()`, `bibliography()` all present).
- [ ] **Suite green:** `pytest -q` (fast + deep) passes on macOS/Linux locally and across the `ci.yml` matrix (OS × py3.10–3.13), on both the numba and `QUIVERLAB_NO_NUMBA=1` engine paths.
- [ ] **Float gate green:** `tests/test_no_floats.py` passes; no float/complex literal anywhere under `src/` (including `src/quiverlab/qpa/`).
- [ ] **Markers:** `fast`/`deep`/`qpa` are a **proven partition** of the collected suite (pairwise disjoint + exhaustive, `test_buckets_partition_the_suite`); `slow ⊆ deep`; no `PytestUnknownMarkWarning`.
- [ ] **QPA cross-check:** `qpa.yml` installs `passagemath-gap[qpa]`, `gap_available()` is true there, and `pytest -m qpa` passes (the `kA_2/GF(2)` fixture returns `[1, 0, 0]` from the enveloping-algebra route); the same tests **skip** cleanly with no GAP.
- [ ] **Packaging modern:** `pyproject.toml` uses the PEP 639 SPDX `license = "MIT"` (no `{text=...}`, no `License ::` classifier); `python -m build` emits no license/classifier deprecation; `[qpa]`/`[docs]`/`[fast]`/`[dev]` extras and project URLs present.
- [ ] **Docs live:** `mkdocs build --strict` exits 0 (tutorials execute, API reference auto-generates, internals chapters render); `docs.yml` deploys to `https://marcoarmenta.github.io/quiverlab/`; the URL is identical in `mkdocs.yml`, `README.md`, and `CITATION.cff`.
- [ ] **Paper compiles:** `paper.yml` builds `paper/paper.md` to PDF via `openjournals/openjournals-draft-action`; all `@cite` keys resolve in the single **packaged** `src/quiverlab/citations/references.bib` (real Plan-06 ids: `ChouhySolotar2015`, `Volkov2019`, …); body is 750–1750 words; the 2026 JOSS section set is present and ordered.
- [ ] **Build installable in a fresh venv:** `python -m build` + `twine check` PASSED; the wheel installs into a clean venv with only `numpy`/`sympy`/`matplotlib` and the README three-line quickstart runs.
- [ ] **Release wired, not fired:** `release.yml` is tag-gated (`v*`), uses OIDC trusted publishing (`id-token: write`, `pypa/gh-action-pypi-publish@release/v1`, no committed token), and guards tag == pyproject version. **No package has been uploaded** — Marco pushes the tag and configures the PyPI pending publisher.
- [ ] **Front door + community:** README badges/quickstart/links; `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `CITATION.cff` present (JOSS review checklist); `CHANGELOG.md` seeded.

---

## Frozen contract for Plan 09 (web)

Plan 09 consumes these Plan-08 outputs verbatim:

- **`QLWEB_DOCS_URL = https://marcoarmenta.github.io/quiverlab/`** — the deployed docs site (web-spec §3 docs links, §16 `/literature`). Single source: `mkdocs.yml:site_url`.
- **`src/quiverlab/citations/references.bib`** (packaged; `quiverlab.citations.references_bib_path()`) — the single canonical bibliography; Plan 09's `/literature` page renders it via `quiverlab.bibliography()` (Plan 06). Plan 08 appends only software keys (`qpa`, `gap4`, `sagemath`, `quiverlab`) to it; it commits no second `.bib` and no `docs/`-tree copy.
- **Marker vocabulary + `tests/conftest.py`** — Plan 09's `tests/webapp/` is auto-marked `fast` and rides the same `ci.yml` matrix rules.
- **`Algebra.crosscheck(...)`** — available to the web worker only where `[qpa]` is installed; the web app never requires it.

## Boundary notes

- **No real publish, ever, in this plan.** Every step goes to the tag; `release.yml` fires only on Marco's pushed tag with a PyPI pending publisher he configures. No tokens in the repo.
- **VERIFY AT EXECUTION (carried from research, flagged inline):** (1) the `passagemath-gap[qpa]` extra name (fallback `[full]`) and its exact pin; (2) the two QPA call surfaces in `qpa/scripts.py` — `AlgebraAsModuleOverEnvelopingAlgebra` (the algebra-as-`A^e`-module constructor) and the `ExtAlgebraGenerators(M, top)[1]` degreewise-dimension read (GAP 1-indexed) — pinned against the QPA manual on first `qpa.yml` run, with the `[1,0,0]` fixture as oracle; (3) exact minor versions of docs/build actions and mkdocs plugins (they drift); (4) whether the deprecated `license` table / `License ::` classifier is now a hard error vs. warning on the installed setuptools; (5) Marco's real ORCID/affiliation in `paper.md` and `CITATION.cff`; (6) re-run the PyPI name check immediately before registering. The *mathematics* (`kA_2` HH `= [1,0,0]`; the commutative square `= [1,0,0]`) is certain.
- **Trace (Plan 07) / bibliography (Plan 06)** are consumed, not built here; the freshness gate STOPs if they are absent.

---

## Self-review

**Spec coverage — §9 (packaging/CI):**
- Pure-Python wheel, hard deps `numpy`/`sympy`/`matplotlib`, extras `[fast]`/`[qpa]` → Task 3 (extras) + Task 10 (wheel/build). ✓
- `[qpa]` = passagemath-gap, macOS/Linux, loud graceful Windows message → Tasks 3–5 (`sys_platform != 'win32'` marker + runtime `_WINDOWS_MSG`). ✓
- CI: GitHub Actions macOS/Linux/Windows × py3.10–3.13; docs build; QPA cross-check job (Linux); float-free lint gate → Tasks 6, 7, 8. ✓
- Versioning semver 0.x; 1.0 at JOSS acceptance → Task 3 (`Development Status :: 4 - Beta`) + Task 12 (release doc). ✓
- **Backlog: the Plan-01-deferred `license = {text}` deprecation** → Task 3 (PEP 639 SPDX fix), with a build-warning check. ✓
- **`.github/` greenfield** → all five workflows written complete. ✓

**Spec coverage — §10 (docs):**
- mkdocs-material, auto API reference (mkdocstrings), CI-executed tutorials, non-coder-first → Task 8 (gen-files/literate-nav recipe, `mkdocs-jupyter execute: true`, `docs/index.md`). ✓
- Internals chapters included on the site ("Under the hood" nav section) → Task 8. ✓
- GitHub Pages hosting (added 2026-07-18) → Task 8 (`docs.yml`, native artifact flow, `site_url`). ✓
- Versioned docs nice-to-have → **decided: defer `mike` to 1.0**, justified (Task 8). ✓

**Spec coverage — §11 (JOSS paper):**
- `paper.md` + single-source `paper.bib` draft with real content (summary, statement of need condensing the §2 gap table, functionality, CS-first claim per Plan 04 honesty **with the quadratic/monomial scope caveat**, acknowledgements, references citing Plan 06's **real** BibTeX ids — `ChouhySolotar2015`, `Volkov2019` (not 2016), `Bardzell1997`, `NegronWitherspoon2016`, `GSZ2001`, `BGMS2005`, `Happel1989`) → Task 9 (full draft; the 2026 section set; State of the field carries the condensed gap table; Software design states the certified scope carefully). Bib is the single **packaged** `src/quiverlab/citations/references.bib` (Plan 06), never a `docs/` copy. ✓
- CI check that the paper compiles (JOSS open-journals action) → Task 9 (`paper.yml`, `openjournals/openjournals-draft-action`). ✓
- Backlog: **CC repr and other Plan-08-tagged minors** — searched `docs/plans/` and the source; the only Plan-08-tagged ledger item found is the `license={text}` deprecation (Task 3). No `Plan 08`-tagged `CC.__repr__` fix survives in the current tree (the grep for "Plan 08" across plans returned only the ROADMAP internals-hosting line and this deprecation). If execution surfaces a Plan-08-tagged minor in a committed `# TODO(Plan 08)` comment, fold it into Task 3/11 as a one-line fix. ⚠ (stated, not hidden)

**Spec §5 component 12 (`qpa`):** AR quiver / property recognizers are named in the spec as `[qpa]` capabilities but are **not** in Plan 08's scope statement (the task scopes `[qpa]` to `A.crosscheck` + the CI job); spec §12 lists native AR-quiver as a **non-goal** (available via `[qpa]` later). Plan 08 ships the `crosscheck` oracle; AR-quiver/recognizer wrappers are a post-v1 `[qpa]` extension. ⚠ (scoped deliberately — matches the roadmap row 08 wording "crosscheck oracle + CI job").

**Placeholder scan:** every `src/` symbol named (`session.py`, `scripts.py`, `crosscheck.py`, `Algebra.crosscheck`, `QpaUnavailableError`) has a complete body. The only deliberately-flagged unknowns are **GAP call surfaces** (sanctioned VERIFY-AT-EXECUTION, GAP not installable locally) and **deployment parameters** (ORCID, affiliation, exact action minors) — none are shipped `src/` placeholders. No `...`/`TODO`/`pass`-stub in any source block.

**Signature/config consistency:**
- Docs URL `https://marcoarmenta.github.io/quiverlab/` identical in `mkdocs.yml:site_url`, README, `CITATION.cff`, `[project.urls]`, and the Plan-09 `QLWEB_DOCS_URL` contract — locked by `test_community.py::test_docs_url_consistent_across_files`. ✓
- Markers `fast`/`deep`/`qpa` are a disjoint+exhaustive partition (the explicit-wins guard includes `fast`, so `fast ∩ qpa = ∅`; `slow ⇒ deep`); registered in `pyproject.toml`, assigned in `tests/conftest.py`, proven by `test_buckets_partition_the_suite`; every CI `pytest -m` expression (`-m fast`, full, `-m qpa`, `QUIVERLAB_NO_NUMBA=1`) matches a real bucket — Tasks 2, 6, 7. ✓
- `passagemath-gap[qpa]` name appears identically in `pyproject.toml` `[qpa]` and `qpa.yml`; `QUIVERLAB_REQUIRE_QPA=1` in `qpa.yml` matches `session.should_skip_qpa()` (folded into the skip predicate — no dead `setup_module` escalation). ✓
- Trusted-publishing invariants (`id-token: write`, `environment: pypi`, `@release/v1`, tag gate) consistent between `release.yml` and the PyPI pending-publisher instructions. ✓
- Repo owner `MarcoArmenta` (from `git remote`) consistent across every workflow badge, URL, and the trusted-publisher config. ✓

**Task count:** 13 (freshness gate → markers → pyproject → qpa package → qpa tests → CI → qpa CI → docs → paper → release → README → community → acceptance). Within the 10–14 band.

**Suite discipline:** every task is TDD-shaped (failing check → artifact → green), full suite green at each commit, ONE tracked background full suite at Task 13, the two-line `Co-Authored-By` + `Claude-Session` trailer inlined on every commit (repo convention, verified in `git log`), banks untouched, no real PyPI upload.

**Rework applied (round 1, all 8 items):**
1. **Bib seam** — canonical bib repointed to the **packaged** `src/quiverlab/citations/references.bib` (accessor `citations.references_bib_path()`) at all ~8 sites: freshness gate, paper.yml `cp` + `paths:` trigger, `test_paper.py` (reads via accessor; body-only key scan avoids the front-matter email), mkdocs (a `gen_bibliography.py` gen-files hook renders a `References` page from the packaged file — the chosen mechanism, no `docs/` copy), Plan-09 contract, checklist. Paper `@cite` keys rewritten to Plan 06's real ids (`ChouhySolotar2015`/`Bardzell1997`/`NegronWitherspoon2016`/`Volkov2019`/`GSZ2001`/`BGMS2005`/`Happel1989`); four software keys appended to the packaged bib (16→20) with an explicit **coverage-test floor bump** owned by Task 9. Zero `docs/references.bib` references remain (one negative assertion excepted).
2. **QPA escalation** — dead `enforce_required_in_ci`/`setup_module` removed; the `QUIVERLAB_REQUIRE_QPA` gate folded **into** `session.should_skip_qpa()` (skip predicate evaluated at collection, so REQUIRE now makes tests run-and-fail, not green-skip); pip caching added to `qpa.yml`.
3. **Markers** — `fast` added to the explicit-wins guard (kills `fast ∩ qpa`); `slow ⇒ deep`; `_DEEP_DIRS` gains `families`/`batch`; exhaustive-partition test + `slow ⊆ deep` test added; future-heavy-dir policy stated + freshness reminder.
4. **Action pins** — bumped past the Node-20 cutoff: `checkout@v7`, `setup-python@v6`, `upload-artifact@v7`, `download-artifact@v8`, `upload-pages-artifact@v5`, `deploy-pages@v5` (kept `configure-pages@v5`, `@release/v1`, `@master` w/ moving-ref caveat).
5. **Docs nav** — nav-coverage test (`test_nav_covers_internals_and_tutorials`) + `--strict omitted_files` backstop; nav comment instructs listing all existing chapters 08-11.
6. **Paper** — CS quadratic/monomial scope caveat added to Software design.
7. **scripts.py** — invented `DimensionExtNthPower` replaced by `ExtAlgebraGenerators(M, top)[1]`; module cross-check narrowed to self-Ext (the one confirmed idiom); VERIFY flag repointed.
8. **Low** — `MPLBACKEND: Agg` (top-level in `ci.yml`, job-level in `qpa.yml`/`docs.yml`); `CITATION.cff` version == pyproject release-core test; two-line commit trailer aligned to repo `git log`.

