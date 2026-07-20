# quiverlab Plan 07 — Visualization + Worked-Steps Trace Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A user draws any `kQ/I` (`A.draw()` → PNG/SVG, `A.tikz()` → paste-into-paper TikZ sharing the *same* computed coordinates), and every computation — by default (`quiverlab.verbose = True`, D9) — writes a human-readable worked-steps document (`./quiverlab_traces/<name>_<hash>.pdf`, a self-contained no-JS HTML showing TeX source when no LaTeX toolchain, plain text on demand) whose every claim is a golden-file-tested equality with the value the engine actually computed, closing with a **References** section drawn from Plan 06's citation registry.

**Architecture:** Two new pure-Python packages. `src/quiverlab/viz/` computes a layered layout in **exact integers/Fractions** (`layout()` → `LayoutData`), then renders it with matplotlib (`draw()`) or as TikZ (`tikz()`) from the identical coordinates — matplotlib becomes a **hard** dependency here (spec §9). `src/quiverlab/trace/` unifies the typed step-event taxonomy already scattered across `groebner.events` (shipped) and `resolutions_cs.trace` (Plan 04), adds the one genuinely new event `RankStep`, threads a `trace`/`verbose` recorder through the Hochschild engines, and renders the recorded events three ways (LaTeX→PDF as the primary math output, a no-JS HTML fallback showing TeX source, text) with size-eliding rules and a bibliography pulled from Plan 06's `bibliography()` registry.

**Tech Stack:** Python ≥ 3.10; **matplotlib ≥ 3.7 (new hard dep)**; the Plan-01 stack (numpy, sympy); pytest. No floats anywhere (§ the AST-gate decision below): viz layout is computed in `int`/`fractions.Fraction` and handed to matplotlib, which coerces to float *inside the library, never in our source*. LaTeX compilation via `pdflatex`/`tectonic` if on PATH (optional, graceful HTML fallback).

## Global Constraints

- **Repo root:** `/Users/marco/Desktop/HomologicalNetworks/quiverlab`. All paths below are relative to it.
- **Interpreter:** use the project venv **`/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python`** (Python 3.12). The *system* `python`/`python3` is 3.8 and MUST NOT be used — it fails on 3.10+ syntax (`X | None`, `list[int]`).
- **Thread throttle:** prefix **every** test command with `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2` (Marco's machine has crashed under agent fleets; keep thread/memory pressure low). No new parallelism is introduced in this plan.
- **Exact arithmetic only, through the `Domain` protocol.** All coefficient/matrix arithmetic that touches algebra goes through `dom.add/sub/neg/mul/inv/is_zero/eq/coerce/zero/one`. **Viz layout coordinates are computed in `int`/`fractions.Fraction`** and never enter any `Domain` operation. Matrix entries shown in a trace are *rendered from* domain elements (via `str`/a domain-aware LaTeX helper) but never combined.
- **Float ban — AST gate, and the DECISION this plan makes about it.** `tests/test_no_floats.py` scans **all of `src/quiverlab/`** (recursive `rglob("*.py")`) for float/complex literals and bare `float(...)` calls and MUST stay green. It has **no** exemption/skip mechanism today. **Decision (locked): compute all viz layout in exact `int`/`Fraction`; write NO float or complex literal and NO `float()` call anywhere under `src/quiverlab/` including `src/quiverlab/viz/`. The gate is left completely UNCHANGED — no directory exemption, no pragma — so it covers viz with the same strength as the algebra core (this is the strongest guarantee and is exactly in the spirit of D3).** Matplotlib receives `int`/`Fraction` and coerces them to float internally *inside matplotlib* (`Fraction.__float__`), which is outside `src/quiverlab/` and thus outside the gate. Two rules make this airtight and are enforced by tests in Tasks 4–5: (i) never use matplotlib's *string* form `connectionstyle="arc3,rad=0.3"` (it would force `f"...rad={float(x)}"`); use the **object** form `matplotlib.patches.ConnectionStyle.Arc3(rad=Fraction(3, 10))`; (ii) loop/arc angles are **integer degrees**, offsets/curvatures are `Fraction`, sizes are `int`. **Trade-offs considered:** the rejected alternative — exempting `src/quiverlab/viz/` in the gate plus a "no viz value reaches a Domain op" compensating test — would weaken the whole-package guarantee over a directory and rely on a hard-to-make-airtight taint test; the chosen alternative costs only mild discipline (Fractions into matplotlib, object-form connection styles) and keeps the gate absolute. A **compensating test** ships anyway (Task 4): `layout()` output contains only `int`/`Fraction`, proving no float can leak from viz into anything downstream.
- **`verbose` default TRUE (D9), but the test suite runs quiet.** `quiverlab.verbose = True` at import (Plan 09 flips it per request). `tests/conftest.py` gains an **autouse** fixture that forces `quiverlab.verbose = False` for the whole suite so unrelated tests never write `./quiverlab_traces/` files; trace tests opt back in explicitly (`quiverlab.verbose = True` inside the test, or by passing `verbose=True`/`trace=[...]`).
- **Suite pattern (Marco's 2026-07-18 directive).** Per task: run the task's **focused** test files in the **FOREGROUND** with the thread throttle (fast, low pressure). Do **NOT** run the full suite per task. After the final task (Task 13), launch the full suite **once** as a **tracked BACKGROUND** job (`run_in_background`), **await** it, and only then make the final acceptance commit and report. Every per-task commit must have its focused tests green; the single background full-suite run gates the acceptance commit. `tests/test_no_floats.py` is included in *every* task's focused set (the gate must never go red).
- **Green every commit.** No commit lands with its focused tests (incl. the float gate) failing.
- **Read-only banks.** `/Users/marco/Desktop/HomologicalNetworks/HomologicalAlgebra/` and every other bank are never read or written. This plan depends only on merged-in `src/quiverlab/` source.
- **Commits:** conventional prefixes (`feat:`, `test:`, `docs:`, `chore:`); every commit message ends with the trailer line
  `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- **Left-to-right path composition** (Assem–Simson–Skowroński): `a*b` = "first `a`, then `b`", `target(a) == source(b)`; a word is an arrow-name tuple read left to right. Layout arrows point source→target.
- **Frozen upstream interfaces consumed verbatim** (do not modify their signatures): `Domain` (`coerce/add/sub/neg/mul/inv/is_zero/eq/zero/one/characteristic/name`); `combinat.Quiver` (`.vertices` list, `.arrows` dict name→(src,tgt) insertion-ordered, `.source/.target/.word_source/.word_target/.compose_ok/.is_acyclic`); `combinat.Relation` (`.terms`, `.source`, `.target`, `.is_monomial`, `__repr__` → e.g. `a*b - c*d`); `core.Algebra` (`.domain/.dim/.T/.unit/.basis_labels/.is_unit_adapted`, **public** `.quiver`, **public** `.relations`, `.unit_adapted()`, `.multiply`); `hochschild.table.HHTable(dims, kind, algebra_repr, engine=...)` with `.dims/.kind/.top/.algebra_repr/.engine`; `hochschild.bar.coboundary_matrix/boundary_matrix/hochschild_cohomology_dims/hochschild_homology_dims`; `groebner.events.{Dispatch, ReductionStep}`; `resolutions_cs.trace.{AmbiguityEvent, ResolutionTerm, DifferentialEvent, LiftStep}` (Plan 04); `bibliography()` and `families()` (Plan 06). **This plan makes compatible extensions only:** `HHTable` **already carries** a `.references` attribute and `Algebra.hochschild_cohomology/homology` **already take** `engine` (including `"cs"`), `auto_cs`, and `trace=None` — all merged (Plans 04–06), with the engine methods setting `.references = self.citations()` (family + engine keys). Plan 07 therefore **adds only** the `verbose=None` kwarg and **extends** trace-filling — today only the CS path fills `trace` — to the bar and fast engines as well. It does **not** modify the `.references = self.citations()` contract, and does **not** drop `engine="cs"`/`auto_cs`. `Algebra` also gains `.draw()`/`.tikz()` methods and an appended (line-0-preserving) `__repr__`.

---

## The trace-event taxonomy this plan freezes (read before Task 1)

Five event *names* were requested for the unified system: `ResolutionTerm, Differential, RankStep, LiftStep, Dispatch`. Ground truth from the merged tree forces two reconciliations, both adopted here and pinned by the Task-1 freshness gate:

1. **`Differential` is actually `DifferentialEvent`** — Plan 04 (`resolutions_cs/trace.py`) named the differential event `DifferentialEvent(degree, chain, terms)` and additionally ships `AmbiguityEvent(degree, chain_words)`. Those names are already merged and are consumed by Plan 04's own tests; renaming them is out of scope and would break Plan 04. The unified taxonomy therefore uses `DifferentialEvent` (the requested `Differential`) and folds in `AmbiguityEvent`.
2. **`RankStep` does not exist yet** — no engine emits it; it is the one genuinely new event. The CS `DifferentialEvent(degree, chain, terms)` carries *symbolic bimodule terms*, which cannot represent the bar/fast engines' *numeric differential matrix + its rank over a field*. `RankStep` fills that gap and is what spec §3.8's "each differential as a matrix over the stated field (elided …) … each rank computation" maps onto for every matrix-rank engine.

**The unified taxonomy (7 event dataclasses), single import surface `quiverlab.trace.events`:**

| event | fields | origin | emitted by |
|---|---|---|---|
| `Dispatch` | `route, reason, n_relations` | `groebner.events` (shipped) | **construction** dispatch (monomial vs groebner) *and* **engine/resolution choice** (bar vs fast vs bardzell vs chouhy-solotar). Same dataclass, two contexts, distinguished by `route`'s value and its position in the event list. Its fields fit both: at engine-choice time `route` = engine name, `reason` = why, `n_relations` = `len(A.relations or ())`. This is spec §3.8's "the chosen resolution and why (dispatch decision)". |
| `ReductionStep` | `word, rule_lead, before, after` | `groebner.events` (shipped) | `groebner.reduce_comb` during construction. |
| `AmbiguityEvent` | `degree, chain_words` | `resolutions_cs.trace` (Plan 04) | CS ambiguity chains. |
| `ResolutionTerm` | `degree, n_generators, collapsed_dim` | `resolutions_cs.trace` (Plan 04) | each resolution/complex term; reused by the bar path (`n_generators = collapsed_dim = dim C^n`). |
| `DifferentialEvent` | `degree, chain, terms` | `resolutions_cs.trace` (Plan 04) | CS symbolic differential. |
| `LiftStep` | `degree, kind, detail` | `resolutions_cs.trace` (Plan 04) | CS comparison/lifting (cup/bracket) steps. |
| **`RankStep`** | `degree, side, nrows, ncols, rank, field, matrix, elided, note` | **new (Plan 07)** | every numeric differential-rank computation (bar + fast + CS matrix ranks). |

`quiverlab.trace.events` **re-exports** the five shipped/Plan-04 dataclasses from their home modules (never redefines them) and **defines** `RankStep`. If any upstream event drifts (renamed, field changed, moved), Task 1 fails loudly and execution STOPS until reconciled.

---

## Fixtures (worked out and, where noted, machine-verified)

### F1 — commutative square (layout golden)
`Q = Quiver([1,2,3,4], {"a":(1,2),"b":(2,4),"c":(1,3),"d":(3,4)})`, `A = Q.algebra(relations=["a*b - c*d"], CC)` (dim 9). Layering by longest path on the SCC-condensation (all singletons here): `depth(1)=0`, `depth(2)=depth(3)=1`, `depth(4)=2`. Columns `{1} | {2,3} | {4}`; `x = depth`; within a column of `k` vertices in `Quiver.vertices` order, `y_i = Fraction(k-1, 2) - i`. Exact positions:
```
1 -> (0, 0)      2 -> (1, 1/2)      3 -> (1, -1/2)      4 -> (2, 0)
```
Edges all `straight` (no parallels, no loops). Relation strings: `["a*b - c*d"]`.

### F2 — two loops (loop-routing golden)
`Q = Quiver([1], {"x":(1,1), "y":(1,1)})`. One vertex `1 -> (0, 0)`. Two loops become `LoopRoute`s at integer base angles `90` and `30` degrees (`angle = 90 - 60*j`, `j` = loop index at that vertex).

### F3 — parallel arrows (offset golden)
`Q = Quiver([1,2], {"a":(1,2), "b":(1,2)})`. Positions `1 -> (0,0)`, `2 -> (1,0)`. The two parallel arrows get symmetric `Fraction` bends `bend_i = Fraction(k-1-2i, 4)`: `a -> 1/4`, `b -> -1/4`.

### F4 — kA₂ (TikZ golden)
`Q = Quiver([1,2], {"a":(1,2)})`, `A = Q.algebra(field=CC)` (no relations). Positions `1->(0,0)`, `2->(1,0)`, one straight arrow. Exact TikZ source pinned in Task 6.

### F5 — HH*(k[x]/(x²), top=2) over CC (worked-steps golden — **machine-verified 2026-07-18**)
`A = truncated_polynomial(2, field=CC)` (monomial, dim 2; `repr(A).splitlines()[0]` = `"Algebra of dimension 2 over CC (computing exactly in QQ)"`). Default `engine="auto"` over a non-prime field ⇒ the **normalized bar complex** oracle. `C^n = m·(m-1)^n = 2` for all `n` (`m=2, m-1=1`). Verified coboundary matrices and ranks:

| degree n | `RankStep.side` | matrix (2×2, `str`-rendered entries) | rank |
|---|---|---|---|
| 0 | `cochain` | `[["0","0"],["0","0"]]` | 0 |
| 1 | `cochain` | `[["0","0"],["2","0"]]` | 1 |
| 2 | `cochain` | `[["0","0"],["0","0"]]` | 0 |

Trace event sequence: `Dispatch(route="normalized bar complex", …, n_relations=1)`, then per `n∈{0,1,2}` a `ResolutionTerm(degree=n, n_generators=2, collapsed_dim=2)` and a `RankStep` (above). Derived dims `HH^n = collapsed_dim(n) − rank(n) − rank(n−1)` = `[2−0−0, 2−1−0, 2−0−1] = [2, 1, 1]` — equal to `A.hochschild_cohomology(2).dims` (machine-verified `[2,1,1]`; over `GF(2)` it is `[2,2,2]`). References: provenance key `("bar",)` (a Plan 06 lowercase registry key), which `resolve_references` turns into the `Hochschild1945` entry via `bibliography()`. The golden test asserts (a) the recorded ranks `[0,1,0]`, (b) each recorded matrix equals the independently recomputed `coboundary_matrix`, (c) the renderer-derived dims equal the real `.dims`, (d) the References section lists the `bar` entry's `Hochschild1945` citation, sourced from the registry.

---

### Task 1: Executable freshness gate — pin every upstream dependency, STOP on drift

**Files:**
- Create: `tests/trace/__init__.py` (empty), `tests/trace/test_freshness_gate.py`
- Create: `tests/viz/__init__.py` (empty)

**Interfaces:**
- Consumes (all from merged Plans 03–06): `groebner.events.{Dispatch, ReductionStep}`; `resolutions_cs.trace.{AmbiguityEvent, ResolutionTerm, DifferentialEvent, LiftStep}`; `quiverlab.bibliography`; `quiverlab.families`; `HHTable`; `hochschild.bar` functions; `Algebra.hochschild_cohomology`.
- Produces: no source; a red-on-drift gate. **If this task fails at execution, STOP — a dependency this plan builds on changed shape; reconcile the plan (event fields, bibliography keys, or family catalog) before writing any Plan-07 source.**

- [ ] **Step 1: Write the failing test**

`tests/trace/__init__.py`, `tests/viz/__init__.py`: empty files.

`tests/trace/test_freshness_gate.py`:

```python
"""Plan-07 freshness gate (spec §3.7-3.9). Pins the exact shape of everything
Plan 07 builds ON TOP OF from merged Plans 03-06. Any drift here means STOP and
reconcile the plan before touching Plan-07 source. This test has NO Plan-07
imports on purpose -- it must be writable and runnable before any Plan-07 code
exists."""
from dataclasses import fields as dc_fields

import pytest

import quiverlab as ql


def _field_names(cls):
    return {f.name for f in dc_fields(cls)}


# --- shipped Groebner events (Plan 03) ---------------------------------------
def test_groebner_events_shape():
    from quiverlab.groebner.events import Dispatch, ReductionStep
    assert _field_names(Dispatch) == {"route", "reason", "n_relations"}
    assert _field_names(ReductionStep) == {"word", "rule_lead", "before", "after"}


# --- Chouhy-Solotar trace events (Plan 04) -----------------------------------
def test_cs_trace_events_shape():
    from quiverlab.resolutions_cs.trace import (
        AmbiguityEvent, ResolutionTerm, DifferentialEvent, LiftStep,
    )
    assert _field_names(AmbiguityEvent) == {"degree", "chain_words"}
    assert _field_names(ResolutionTerm) == {"degree", "n_generators", "collapsed_dim"}
    assert _field_names(DifferentialEvent) == {"degree", "chain", "terms"}
    assert _field_names(LiftStep) == {"degree", "kind", "detail"}


# --- HHTable shape (Plan 01, still frozen) -----------------------------------
def test_hhtable_shape():
    from quiverlab.hochschild.table import HHTable
    t = HHTable([1, 0], "HH^", "an algebra")
    assert t.dims == [1, 0]
    assert hasattr(t, "engine") and isinstance(t.engine, str)
    # Plan 07 will attach `.references`; it must NOT already be a slotted/blocked attr.
    t.references = ("X",)
    assert t.references == ("X",)


# --- Hochschild method contract Plan 07 EXTENDS (must not regress) ------------
def test_hochschild_method_contract_is_frozen():
    """Plan 07 only ADDS `verbose=` to hochschild_cohomology/homology; the merged
    `engine='cs'`/`auto_cs`/`trace` routing, the unknown-engine QuiverlabError, and
    the `.references = self.citations()` (family+engine) union are contracts it must
    PRESERVE (Tasks 8 & 11). Pinned here so a regression (dropping 'cs', raising
    ValueError, clobbering .references engine-only, or un-wrapping the fast list)
    fails at Task 1, not silently at Task 13."""
    import inspect
    from quiverlab.errors import QuiverlabError
    from quiverlab.core.algebra import Algebra
    params = inspect.signature(Algebra.hochschild_cohomology).parameters
    assert "auto_cs" in params and "trace" in params, params
    # unknown engine raises QuiverlabError (NOT a bare ValueError; QuiverlabError is
    # not a ValueError subclass), and its guidance still names 'cs' as a valid engine.
    A = ql.truncated_polynomial(2, field=ql.CC)
    with pytest.raises(QuiverlabError) as ei:
        A.hochschild_cohomology(1, engine="definitely-not-an-engine")
    assert "cs" in str(ei.value)
    # .references is the family+engine union (== citations()), not engine-only: use a
    # stamped catalog family (NakayamaAlgebra stamps ("nakayama","assem_book")).
    N = ql.NakayamaAlgebra(n=3, l=2)
    refs = N.hochschild_cohomology(2).references
    assert refs == N.citations() and "nakayama" in refs and "bar" in refs
    # the fast GF(p) path returns an HHTable whose `.engine` is a string (list wrapped).
    fast = ql.truncated_polynomial(2, field=ql.GF(2)).hochschild_cohomology(2)
    assert isinstance(fast.engine, str) and fast.engine


# --- bar-complex entry points Plan 07 instruments ----------------------------
def test_bar_entry_points_exist():
    from quiverlab.hochschild import bar
    for name in ("coboundary_matrix", "hochschild_cohomology_dims",
                 "boundary_matrix", "hochschild_homology_dims"):
        assert hasattr(bar, name), name


# --- Plan-06 citation registry (bibliography()) ------------------------------
def test_bibliography_registry_has_needed_keys():
    """Plan 06's bibliography() returns a Bibliography dataclass: a `.keys` TUPLE
    (lowercase REGISTRY keys) and iteration (__iter__) yielding entry views with
    .key / .formatted / .bibtex_key (and .doi/.arxiv/.topic/.annotation). There is
    NO .keys() method and NO subscripting. Plan-07's provenance map (Task 11) keys
    off these registry names -- if Plan 06 renames them, update
    quiverlab.trace.provenance.ENGINE_REFERENCES and this set together."""
    assert hasattr(ql, "bibliography"), (
        "Plan 06 must export bibliography() (the citation registry Plan 07 renders)")
    bib = ql.bibliography()
    assert hasattr(bib, "keys") and isinstance(bib.keys, tuple), (
        "bibliography().keys must be a TUPLE attribute (not a method)")
    needed = {"bar", "bardzell", "chouhy_solotar"}
    missing = needed - set(bib.keys)
    assert not missing, f"bibliography() missing registry keys Plan 07 needs: {sorted(missing)}"
    # entry views by registry key (the __iter__ protocol Plan 06 is adding this round)
    by_key = {e.key: e for e in bib}
    assert needed <= set(by_key), "bibliography() iteration must yield the needed keys"
    bar = by_key["bar"]
    for attr in ("key", "formatted", "bibtex_key"):
        assert hasattr(bar, attr), f"entry view missing .{attr}"
    # The `bar` registry key is backed by the genuine Hochschild (1945) reference the
    # golden worked-steps trace displays. STOP if Plan 06 lands without it.
    assert bar.bibtex_key == "Hochschild1945", (
        "Plan 06's `bar` entry must carry bibtex_key == 'Hochschild1945' "
        "(G. Hochschild, Ann. of Math. 46 (1945), 58-67); coordinate its addition")
    assert "Hochschild" in bar.formatted


# --- Plan-06 family catalog --------------------------------------------------
def test_family_catalog_exists():
    assert hasattr(ql, "families"), "Plan 06 must export families() (catalog discovery)"
    cat = ql.families()
    # FamilyListing is always truthy; assert the catalog actually enumerates names.
    assert cat.names(), "families() catalog is empty (FamilyListing.names() returned nothing)"


# --- verbose flag is a plain module attribute Plan 07/09 toggle --------------
def test_verbose_flag_is_settable():
    # Before Plan 07 the attribute may be absent; after Task 2 it exists and defaults True.
    # This assertion documents the contract; it is xfail-until-Task-2.
    if not hasattr(ql, "verbose"):
        pytest.xfail("quiverlab.verbose is introduced in Plan-07 Task 2")
    prev = ql.verbose
    try:
        ql.verbose = False
        assert ql.verbose is False
    finally:
        ql.verbose = prev
```

- [ ] **Step 2: Run test to verify it (mostly) fails now, and pins drift at execution**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_freshness_gate.py -q`
Expected **at execution time (Plans 03–06 merged)**: all pass except `test_verbose_flag_is_settable` which **xfails** (until Task 2). **If any of the other tests FAIL, STOP** — an upstream contract drifted; reconcile the taxonomy table / bibliography keys / family catalog references in this plan before proceeding. (Do not weaken the gate to make it pass.)

- [ ] **Step 3: (no implementation — this task is a gate)**

If drift is found, the fix is to *this plan* (and possibly `quiverlab.trace.provenance` in Task 11 / the event table above), not to weaken the assertions. Record what drifted in the commit message.

- [ ] **Step 4: Run the focused suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_freshness_gate.py tests/test_no_floats.py -q`
Expected: pass (with the one documented xfail); float gate green.

- [ ] **Step 5: Commit**

```bash
git add tests/trace/__init__.py tests/viz/__init__.py tests/trace/test_freshness_gate.py
git commit -m "test(trace): Plan-07 freshness gate pinning upstream event/bibliography/family contracts

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: Unified trace event taxonomy + `quiverlab.verbose` flag + quiet-suite conftest

**Files:**
- Create: `src/quiverlab/trace/__init__.py`, `src/quiverlab/trace/events.py`, `tests/trace/test_events.py`
- Modify: `src/quiverlab/__init__.py` (add `verbose = True` module attribute + export nothing else new yet)
- Create: `tests/conftest.py`

**Interfaces:**
- Consumes: the frozen upstream events (Task 1).
- Produces: `quiverlab.trace.events` re-exporting `Dispatch, ReductionStep` (from `groebner.events`), `AmbiguityEvent, ResolutionTerm, DifferentialEvent, LiftStep` (from `resolutions_cs.trace`), and **defining** `@dataclass RankStep(degree:int, side:str, nrows:int, ncols:int, rank:int, field:str, matrix=None, elided:bool=False, note:str="")`. `quiverlab.verbose: bool = True`.

- [ ] **Step 1: Write the failing test**

`tests/trace/test_events.py`:

```python
"""Unified trace-event taxonomy (spec §3.8; Plan-07 reconciliation)."""
from dataclasses import fields as dc_fields

import quiverlab
from quiverlab.trace.events import (
    Dispatch, ReductionStep, AmbiguityEvent, ResolutionTerm,
    DifferentialEvent, LiftStep, RankStep,
)


def _names(cls):
    return {f.name for f in dc_fields(cls)}


def test_reexports_are_the_shipped_classes_not_copies():
    import quiverlab.groebner.events as ge
    import quiverlab.resolutions_cs.trace as ct
    assert Dispatch is ge.Dispatch
    assert ReductionStep is ge.ReductionStep
    assert AmbiguityEvent is ct.AmbiguityEvent
    assert ResolutionTerm is ct.ResolutionTerm
    assert DifferentialEvent is ct.DifferentialEvent
    assert LiftStep is ct.LiftStep


def test_rankstep_is_new_and_shaped():
    assert _names(RankStep) == {
        "degree", "side", "nrows", "ncols", "rank", "field", "matrix", "elided", "note"}
    r = RankStep(degree=1, side="cochain", nrows=2, ncols=2, rank=1, field="CC",
                 matrix=[["0", "0"], ["2", "0"]])
    assert r.rank == 1 and r.elided is False and r.note == "" and r.matrix[1][0] == "2"


def test_verbose_flag_defaults_true():
    assert quiverlab.verbose is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_events.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'quiverlab.trace'`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/trace/events.py`:

```python
"""The unified worked-steps step-event taxonomy (spec §3.8).

This module is the SINGLE import surface for trace events. It RE-EXPORTS the
event dataclasses shipped by Plan 03 (groebner) and Plan 04 (Chouhy-Solotar)
from their home modules -- it never redefines them -- and DEFINES the one new
event, RankStep, which carries a numeric differential matrix over a field plus
its rank (the bar/fast engines' notion of a differential; the CS
DifferentialEvent carries symbolic bimodule terms instead).

Reconciliation notes (see the plan's taxonomy table):
  * The requested `Differential` is Plan 04's DifferentialEvent.
  * AmbiguityEvent (Plan 04) is folded in.
  * Dispatch (Plan 03) is reused for BOTH the construction route (monomial vs
    groebner) AND the engine/resolution choice (bar/fast/bardzell/chouhy-solotar);
    same dataclass, distinguished by the value of `route`.
"""
from dataclasses import dataclass

from quiverlab.groebner.events import Dispatch, ReductionStep  # noqa: F401
from quiverlab.resolutions_cs.trace import (  # noqa: F401
    AmbiguityEvent, ResolutionTerm, DifferentialEvent, LiftStep,
)


@dataclass
class RankStep:
    """One rank computation of a differential matrix over a stated field.

    `matrix` is a list[list[str]] of domain-element string renderings (kept small),
    or None when `elided` is True (matrix larger than the elision threshold -- only
    shape + rank are retained). `side` is "cochain" (d^n : C^n -> C^{n+1}) or
    "chain" (b_n : C_n -> C_{n-1}). `nrows`/`ncols` are the matrix dimensions
    (= dim of the target/source cochain space)."""
    degree: int
    side: str
    nrows: int
    ncols: int
    rank: int
    field: str
    matrix: object = None
    elided: bool = False
    note: str = ""


__all__ = [
    "Dispatch", "ReductionStep", "AmbiguityEvent", "ResolutionTerm",
    "DifferentialEvent", "LiftStep", "RankStep",
]
```

`src/quiverlab/trace/__init__.py`:

```python
"""quiverlab.trace: the worked-steps trace subsystem (spec §3.8, component 11).

Typed step events (events.py), a recording buffer with size-eliding rules
(recorder.py), engine->citation provenance (provenance.py), and three renderers
(render_text/render_latex/render_html) driven by writer.py. Default-on per D9
(quiverlab.verbose = True)."""
from quiverlab.trace.events import (  # noqa: F401
    Dispatch, ReductionStep, AmbiguityEvent, ResolutionTerm,
    DifferentialEvent, LiftStep, RankStep,
)

__all__ = [
    "Dispatch", "ReductionStep", "AmbiguityEvent", "ResolutionTerm",
    "DifferentialEvent", "LiftStep", "RankStep",
]
```

Add to `src/quiverlab/__init__.py` (after `__version__`, before the imports block so Plan 09's `getattr(ql, "verbose", …)` always resolves):

```python
# Worked-steps traces are ON by default (spec D9). Flip per-call via
# A.hochschild_cohomology(..., verbose=False) or globally via quiverlab.verbose.
verbose = True
```

and add `"verbose"` to `__all__`.

`tests/conftest.py`:

```python
"""Global test fixtures. The worked-steps trace subsystem is ON by default
(quiverlab.verbose = True, spec D9); force it OFF for the whole suite so that
unrelated computations do not write ./quiverlab_traces/ files. Trace tests that
need it opt back in explicitly (set quiverlab.verbose = True or pass
verbose=True / trace=[...])."""
import pytest

import quiverlab


@pytest.fixture(autouse=True)
def _quiet_traces():
    prev = getattr(quiverlab, "verbose", True)
    quiverlab.verbose = False
    try:
        yield
    finally:
        quiverlab.verbose = prev
```

- [ ] **Step 4: Run the focused suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_events.py tests/trace/test_freshness_gate.py tests/test_no_floats.py -q`
Expected: pass (the freshness gate's `test_verbose_flag_is_settable` no longer takes its xfail branch — `quiverlab.verbose` now exists, so it runs the assertions and passes normally).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/trace tests/trace/test_events.py tests/conftest.py src/quiverlab/__init__.py
git commit -m "feat(trace): unified event taxonomy (RankStep + re-exports) and verbose flag

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: matplotlib hard dependency + viz package skeleton + `print(A)` polish

**Files:**
- Modify: `pyproject.toml` (add `matplotlib>=3.7` to `dependencies`)
- Create: `src/quiverlab/viz/__init__.py`, `tests/viz/test_skeleton.py`
- Modify: `src/quiverlab/core/algebra.py` (extend `__repr__` to append vertices/arrows/relations, **preserving line 0**)
- Create: `tests/core/test_repr_polish.py`

**Interfaces:**
- Consumes: `matplotlib` (new hard dep), `Algebra.quiver`, `Algebra.relations`.
- Produces: `quiverlab.viz` importable with a working `matplotlib` (Agg backend forced for headless CI); `repr(A)` line 0 unchanged, with appended `vertices:` / `arrows:` / `relations:` lines when `A.quiver` is present (spec §3.7 "plain-text: vertices, arrows, relations").

- [ ] **Step 1: Write the failing tests**

`tests/viz/test_skeleton.py`:

```python
"""viz package skeleton + matplotlib hard-dependency smoke test (spec §5 c.10)."""


def test_matplotlib_is_a_hard_dependency():
    import matplotlib  # must import without any extra
    assert hasattr(matplotlib, "__version__")
    major, minor = (int(x) for x in matplotlib.__version__.split(".")[:2])
    assert (major, minor) >= (3, 7)   # the pinned floor (pyproject: matplotlib>=3.7)


def test_viz_package_imports():
    import quiverlab.viz as viz
    assert hasattr(viz, "__all__")   # the package exposes its public surface
```

`tests/core/test_repr_polish.py`:

```python
"""print(A) shows vertices, arrows, relations (spec §3.7) without changing the
HHTable title contract (repr line 0 is used verbatim as the table title)."""
from quiverlab import Quiver, CC


def test_repr_line0_is_unchanged():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(field=CC)
    assert repr(A).splitlines()[0] == f"Algebra of dimension {A.dim} over {A.domain.name}"


def test_repr_appends_quiver_summary():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    A = Q.algebra(relations=["a*b - c*d"], field=CC)
    text = repr(A)
    assert "vertices: 1, 2, 3, 4" in text
    assert "a: 1 -> 2" in text and "d: 3 -> 4" in text
    assert "relations:" in text and "a*b - c*d" in text


def test_repr_without_quiver_has_no_summary():
    # An Algebra built without a quiver (e.g. the structure-constant escape hatch)
    # has .quiver = None, so the vertices/arrows/relations block is NOT appended.
    from quiverlab.core.algebra import Algebra
    from quiverlab.fields import QQ
    one = QQ.one()
    A = Algebra(QQ, [[[one]]], [one])   # 1-dim algebra k; _quiver defaults to None
    assert A.quiver is None
    assert "vertices:" not in repr(A)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/viz/test_skeleton.py tests/core/test_repr_polish.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'quiverlab.viz'` and the appended-summary assertions fail.

- [ ] **Step 3: Write minimal implementation**

In `pyproject.toml`, change the dependencies line to:

```toml
dependencies = ["numpy>=1.21", "sympy>=1.12", "matplotlib>=3.7"]
```

**Coordination note (Plan 08 packaging).** matplotlib is a **hard runtime dependency** and must stay in `[project].dependencies` — never demoted to an optional extra. Plan 08's packaging task must keep `matplotlib>=3.7` here. The `>=3.7` floor is deliberate and the code is verified float-free against matplotlib up to 3.11 (the axis-limit `math.floor/ceil` snapping in Task 5 is what keeps 3.11's `np.isfinite` limit check happy).

Reinstall so the venv has matplotlib:

```bash
/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pip install -e '.[dev]'
```

`src/quiverlab/viz/__init__.py`:

```python
"""quiverlab.viz: draw() and tikz() for quivers with relations (spec §5 c.10,
§3.7). A layered layout is computed ONCE in exact int/Fraction coordinates
(layout.py) and rendered either with matplotlib (draw.py) or as TikZ (tikz.py)
from the identical coordinates.

FLOAT POLICY: this package writes NO float/complex literal and NO float() call.
All layout arithmetic is int/Fraction; matplotlib coerces to float internally
(outside src/), so tests/test_no_floats.py covers viz with no exemption (D3).

BACKEND POLICY: importing this package does NOT mutate matplotlib's global backend
(no matplotlib.use(...) call) -- that would hijack a user's interactive session.
draw.py builds a Figure with its own Agg canvas for file export and returns it; a
user who wants an interactive figure gets the returned Figure and shows it with
their own pyplot, untouched by us."""

__all__ = []  # draw/tikz are Algebra methods; layout() is imported directly when needed
```

In `src/quiverlab/core/algebra.py`, replace the existing `__repr__` (currently returns the dimension line + optional basis line) with a version that keeps line 0 and appends a quiver summary when present. Locate the method and replace its body with:

```python
    def __repr__(self):
        base = f"Algebra of dimension {self.dim} over {self.domain.name}"
        lines = [base]
        if self.basis_labels:
            lines.append("basis: " + ", ".join(self.basis_labels))
        q = self.quiver
        if q is not None:  # spec 3.7: plain-text shows vertices, arrows, relations
            lines.append("vertices: " + ", ".join(str(v) for v in q.vertices))
            lines.append("arrows: " + "; ".join(
                f"{n}: {s} -> {t}" for n, (s, t) in q.arrows.items()))
            rels = self.relations or []
            if rels:
                lines.append("relations: " + "; ".join(repr(r) for r in rels))
        return "\n".join(lines)
```

(Line 0 is byte-for-byte the previous first line, so `repr(A).splitlines()[0]` — the HHTable title — is unchanged. `arrows:` renders `name: s -> t`; `relations:` renders each `Relation.__repr__`, e.g. `a*b - c*d`.)

- [ ] **Step 4: Run the focused suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/viz/test_skeleton.py tests/core/test_repr_polish.py tests/hochschild -q tests/test_no_floats.py`
Expected: pass (the `tests/hochschild` sweep confirms the HHTable title still reads correctly off the unchanged line 0).

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml src/quiverlab/viz/__init__.py tests/viz/test_skeleton.py src/quiverlab/core/algebra.py tests/core/test_repr_polish.py
git commit -m "feat(viz): matplotlib hard dep, viz skeleton, print(A) shows vertices/arrows/relations

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: Layout engine — exact int/Fraction coordinates + golden coordinate tests

**Files:**
- Create: `src/quiverlab/viz/layout.py`, `tests/viz/test_layout.py`

**Interfaces:**
- Consumes: `Quiver` (`.vertices`, `.arrows`, `.source`, `.target`).
- Produces (all coordinates `int`/`fractions.Fraction`, never `float`):
  - `@dataclass(frozen=True) EdgeRoute(name, src, tgt, kind, bend)` — `kind ∈ {"straight","parallel"}`, `bend: Fraction`.
  - `@dataclass(frozen=True) LoopRoute(name, at, angle_deg:int)`.
  - `@dataclass(frozen=True) LayoutData(positions:dict, edges:tuple, loops:tuple, relations:tuple, columns:tuple)` — `positions: dict[vertex,(x,y)]` with `x:int` (depth), `y:Fraction`.
  - `layout(quiver, relations=()) -> LayoutData` (relations = iterable rendered to strings via `repr`).
  - `layer(quiver) -> dict[vertex,int]` (longest-path depth on the SCC-condensation; handles loops/cycles).

- [ ] **Step 1: Write the failing test**

`tests/viz/test_layout.py`:

```python
"""Layered layout in exact int/Fraction coordinates (spec §3.7). We golden-test
the LAYOUT DATA (coordinates + routing), never pixels."""
from fractions import Fraction
from numbers import Integral

from quiverlab import Quiver
from quiverlab.viz.layout import layout, layer, LayoutData, EdgeRoute, LoopRoute


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})


def test_layer_depths_by_longest_path():
    assert layer(_square()) == {1: 0, 2: 1, 3: 1, 4: 2}


def test_layer_handles_loops_single_column():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    assert layer(Q) == {1: 0}


def test_square_positions_are_exact():
    L = layout(_square(), relations=["a*b - c*d"])
    assert L.positions == {
        1: (0, Fraction(0)),
        2: (1, Fraction(1, 2)),
        3: (1, Fraction(-1, 2)),
        4: (2, Fraction(0)),
    }
    assert L.relations == ("a*b - c*d",)
    assert L.loops == ()
    kinds = {e.name: e.kind for e in L.edges}
    assert kinds == {"a": "straight", "b": "straight", "c": "straight", "d": "straight"}


def test_all_coordinates_are_int_or_fraction_never_float():
    """Compensating test for the AST-gate decision: no float can leak out of viz."""
    for Q, rels in [
        (_square(), ["a*b - c*d"]),
        (Quiver([1], {"x": (1, 1), "y": (1, 1)}), []),
        (Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}), []),
    ]:
        L = layout(Q, relations=rels)
        for (x, y) in L.positions.values():
            assert isinstance(x, Integral) and isinstance(y, Fraction)
            assert not isinstance(x, float) and not isinstance(y, float)
        for e in L.edges:
            assert isinstance(e.bend, Fraction) and not isinstance(e.bend, float)
        for lp in L.loops:
            assert isinstance(lp.angle_deg, Integral)


def test_two_loops_become_loop_routes_at_integer_angles():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    L = layout(Q)
    got = sorted((lp.name, lp.at, lp.angle_deg) for lp in L.loops)
    assert got == [("x", 1, 90), ("y", 1, 30)]
    assert L.edges == ()  # loops are not edges


def test_parallel_arrows_get_symmetric_fraction_bends():
    Q = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)})
    L = layout(Q)
    bends = {e.name: e.bend for e in L.edges}
    assert bends == {"a": Fraction(1, 4), "b": Fraction(-1, 4)}
    assert all(e.kind == "parallel" for e in L.edges)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/viz/test_layout.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'quiverlab.viz.layout'`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/viz/layout.py`:

```python
"""Layered layout for a quiver with relations, in EXACT int/Fraction coordinates
(spec §3.7). Vertices are placed in columns by their longest-path depth on the
strongly-connected-component condensation (so loops and cycles are handled --
all members of a cycle share a column); within a column they are centered on
integer/half-integer rows. Parallel arrows fan out with symmetric Fraction
bends; loops become LoopRoutes at integer base angles. NO floats: matplotlib
receives these exact numbers and coerces them itself (see viz/__init__.py)."""
from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class EdgeRoute:
    name: str
    src: object
    tgt: object
    kind: str          # "straight" | "parallel"
    bend: Fraction     # 0 for straight; symmetric offsets for parallel bundles


@dataclass(frozen=True)
class LoopRoute:
    name: str
    at: object
    angle_deg: int     # integer base angle of the self-arc


@dataclass(frozen=True)
class LayoutData:
    positions: dict     # vertex -> (x:int, y:Fraction)
    edges: tuple        # tuple[EdgeRoute, ...]  (non-loop arrows)
    loops: tuple        # tuple[LoopRoute, ...]
    relations: tuple    # tuple[str, ...]
    columns: tuple      # tuple[tuple[vertex, ...], ...] by ascending depth


def _sccs(quiver):
    """Strongly-connected components via iterative Tarjan; returns
    (comp_id: dict[vertex,int], order: list of components in reverse finish order)."""
    index = {}
    low = {}
    onstack = {}
    stack = []
    comp = {}
    counter = [0]
    ncomp = [0]
    adj = {v: [] for v in quiver.vertices}
    for _n, (s, t) in quiver.arrows.items():
        adj[s].append(t)
    for root in quiver.vertices:
        if root in index:
            continue
        work = [(root, 0)]                 # (vertex, next-neighbor-index)
        while work:
            v, i = work[-1]
            if i == 0:
                index[v] = low[v] = counter[0]
                counter[0] += 1
                stack.append(v)
                onstack[v] = True
            recursed = False
            neigh = adj[v]
            while i < len(neigh):
                w = neigh[i]
                i += 1
                if w not in index:
                    work[-1] = (v, i)
                    work.append((w, 0))
                    recursed = True
                    break
                if onstack.get(w):
                    low[v] = min(low[v], index[w])
            if recursed:
                continue
            work[-1] = (v, i)
            if low[v] == index[v]:
                while True:
                    w = stack.pop()
                    onstack[w] = False
                    comp[w] = ncomp[0]
                    if w == v:
                        break
                ncomp[0] += 1
            work.pop()
            if work:
                p, _pi = work[-1]
                low[p] = min(low[p], low[v])
    return comp


def layer(quiver):
    """Longest-path depth per vertex on the SCC condensation (handles cycles/loops)."""
    comp = _sccs(quiver)
    ncomp = max(comp.values()) + 1 if comp else 0
    cadj = {c: set() for c in range(ncomp)}
    cindeg = {c: 0 for c in range(ncomp)}
    for _name, (s, t) in quiver.arrows.items():
        cs, ct = comp[s], comp[t]
        if cs != ct and ct not in cadj[cs]:
            cadj[cs].add(ct)
            cindeg[ct] += 1
    depth = {c: 0 for c in range(ncomp)}
    queue = [c for c in range(ncomp) if cindeg[c] == 0]
    order = []
    while queue:
        c = queue.pop()
        order.append(c)
        for d in cadj[c]:
            if depth[c] + 1 > depth[d]:
                depth[d] = depth[c] + 1
            cindeg[d] -= 1
            if cindeg[d] == 0:
                queue.append(d)
    return {v: depth[comp[v]] for v in quiver.vertices}


def layout(quiver, relations=()):
    depth = layer(quiver)
    maxd = max(depth.values(), default=0)
    columns = []
    positions = {}
    for d in range(maxd + 1):
        col = [v for v in quiver.vertices if depth[v] == d]  # Quiver.vertices order
        columns.append(tuple(col))
        k = len(col)
        for i, v in enumerate(col):
            positions[v] = (d, Fraction(k - 1, 2) - i)
    # bundle arrows by (src, tgt); loops split off; parallels get symmetric bends
    bundles = {}
    loops = []
    for name, (s, t) in quiver.arrows.items():
        if s == t:
            loops.append(name)
        else:
            bundles.setdefault((s, t), []).append(name)
    loop_routes = []
    loop_counter = {}
    for name in loops:
        s = quiver.source(name)
        j = loop_counter.get(s, 0)
        loop_counter[s] = j + 1
        loop_routes.append(LoopRoute(name=name, at=s, angle_deg=90 - 60 * j))
    edges = []
    for (s, t), names in bundles.items():
        k = len(names)
        for i, name in enumerate(names):
            if k == 1:
                edges.append(EdgeRoute(name=name, src=s, tgt=t, kind="straight",
                                       bend=Fraction(0)))
            else:
                edges.append(EdgeRoute(name=name, src=s, tgt=t, kind="parallel",
                                       bend=Fraction(k - 1 - 2 * i, 4)))
    return LayoutData(
        positions=positions,
        edges=tuple(edges),
        loops=tuple(loop_routes),
        relations=tuple(repr(r) if not isinstance(r, str) else r for r in relations),
        columns=tuple(columns),
    )
```

- [ ] **Step 4: Run the focused suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/viz/test_layout.py tests/test_no_floats.py -q`
Expected: pass (exact coordinates for F1; loop angles 90/30 for F2; bends ±1/4 for F3; **the compensating "no float leaks" test passes**; the float gate stays green because `layout.py` uses only `int`/`Fraction`).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/viz/layout.py tests/viz/test_layout.py
git commit -m "feat(viz): exact int/Fraction layered layout (SCC-condensation depths)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: `A.draw()` — matplotlib rendering + PNG/SVG export

**Files:**
- Create: `src/quiverlab/viz/draw.py`, `tests/viz/test_draw.py`
- Modify: `src/quiverlab/core/algebra.py` (add `draw()` method)

**Interfaces:**
- Consumes: `LayoutData`, `layout()`, matplotlib, `Algebra.quiver`, `Algebra.relations`.
- Produces:
  - `draw_quiver(quiver, relations=(), file=None) -> matplotlib.figure.Figure` — vertices as circles at `positions`, straight/parallel arrows as `FancyArrowPatch` (parallel bends via `ConnectionStyle.Arc3(rad=<Fraction>)`), loops as `patches.Arc` (integer angles), the relation list rendered as a text block below the diagram; `file=` (str/Path) saves PNG or SVG by extension and returns the Figure.
  - `Algebra.draw(self, file=None)` delegating to `draw_quiver(self.quiver, self.relations, file=file)`.
- Test the **structure** of the returned Figure (artist counts, that positions match `layout()`), never pixels.

- [ ] **Step 1: Write the failing test**

`tests/viz/test_draw.py`:

```python
"""A.draw(): matplotlib rendering from the exact layout. Structure-only asserts.
No matplotlib.use() here -- draw_quiver builds its own Agg canvas, so no global
backend is needed or mutated."""
from fractions import Fraction

from quiverlab import Quiver, CC
from quiverlab.viz.draw import draw_quiver
from quiverlab.viz.layout import layout


def _square():
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})


def test_draw_returns_figure_with_one_axes():
    fig = draw_quiver(_square(), relations=["a*b - c*d"])
    assert fig.__class__.__name__ == "Figure"
    assert len(fig.axes) == 1


def test_draw_places_one_label_per_vertex_at_layout_coords():
    Q = _square()
    fig = draw_quiver(Q, relations=["a*b - c*d"])
    ax = fig.axes[0]
    L = layout(Q, relations=["a*b - c*d"])
    texts = {t.get_text(): t.get_position() for t in ax.texts}
    for v, (x, y) in L.positions.items():
        assert str(v) in texts
        px, py = texts[str(v)]
        # matplotlib stores float positions; compare to the exact layout coerced.
        assert (px, py) == (float(x), float(y))


def test_draw_renders_relation_list_text():
    fig = draw_quiver(_square(), relations=["a*b - c*d"])
    ax = fig.axes[0]
    assert any("a*b - c*d" in t.get_text() for t in ax.texts)


def test_draw_one_patch_per_arrow_plus_loops():
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    fig = draw_quiver(Q)
    ax = fig.axes[0]
    # two loops -> two Arc patches; no straight/parallel edges here.
    arcs = [p for p in ax.patches if p.__class__.__name__ == "Arc"]
    assert len(arcs) == 2


def test_draw_file_export_png_and_svg(tmp_path):
    Q = _square()
    png = tmp_path / "A.png"
    svg = tmp_path / "A.svg"
    draw_quiver(Q, relations=["a*b - c*d"], file=str(png))
    draw_quiver(Q, relations=["a*b - c*d"], file=str(svg))
    assert png.exists() and png.stat().st_size > 0
    assert svg.exists() and svg.read_text().lstrip().startswith("<?xml")


def test_algebra_draw_method():
    A = _square().algebra(relations=["a*b - c*d"], field=CC)
    fig = A.draw()
    assert fig.__class__.__name__ == "Figure"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/viz/test_draw.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'quiverlab.viz.draw'`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/viz/draw.py`:

```python
"""A.draw(): render the exact layered layout with matplotlib (spec §3.7, §5 c.10).

FLOAT POLICY: no float/complex literal and no float() call in this file. Vertex
radius, curvature and offsets are int/Fraction; loop angles are integer degrees;
parallel bends use ConnectionStyle.Arc3(rad=<Fraction>) (the OBJECT form -- never
the "arc3,rad=0.3" string form, which would require formatting a float). Matplotlib
coerces the Fractions to float internally, outside src/quiverlab/.

Axis LIMITS are the one place matplotlib rejects Fraction bounds: matplotlib >= 3.11
runs np.isfinite on the limits, which raises TypeError on a Fraction (even an
integer-valued one). So the limits are snapped to ints with math.floor/math.ceil --
still float()-free (floor/ceil return int). Artist geometry (circles, arcs, arrows,
text) takes Fractions fine; only set_xlim/set_ylim need ints.

BACKEND: a Figure with its own FigureCanvasAgg is built directly -- no pyplot, no
matplotlib.use() -- so drawing never mutates the user's global backend. The returned
Figure renders inline (its Agg canvas) and fig.savefig handles PNG (Agg) and SVG
(matplotlib swaps to an SVG canvas by file extension)."""
import math
from fractions import Fraction

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.patches import Arc, Circle, ConnectionStyle, FancyArrowPatch

from quiverlab.viz.layout import layout

_VR = Fraction(1, 6)     # vertex circle radius
_LOOP_W = Fraction(1, 2)  # loop arc width/height


def draw_quiver(quiver, relations=(), file=None):
    L = layout(quiver, relations=relations)
    fig = Figure(figsize=(7, 5))
    FigureCanvasAgg(fig)                 # attach an Agg canvas WITHOUT touching pyplot
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect("equal")
    ax.axis("off")

    for v, (x, y) in L.positions.items():
        ax.add_patch(Circle((x, y), _VR, fill=False, linewidth=1, zorder=2))
        ax.text(x, y, str(v), ha="center", va="center", zorder=3)

    for e in L.edges:
        sx, sy = L.positions[e.src]
        tx, ty = L.positions[e.tgt]
        style = ConnectionStyle.Arc3(rad=e.bend)  # bend is a Fraction (0 for straight)
        arrow = FancyArrowPatch(
            (sx, sy), (tx, ty), connectionstyle=style, arrowstyle="-|>",
            mutation_scale=15, shrinkA=12, shrinkB=12, linewidth=1, zorder=1)
        ax.add_patch(arrow)
        mx, my = (sx + tx) * Fraction(1, 2), (sy + ty) * Fraction(1, 2) + e.bend
        ax.text(mx, my, e.name, ha="center", va="center", zorder=3)

    for lp in L.loops:
        cx, cy = L.positions[lp.at]
        # place the self-arc centered a little off the vertex along the base angle
        ax.add_patch(Arc((cx, cy + _VR), _LOOP_W, _LOOP_W,
                         angle=lp.angle_deg, theta1=200, theta2=520, linewidth=1, zorder=1))
        ax.text(cx, cy + _VR + _LOOP_W, lp.name, ha="center", va="bottom", zorder=3)

    xs = [x for (x, _y) in L.positions.values()]
    ys = [y for (_x, y) in L.positions.values()]
    pad = 1
    # Snap axis limits to ints: matplotlib>=3.11 rejects Fraction bounds (np.isfinite).
    x_lo, x_hi = math.floor(min(xs)) - pad, math.ceil(max(xs)) + pad
    y_lo, y_hi = math.floor(min(ys)) - pad, math.ceil(max(ys)) + pad + 1
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)

    if L.relations:
        ax.text(x_lo, y_lo, "relations:  " + ";  ".join(L.relations),
                ha="left", va="bottom", zorder=3)

    if file is not None:
        fig.savefig(str(file), bbox_inches="tight")
    return fig
```

In `src/quiverlab/core/algebra.py`, add the method to `class Algebra` (near the other user-facing methods):

```python
    def draw(self, file=None):
        """Draw the quiver (matplotlib): vertices by depth, loops as self-arcs,
        parallel arrows fanned out, the relation list below (spec §3.7). Returns
        the Figure; pass file="A.png"/"A.svg" to also save it."""
        if self.quiver is None:
            from quiverlab.errors import QuiverlabError
            raise QuiverlabError(
                "this Algebra has no quiver to draw",
                hint="build it via Quiver.algebra(...) rather than from_structure_constants")
        from quiverlab.viz.draw import draw_quiver
        return draw_quiver(self.quiver, self.relations or [], file=file)
```

- [ ] **Step 4: Run the focused suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/viz/test_draw.py tests/test_no_floats.py -q`
Expected: pass (Figure structure matches the layout; PNG + SVG export; two Arc patches for the two-loop quiver; float gate green — no float literal/`float()` in `draw.py`; the test's own `float(x)` call lives in the *test*, which the gate does not scan).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/viz/draw.py tests/viz/test_draw.py src/quiverlab/core/algebra.py
git commit -m "feat(viz): A.draw() matplotlib rendering with PNG/SVG export

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 6: `A.tikz()` — TikZ from the same layout + kA₂ golden file

**Files:**
- Create: `src/quiverlab/viz/tikz.py`, `tests/viz/test_tikz.py`, `tests/viz/golden/ka2.tikz`
- Modify: `src/quiverlab/core/algebra.py` (add `tikz()` method)

**Interfaces:**
- Consumes: `layout()`, `LayoutData`, `Algebra.quiver`, `Algebra.relations`.
- Produces:
  - `tikz_quiver(quiver, relations=()) -> str` — a `tikzpicture` environment placing `\node` per vertex at the exact layout coordinate (integer → `k`, `Fraction p/q` → pgfmath `{p/q}`), a `\draw[->]` per arrow (straight `--`; parallel `to[bend left=…]`; loops `to[loop …]`), and a trailing relations node when relations exist. Shares `layout()` with `draw()` (identical coordinates — asserted).
  - `Algebra.tikz(self) -> str`.

- [ ] **Step 1: Write the failing test + golden file**

`tests/viz/golden/ka2.tikz` (exact expected output for F4 — `Quiver([1,2],{"a":(1,2)}).algebra()`):

```
\begin{tikzpicture}[>=stealth]
  \node[draw, circle] (v1) at (0, 0) {$1$};
  \node[draw, circle] (v2) at (1, 0) {$2$};
  \draw[->] (v1) -- (v2) node[midway, above] {$a$};
\end{tikzpicture}
```

`tests/viz/test_tikz.py`:

```python
"""A.tikz(): TikZ from the SAME layout coordinates as draw() (spec §3.7)."""
import pathlib
from fractions import Fraction

from quiverlab import Quiver, CC
from quiverlab.viz.tikz import tikz_quiver
from quiverlab.viz.layout import layout

GOLDEN = pathlib.Path(__file__).parent / "golden" / "ka2.tikz"


def test_ka2_matches_golden_exactly():
    Q = Quiver([1, 2], {"a": (1, 2)})
    assert tikz_quiver(Q, relations=[]) == GOLDEN.read_text()


def test_tikz_uses_the_same_layout_as_draw():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    L = layout(Q, relations=["a*b - c*d"])
    src = tikz_quiver(Q, relations=["a*b - c*d"])
    # every vertex node uses its exact layout coordinate (half-integers via pgf {p/q})
    assert "(v1) at (0, 0)" in src
    assert "(v2) at (1, {1/2})" in src
    assert "(v3) at (1, {-1/2})" in src
    assert "(v4) at (2, 0)" in src


def test_tikz_emits_relations_node_when_present():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    src = tikz_quiver(Q, relations=["a*b - c*d"])
    assert "relations: $a*b - c*d$" in src


def test_algebra_tikz_method():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(field=CC)
    assert A.tikz() == (pathlib.Path(__file__).parent / "golden" / "ka2.tikz").read_text()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/viz/test_tikz.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'quiverlab.viz.tikz'`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/viz/tikz.py`:

```python
"""A.tikz(): the SAME layered layout as draw(), emitted as TikZ (spec §3.7).
Coordinates are exact: an integer prints as itself; a Fraction p/q prints as the
pgfmath expression {p/q}, which pgf evaluates -- so the emitted SOURCE contains
no float literal (this file is float-free like all of viz)."""
from fractions import Fraction

from quiverlab.viz.layout import layout


def _coord(z):
    z = Fraction(z)
    if z.denominator == 1:
        return str(z.numerator)
    return "{%d/%d}" % (z.numerator, z.denominator)


def tikz_quiver(quiver, relations=()):
    L = layout(quiver, relations=relations)
    lines = [r"\begin{tikzpicture}[>=stealth]"]
    for v in quiver.vertices:
        x, y = L.positions[v]
        lines.append(r"  \node[draw, circle] (v%s) at (%s, %s) {$%s$};"
                     % (v, _coord(x), _coord(y), v))
    for e in L.edges:
        if e.kind == "straight":
            lines.append(r"  \draw[->] (v%s) -- (v%s) node[midway, above] {$%s$};"
                         % (e.src, e.tgt, e.name))
        else:  # parallel: bend proportionally to the Fraction offset (integer degrees)
            deg = int(e.bend * 60)
            side = "left" if deg >= 0 else "right"
            lines.append(r"  \draw[->] (v%s) to[bend %s=%d] node[midway, above] {$%s$} (v%s);"
                         % (e.src, side, abs(deg), e.name, e.tgt))
    for lp in L.loops:
        lines.append(r"  \draw[->] (v%s) to[loop, in=%d, out=%d] node {$%s$} (v%s);"
                     % (lp.at, lp.angle_deg - 20, lp.angle_deg + 20, lp.name, lp.at))
    if L.relations:
        lines.append(r"  \node[align=left, below] at (current bounding box.south) "
                     r"{relations: %s};" % ";  ".join("$%s$" % r for r in L.relations))
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines) + "\n"
```

In `src/quiverlab/core/algebra.py`, add:

```python
    def tikz(self):
        """TikZ source for the quiver, sharing draw()'s exact layout coordinates
        (spec §3.7). Paste into a LaTeX document."""
        if self.quiver is None:
            from quiverlab.errors import QuiverlabError
            raise QuiverlabError(
                "this Algebra has no quiver to render as TikZ",
                hint="build it via Quiver.algebra(...) rather than from_structure_constants")
        from quiverlab.viz.tikz import tikz_quiver
        return tikz_quiver(self.quiver, self.relations or [])
```

- [ ] **Step 4: Run the focused suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/viz/test_tikz.py tests/test_no_floats.py -q`
Expected: pass (kA₂ matches the golden byte-for-byte; the square's half-integer coords print as `{1/2}`/`{-1/2}`; draw() and tikz() share `layout()`).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/viz/tikz.py tests/viz/test_tikz.py tests/viz/golden/ka2.tikz src/quiverlab/core/algebra.py
git commit -m "feat(viz): A.tikz() from shared layout + kA_2 golden file

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 7: Trace recorder — buffering, size-eliding rules, verbose resolution

**Files:**
- Create: `src/quiverlab/trace/recorder.py`, `tests/trace/test_recorder.py`

**Interfaces:**
- Consumes: `trace.events.RankStep`, `Domain`.
- Produces:
  - Constants `MAX_EVENTS = 5000`, `MATRIX_ELISION_CELLS = 400`, `TERMS_ELISION = 100`.
  - `class Trace` — list-compatible (`.append` alias for `.record`, `__iter__`, `__len__`, `.events`), caps at `MAX_EVENTS` (further events dropped, `.elided_events` counts them), `.elision_notes` list. The `.append` alias makes it drop-in where the shipped `groebner.reduce_comb(trace=...)` expects a plain list.
  - `rankstep(degree, side, D, ncols, nrows, rank, dom) -> RankStep` — renders the matrix `D` (list-of-lists of domain elements) to `list[list[str]]` via `str`, eliding to `matrix=None, elided=True, note=…` when `nrows*ncols > MATRIX_ELISION_CELLS`.
  - `resolve_verbose(per_call, global_flag) -> bool` (per-call overrides global).

- [ ] **Step 1: Write the failing test**

`tests/trace/test_recorder.py`:

```python
"""Trace recorder: buffering caps, matrix elision, verbose resolution (spec §3.8
'performance guard: verbose must not blow up long computations')."""
from quiverlab.fields import CC
from quiverlab.trace.recorder import (
    Trace, rankstep, resolve_verbose, MAX_EVENTS, MATRIX_ELISION_CELLS,
)
from quiverlab.trace.events import RankStep, Dispatch


def test_trace_is_list_compatible_for_groebner():
    tr = Trace()
    tr.append(Dispatch(route="monomial", reason="x", n_relations=1))  # .append == .record
    assert len(tr) == 1 and list(tr)[0].route == "monomial"


def test_buffer_cap_drops_and_counts_overflow():
    tr = Trace(max_events=3)
    for i in range(10):
        tr.append(Dispatch(route=str(i), reason="", n_relations=0))
    assert len(tr) == 3
    assert tr.elided_events == 7
    assert any("elided" in n.lower() for n in tr.elision_notes)


def test_rankstep_keeps_small_matrix():
    dom = CC.make_domain([CC.parse_entry(0), CC.parse_entry(1), CC.parse_entry(2)])
    D = [[dom.zero(), dom.zero()], [dom.coerce(2), dom.zero()]]
    r = rankstep(1, "cochain", D, 2, 2, 1, dom)
    assert isinstance(r, RankStep)
    assert r.elided is False
    assert r.matrix == [["0", "0"], ["2", "0"]]
    assert r.field == dom.name


def test_rankstep_elides_large_matrix():
    dom = CC.make_domain([CC.parse_entry(0), CC.parse_entry(1)])
    n = 30  # 30*30 = 900 > MATRIX_ELISION_CELLS (400)
    D = [[dom.zero()] * n for _ in range(n)]
    r = rankstep(4, "cochain", D, n, n, 0, dom)
    assert r.elided is True and r.matrix is None
    assert str(MATRIX_ELISION_CELLS) in r.note or "elided" in r.note.lower()
    assert r.nrows == n and r.ncols == n and r.rank == 0  # shape + rank always kept


def test_resolve_verbose_per_call_overrides_global():
    assert resolve_verbose(None, True) is True
    assert resolve_verbose(None, False) is False
    assert resolve_verbose(False, True) is False
    assert resolve_verbose(True, False) is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_recorder.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'quiverlab.trace.recorder'`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/trace/recorder.py`:

```python
"""Trace recorder: a bounded, list-compatible event buffer with size-eliding rules
(spec §3.8 performance guard). verbose must NOT blow up a long computation:

  * the buffer is capped at MAX_EVENTS (5000); beyond that, events are dropped and
    counted (one elision note), so memory is O(MAX_EVENTS) regardless of depth;
  * a differential matrix larger than MATRIX_ELISION_CELLS (400 = 20x20) is not
    stored -- only its (shape, rank) survive, with an explicit elision note (so a
    deep resolution records at most one small ResolutionTerm + one RankStep per
    degree, i.e. O(top) small records, never a giant matrix dump);
  * DifferentialEvent term lists over TERMS_ELISION (100) are the renderers'
    responsibility to truncate (Task 9/10) -- the recorder keeps the cap constant
    available to them.

Concretely: a top=40 monomial resolution records ~41 ResolutionTerm + ~41 RankStep
(each RankStep either a <=400-cell matrix or a one-line elision note) plus the
capped construction ReductionSteps -- well under MAX_EVENTS, bounded memory."""
from quiverlab.trace.events import RankStep

MAX_EVENTS = 5000
MATRIX_ELISION_CELLS = 400
TERMS_ELISION = 100


class Trace:
    """A list-like event sink. `.append` is an alias for `.record` so this drops
    in wherever a plain list `trace` is expected (e.g. groebner.reduce_comb)."""

    def __init__(self, max_events=MAX_EVENTS):
        self.events = []
        self.max_events = max_events
        self.elided_events = 0
        self.elision_notes = []

    def record(self, event):
        if len(self.events) >= self.max_events:
            if self.elided_events == 0:
                self.elision_notes.append(
                    f"event buffer full at {self.max_events}; further steps elided")
            self.elided_events += 1
            return
        self.events.append(event)

    append = record  # list-compatible

    def __iter__(self):
        return iter(self.events)

    def __len__(self):
        return len(self.events)


def rankstep(degree, side, D, ncols, nrows, rank, dom):
    """Build a RankStep from a domain-element matrix D (list of lists). Elide the
    matrix body (keep shape + rank) when it exceeds MATRIX_ELISION_CELLS."""
    cells = nrows * ncols
    if cells > MATRIX_ELISION_CELLS:
        return RankStep(
            degree=degree, side=side, nrows=nrows, ncols=ncols, rank=rank,
            field=dom.name, matrix=None, elided=True,
            note=(f"{nrows}x{ncols} matrix over {dom.name} elided "
                  f"(> {MATRIX_ELISION_CELLS} cells); rank {rank} recorded"))
    rendered = [[str(D[i][j]) for j in range(ncols)] for i in range(nrows)]
    return RankStep(
        degree=degree, side=side, nrows=nrows, ncols=ncols, rank=rank,
        field=dom.name, matrix=rendered, elided=False, note="")


def resolve_verbose(per_call, global_flag):
    """Per-call verbose (True/False) overrides the global flag; None defers to it."""
    return bool(global_flag) if per_call is None else bool(per_call)
```

- [ ] **Step 4: Run the focused suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_recorder.py tests/test_no_floats.py -q`
Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/trace/recorder.py tests/trace/test_recorder.py
git commit -m "feat(trace): bounded recorder with matrix elision + verbose resolution

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 8: Emit events from the Hochschild engines + thread `verbose`/`trace`

**Files:**
- Modify: `src/quiverlab/hochschild/bar.py` (emit `ResolutionTerm` + `RankStep` per degree; add `trace=None`)
- Modify: `src/quiverlab/core/algebra.py` (`hochschild_cohomology`/`hochschild_homology` gain **only** `verbose=None` — `engine` incl. `"cs"`, `auto_cs`, and `trace` are ALREADY merged and are KEPT; emit the engine-choice `Dispatch` on the bar/fast paths; **preserve** the `_route_to_cs` CS branch, the `QuiverlabError` unknown-engine guard, and `.references = self.citations()`; **no file writing yet** — Task 11 wires `write_trace`)
- Create: `tests/trace/test_emit.py`

**Interfaces:**
- Consumes: `Trace`, `rankstep`, `resolve_verbose`, `ResolutionTerm`, `Dispatch`, `quiverlab.verbose`.
- Produces:
  - `hochschild_cohomology_dims(A, top, max_cells=..., trace=None)` / `hochschild_homology_dims(A, top, max_cells=..., trace=None)` — when `trace is not None`, record one `ResolutionTerm(degree=n, n_generators=cn, collapsed_dim=cn)` and one `RankStep` per degree (side `"cochain"`/`"chain"`).
  - `Algebra.hochschild_cohomology(self, top, max_cells=..., engine="auto", auto_cs=False, verbose=None, trace=None)` (and `_homology` mirror): KEEP the merged `("auto","bar","fast","cs")` `QuiverlabError` guard, the `_route_to_cs`/`_use_fast_engine` routing, and `.references = self.citations()`; ADD `verbose` resolution, build/receive a recorder, record the engine-choice `Dispatch` on the bar/fast paths, thread `trace=rec` into the bar/fast/CS engines, return the `HHTable`. (Writing the file is Task 11; `.references` stays `self.citations()`.)

- [ ] **Step 1: Write the failing test**

`tests/trace/test_emit.py`:

```python
"""Engines emit typed events; claims equal computed values (spec §3.8, §8)."""
import quiverlab
from quiverlab import truncated_polynomial, CC, GF
from quiverlab.trace.events import Dispatch, ResolutionTerm, RankStep
from quiverlab.trace.recorder import Trace
from quiverlab.hochschild.bar import coboundary_matrix, hochschild_cohomology_dims
from quiverlab.fields.linalg import rank


def test_bar_dims_emit_terms_and_ranks():
    A = truncated_polynomial(2, field=CC)
    tr = Trace()
    table = hochschild_cohomology_dims(A, 2, trace=tr)
    terms = [e for e in tr if isinstance(e, ResolutionTerm)]
    ranks = [e for e in tr if isinstance(e, RankStep)]
    assert [t.degree for t in terms] == [0, 1, 2]
    assert [t.collapsed_dim for t in terms] == [2, 2, 2]
    assert [r.rank for r in ranks] == [0, 1, 0]           # machine-verified (F5)
    assert ranks[1].matrix == [["0", "0"], ["2", "0"]]     # d^1 over CC
    # claims == computed: recorded rank equals an independent recomputation
    B = A.unit_adapted()
    D, nc, nr = coboundary_matrix(B, 1, 4_000_000)
    assert ranks[1].rank == rank(D, B.domain)
    assert table.dims == [2, 1, 1]


def test_algebra_records_engine_dispatch_via_trace_param():
    A = truncated_polynomial(2, field=CC)
    tr = []
    A.hochschild_cohomology(2, trace=tr)
    disp = [e for e in tr if isinstance(e, Dispatch)]
    assert disp and disp[0].route == "normalized bar complex"
    assert disp[0].n_relations == 1                        # x^2


def test_prime_field_uses_fast_engine_dispatch():
    A = truncated_polynomial(2, field=GF(2))
    tr = []
    A.hochschild_cohomology(2, trace=tr)
    routes = [e.route for e in tr if isinstance(e, Dispatch)]
    assert routes and "fast" in routes[0].lower()


def test_explicit_trace_does_not_flip_global_verbose(tmp_path, monkeypatch):
    # passing trace=[] is programmatic: it must NOT write a file even if verbose is on
    monkeypatch.chdir(tmp_path)
    quiverlab.verbose = True
    try:
        A = truncated_polynomial(2, field=CC)
        A.hochschild_cohomology(2, trace=[])
    finally:
        quiverlab.verbose = False
    assert not (tmp_path / "quiverlab_traces").exists()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_emit.py -q`
Expected: FAIL — `hochschild_cohomology_dims()` takes no `trace` kwarg / no events recorded.

- [ ] **Step 3: Write minimal implementation**

In `src/quiverlab/hochschild/bar.py`, add the imports at the top:

```python
from quiverlab.trace.events import ResolutionTerm
from quiverlab.trace.recorder import rankstep
```

Modify `hochschild_cohomology_dims` to thread `trace` and emit per degree (keep the existing dim arithmetic; add only the emission):

```python
def hochschild_cohomology_dims(A, top, max_cells=4_000_000, trace=None):
    B = A.unit_adapted()
    dom = B.domain
    m = B.dim
    ranks = []
    dims = []
    prev = 0
    for n in range(top + 1):
        D, ncols, nrows = coboundary_matrix(B, n, max_cells)
        r = rank(D, dom) if nrows and ncols else 0
        cn = m * (m - 1) ** n
        dims.append(cn - r - prev)
        if trace is not None:
            trace.append(ResolutionTerm(degree=n, n_generators=cn, collapsed_dim=cn))
            trace.append(rankstep(n, "cochain", D, ncols, nrows, r, dom))
        ranks.append(r)
        prev = r
    return HHTable(dims, "HH^", repr(A).splitlines()[0])
```

Mirror in `hochschild_homology_dims` (side `"chain"`, `boundary_matrix`, its existing `ranks=[0]`/dim arithmetic preserved; emit `ResolutionTerm(degree=n, n_generators=cn, collapsed_dim=cn)` and `trace.append(rankstep(n, "chain", D, ncols, nrows, r, dom))` inside the loop).

In `src/quiverlab/core/algebra.py`, **start from the CURRENT merged bodies** of `hochschild_cohomology`/`hochschild_homology` (which already route `engine="cs"`/`auto_cs` to Chouhy–Solotar, wrap the fast `engine_cohomology_dims` `list[int]` in an `HHTable`, and set `.references = self.citations()`) and **add only** `verbose=None` plus engine-side trace recording. Do **not** drop the `_route_to_cs` CS branch, the `QuiverlabError` unknown-engine guard, the `auto_cs` kwarg, or the `self.citations()` references — all merged. Cohomology becomes:

```python
    def hochschild_cohomology(self, top, max_cells=4_000_000, engine="auto",
                              auto_cs=False, verbose=None, trace=None):
        import quiverlab
        from quiverlab.hochschild.bar import hochschild_cohomology_dims
        from quiverlab.hochschild.table import HHTable
        from quiverlab.trace.events import Dispatch
        from quiverlab.trace.recorder import Trace, resolve_verbose
        if engine not in ("auto", "bar", "fast", "cs"):
            raise QuiverlabError(f"unknown engine {engine!r}",
                                 hint="choose 'auto', 'bar', 'fast', or 'cs'")
        want = resolve_verbose(verbose, quiverlab.verbose)
        rec = trace if trace is not None else (Trace() if want else None)
        if self._route_to_cs(engine, auto_cs):
            from quiverlab.resolutions_cs.homology import cs_cohomology_dims
            table = cs_cohomology_dims(self, top, max_cells=max_cells, trace=rec)  # CS fills rec
        elif self._use_fast_engine(engine):
            from quiverlab.engine.adapter import engine_cohomology_dims
            if rec is not None:
                rec.append(Dispatch(
                    route="hanlab fast GF(p) rank",
                    reason="domain is a prime field; the exact mod-p rank engine applies",
                    n_relations=len(self.relations or ())))
            dims = engine_cohomology_dims(self, top, max_cells=max_cells)   # plain list[int]
            table = HHTable(dims, "HH^", repr(self).splitlines()[0],
                            engine="hanlab engine (F_p fast rank)")         # WRAP the list
        else:
            if rec is not None:
                rec.append(Dispatch(
                    route="normalized bar complex",
                    reason="domain is not a prime field; the exact bar oracle is used",
                    n_relations=len(self.relations or ())))
            table = hochschild_cohomology_dims(self, top, max_cells=max_cells, trace=rec)
        table.references = self.citations()   # FROZEN contract (family+engine keys); Task 11 must NOT change it
        # NOTE: the verbose worked-steps FILE is written in Task 11 (write_trace); here we only record/return.
        return table
```

`hochschild_homology` mirrors this exactly with `HH_`/`hochschild_homology_dims`/`cs_homology_dims`/`engine_homology_dims`: the same `("auto","bar","fast","cs")` guard raising `QuiverlabError`; the same `_route_to_cs` CS branch (`cs_homology_dims(self, top, max_cells=max_cells, trace=rec)`); the fast branch recording the fast `Dispatch` then wrapping `engine_homology_dims(self, top, max_cells=max_cells)` in `HHTable(dims, "HH_", repr(self).splitlines()[0], engine="hanlab engine (F_p fast rank)")`; the bar branch recording the bar `Dispatch` and calling `hochschild_homology_dims(self, top, max_cells=max_cells, trace=rec)`; then `table.references = self.citations()` and `return table`.

**Fast path is Dispatch-only, on purpose.** The `engine.adapter` fast GF(p) engine is left UNCHANGED (no new `trace=` kwarg): it works over several primes at once and surfaces no single exact matrix, so recording a per-degree matrix there would be dishonest. The engine-choice `Dispatch` (recorded above, in `algebra.py`) is enough for `table.references` and a short prime-field trace; the *rich* worked example (terms + matrices + ranks) is the exact bar/CS path, which is what the golden fixture (F5) and every worked-steps test exercise.

- [ ] **Step 4: Run the focused suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_emit.py tests/hochschild tests/resolutions_cs/test_dispatch.py tests/engine/test_adapter_dispatch.py tests/test_no_floats.py -q`
Expected: pass (bar path emits terms + ranks `[0,1,0]`, matrix for `d^1` = `[["0","0"],["2","0"]]`, dims `[2,1,1]`; prime-field path routes fast; `trace=[]` writes no file; existing hochschild tests unaffected). The added `tests/resolutions_cs/test_dispatch.py` and `tests/engine/test_adapter_dispatch.py` confirm the PRESERVED merged contracts in-task — `engine="cs"`/`auto_cs` CS routing, the unknown-engine `QuiverlabError`, and the fast-engine HHTable wrap — instead of only at Task 13.

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/hochschild/bar.py src/quiverlab/core/algebra.py tests/trace/test_emit.py
git commit -m "feat(trace): Hochschild engines emit ResolutionTerm/RankStep/Dispatch behind verbose/trace

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 9: Text renderer + the golden worked-steps trace for HH*(k[x]/(x²))

**Files:**
- Create: `src/quiverlab/trace/render_text.py`, `tests/trace/test_render_text.py`
- Create: `tests/trace/golden/hh_kx2.txt`

**Interfaces:**
- Consumes: all seven events; `trace.recorder.Trace`.
- Produces:
  - `render_text(events, title="", references=()) -> str` — a plain-text worked-steps document: title, the `Dispatch` (chosen resolution + why), per-degree `ResolutionTerm` ("term of degree n: cn generators"), each `RankStep` (matrix or elision note + "rank = k over <field>"), the **derived** resulting dimensions (`HH^n = collapsed_dim − rank_n − rank_{n-1}`, computed *from the recorded events*, not echoed), any `elision_notes`, and a **References** section listing the passed `references` keys. `DifferentialEvent`/`LiftStep`/`AmbiguityEvent` render as their own lines (for CS traces); `ReductionStep` renders compactly (capped at `TERMS_ELISION`).
  - `derive_dims(events) -> list[int]` (the binding computation reused by all three renderers).

- [ ] **Step 1: Write the failing test + golden file**

`tests/trace/golden/hh_kx2.txt` (exact expected text for F5; the citation line is spliced from Plan 06's `bar` registry entry at test time — see `_bar_reference()` below — so only the `__BAR_FORMATTED__` placeholder stands in for Plan 06's exact wording):

```
Worked steps: HH^ of Algebra of dimension 2 over CC (computing exactly in QQ)

Chosen resolution: normalized bar complex
  reason: domain is not a prime field; the exact bar oracle is used
  defining relations: 1

Degree 0: term with 2 generators (dim C = 2)
  differential d^0 (cochain), 2 x 2 over CC (computing exactly in QQ):
    [ 0  0 ]
    [ 0  0 ]
  rank = 0

Degree 1: term with 2 generators (dim C = 2)
  differential d^1 (cochain), 2 x 2 over CC (computing exactly in QQ):
    [ 0  0 ]
    [ 2  0 ]
  rank = 1

Degree 2: term with 2 generators (dim C = 2)
  differential d^2 (cochain), 2 x 2 over CC (computing exactly in QQ):
    [ 0  0 ]
    [ 0  0 ]
  rank = 0

Result: HH^0 = 2   HH^1 = 1   HH^2 = 1

References:
  [Hochschild1945] __BAR_FORMATTED__
```

(The `[Hochschild1945]` label is the stable BibTeX id the freshness gate pins to Plan 06's `bar` entry; the `__BAR_FORMATTED__` token is a placeholder the test substitutes with the `bar` entry's `.formatted` string from `bibliography()`, so Plan 06's exact wording is never duplicated here.)

`tests/trace/test_render_text.py`:

```python
"""Plain-text worked-steps renderer + the HH*(k[x]/x^2) golden trace (spec §3.8, §8).
The binding discipline: every claim in the rendered text equals what the engine
computed (recomputed independently here)."""
import pathlib

import quiverlab
from quiverlab import truncated_polynomial, CC
from quiverlab.trace.recorder import Trace
from quiverlab.trace.render_text import render_text, derive_dims
from quiverlab.trace.events import RankStep

GOLDEN = pathlib.Path(__file__).parent / "golden" / "hh_kx2.txt"


def _record():
    A = truncated_polynomial(2, field=CC)
    tr = Trace()
    table = A.hochschild_cohomology(2, trace=tr)
    return A, tr, table


def _bar_reference():
    """The bar engine's citation, taken from Plan 06's registry (never hardcoded):
    the `bar` registry key's (bibtex_key, formatted)."""
    bib = quiverlab.bibliography()
    entry = {e.key: e for e in bib}["bar"]
    return entry.bibtex_key, entry.formatted


def test_derive_dims_matches_computed():
    A, tr, table = _record()
    assert derive_dims(list(tr)) == table.dims == [2, 1, 1]


def test_golden_text_trace_matches_with_registry_reference():
    A, tr, table = _record()
    title = "HH^ of " + repr(A).splitlines()[0]
    bibtex_key, formatted = _bar_reference()
    assert bibtex_key == "Hochschild1945"          # stable BibTeX id Plan 06 backs `bar` with
    refs = ((bibtex_key, formatted),)
    text = render_text(list(tr), title=title, references=refs)
    # every line but the registry-sourced citation is golden-fixed; splice the
    # formatted line in from the registry so Plan 06's exact wording is not duplicated.
    expected = GOLDEN.read_text().replace("__BAR_FORMATTED__", formatted)
    assert text == expected
    assert "Hochschild" in formatted


def test_claims_equal_computed_ranks():
    A, tr, table = _record()
    ranks = [e.rank for e in tr if isinstance(e, RankStep)]
    assert ranks == [0, 1, 0]
    # the "rank = k" lines in the text are exactly these
    text = render_text(list(tr))
    assert text.count("rank = 0") == 2 and text.count("rank = 1") == 1


def test_elision_note_rendered(tmp_path):
    from quiverlab.trace.events import RankStep as RS
    ev = [RS(degree=5, side="cochain", nrows=30, ncols=30, rank=3, field="CC",
             matrix=None, elided=True, note="30x30 matrix over CC elided (> 400 cells); rank 3 recorded")]
    text = render_text(ev)
    assert "elided" in text and "rank = 3" in text
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_render_text.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'quiverlab.trace.render_text'`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/trace/render_text.py`:

```python
"""Plain-text worked-steps renderer (spec §3.8). The resulting dimensions are
DERIVED from the recorded ResolutionTerm + RankStep events (never echoed) -- that
is the binding discipline: the golden tests assert these derived numbers equal the
engine's own .dims, so a trace can never claim something the engine did not compute."""
from quiverlab.trace.events import (
    Dispatch, ResolutionTerm, RankStep, DifferentialEvent, LiftStep,
    AmbiguityEvent, ReductionStep,
)
from quiverlab.trace.recorder import TERMS_ELISION


def derive_dims(events):
    """HH^n / HH_n = collapsed_dim(n) - rank(n) - rank(n-1), from the events."""
    cn = {e.degree: e.collapsed_dim for e in events if isinstance(e, ResolutionTerm)}
    rk = {e.degree: e.rank for e in events if isinstance(e, RankStep)}
    dims = []
    for n in sorted(cn):
        dims.append(cn[n] - rk.get(n, 0) - rk.get(n - 1, 0))
    return dims


def _matrix_block(rs):
    if rs.elided or rs.matrix is None:
        return ["  (%s)" % rs.note]
    widths = [max(len(rs.matrix[i][j]) for i in range(rs.nrows)) for j in range(rs.ncols)]
    out = []
    for i in range(rs.nrows):
        row = "  ".join(rs.matrix[i][j].rjust(widths[j]) for j in range(rs.ncols))
        out.append("    [ %s ]" % row)
    return out


def render_text(events, title="", references=()):
    events = list(events)
    lines = []
    if title:
        lines.append("Worked steps: " + title)
        lines.append("")
    for e in events:
        if isinstance(e, Dispatch):
            lines.append("Chosen resolution: " + e.route)
            lines.append("  reason: " + e.reason)
            lines.append("  defining relations: %d" % e.n_relations)
            lines.append("")
    terms = {e.degree: e for e in events if isinstance(e, ResolutionTerm)}
    ranks = {e.degree: e for e in events if isinstance(e, RankStep)}
    for n in sorted(terms):
        t = terms[n]
        lines.append("Degree %d: term with %d generators (dim C = %d)"
                     % (n, t.n_generators, t.collapsed_dim))
        if n in ranks:
            rs = ranks[n]
            side = "d^%d" % n if rs.side == "cochain" else "b_%d" % n
            lines.append("  differential %s (%s), %d x %d over %s:"
                         % (side, rs.side, rs.nrows, rs.ncols, rs.field))
            lines.extend(_matrix_block(rs))
            lines.append("  rank = %d" % rs.rank)
        lines.append("")
    # CS symbolic differentials / lifts / ambiguities (present in CS traces)
    for e in events:
        if isinstance(e, DifferentialEvent):
            shown = e.terms[:TERMS_ELISION]
            lines.append("Symbolic differential (degree %d): %d term(s)%s"
                         % (e.degree, len(e.terms),
                            "" if len(e.terms) <= TERMS_ELISION
                            else " (%d shown)" % TERMS_ELISION))
        elif isinstance(e, LiftStep):
            lines.append("Lift/comparison step (degree %d): %s" % (e.degree, e.kind))
        elif isinstance(e, AmbiguityEvent):
            lines.append("Ambiguity chain (degree %d): %d word(s)"
                         % (e.degree, len(e.chain_words)))
    dims = derive_dims(events)
    if dims:
        kind = "HH^" if any(getattr(e, "side", "") == "cochain"
                            for e in events if isinstance(e, RankStep)) else "HH_"
        cells = "   ".join("%s%d = %d" % (kind, i, d) for i, d in enumerate(dims))
        lines.append("Result: " + cells)
        lines.append("")
    if references:
        lines.append("References:")
        for key, entry in references:
            lines.append("  [%s] %s" % (key, entry))
    return "\n".join(lines).rstrip("\n") + "\n"
```

(The golden `hh_kx2.txt` shows the exact spacing this function produces for F5: a 2-wide right-justified matrix, `rank = k` lines, the derived `Result: HH^0 = 2 ...`, and the single reference.)

- [ ] **Step 4: Run the focused suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_render_text.py tests/test_no_floats.py -q`
Expected: pass (golden text matches byte-for-byte; derived dims `[2,1,1]` == `.dims`; elision note rendered). If the byte-for-byte diff fails, adjust the golden file to the renderer's exact spacing — do NOT loosen the binding assertions (`derive_dims == .dims`, ranks `[0,1,0]`).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/trace/render_text.py tests/trace/test_render_text.py tests/trace/golden/hh_kx2.txt
git commit -m "feat(trace): plain-text renderer + golden HH*(k[x]/x^2) worked-steps trace

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 10: LaTeX→PDF (primary) + no-JS HTML fallback renderers, renderer selection, output path, printed one-liner

**Files:**
- Create: `src/quiverlab/trace/render_latex.py`, `src/quiverlab/trace/render_html.py`, `src/quiverlab/trace/writer.py`, `tests/trace/test_writer.py`

**Interfaces:**
- Consumes: `render_text.derive_dims`, all events, `Trace`.
- Produces:
  - `render_latex(events, title="", references=()) -> str` (standalone `article`; matrices as `pmatrix`; References as `thebibliography`).
  - `render_html(events, title="", references=()) -> str` — **NO JavaScript, NO external resource** (Marco's decision): a self-contained, offline `.html` whose formulas are shown as **plain TeX source** (matrices as `\begin{pmatrix}` text inside `<pre>`/`<code>`, NOT typeset); References as an ordered list. LaTeX→PDF is the primary math output; this HTML is the no-toolchain fallback.
  - `have_latex() -> str | None` (`"pdflatex"`/`"tectonic"` via `shutil.which`, else `None`).
  - `write_trace(events, table, algebra, kind, top, references=(), out_dir=None) -> str` — chooses PDF (compile LaTeX with the found engine) else HTML (loud one-line explanation), writes to `<out_dir or ./quiverlab_traces>/<kind>_<hash>.<ext>`, **prints** the spec §3.8 one-liner, returns the path.

- [ ] **Step 1: Write the failing test**

`tests/trace/test_writer.py`:

```python
"""LaTeX/HTML renderers + renderer selection + output-path contract (spec §3.8)."""
import pathlib

import pytest

from quiverlab import truncated_polynomial, CC
from quiverlab.trace.recorder import Trace
from quiverlab.trace import writer as W
from quiverlab.trace.render_latex import render_latex
from quiverlab.trace.render_html import render_html


def _events():
    A = truncated_polynomial(2, field=CC)
    tr = Trace()
    table = A.hochschild_cohomology(2, trace=tr)
    return A, list(tr), table


# A generic (bibtex_key, formatted) fixture; writer/renderer tests do not assert on
# any specific citation, so they stay decoupled from Plan 06's registry.
REFS = (("Refkey2020", "A. Author, A Journal 1 (2020), 1-2."),)


def test_latex_has_matrices_dims_and_bibliography():
    A, ev, table = _events()
    src = render_latex(ev, title="HH", references=REFS)
    assert r"\begin{pmatrix}" in src
    assert "HH^{0} = 2" in src or "HH^0 = 2" in src
    assert r"\begin{thebibliography}" in src and "Refkey2020" in src
    assert src.startswith(r"\documentclass")


def test_html_is_self_contained_no_js_tex_source():
    A, ev, table = _events()
    html = render_html(ev, title="HH", references=REFS)
    assert "<html" in html.lower()
    # Marco's decision: the HTML fallback is JavaScript-free and self-contained.
    assert "<script" not in html.lower() and "mathjax" not in html.lower()
    assert "polyfill.io" not in html and "jsdelivr" not in html
    # math is shown as readable TeX source (not typeset)
    assert r"\begin{pmatrix}" in html and r"\operatorname{rank}" in html
    assert "HH" in html and "Refkey2020" in html


def test_selection_prefers_pdf_when_toolchain_present(tmp_path, monkeypatch):
    A, ev, table = _events()
    monkeypatch.setattr(W, "have_latex", lambda: "tectonic")
    compiled = {}
    def fake_compile(tex, out_pdf, engine):
        pathlib.Path(out_pdf).write_bytes(b"%PDF-1.5 fake\n%%EOF\n")
        compiled["engine"] = engine
        return 1  # page count
    monkeypatch.setattr(W, "_compile_pdf", fake_compile)
    printed = {}
    monkeypatch.setattr("builtins.print", lambda *a, **k: printed.setdefault("line", " ".join(map(str, a))))
    path = W.write_trace(ev, table, algebra=A, kind="HH", top=2, references=REFS, out_dir=str(tmp_path))
    assert path.endswith(".pdf") and pathlib.Path(path).exists()
    assert compiled["engine"] == "tectonic"
    assert printed["line"].startswith("Worked steps: ") and ".pdf" in printed["line"]


def test_selection_falls_back_to_html_with_loud_message(tmp_path, monkeypatch):
    A, ev, table = _events()
    monkeypatch.setattr(W, "have_latex", lambda: None)   # no toolchain
    printed = {}
    monkeypatch.setattr("builtins.print", lambda *a, **k: printed.setdefault("line", " ".join(map(str, a))))
    path = W.write_trace(ev, table, algebra=A, kind="HH", top=2, references=REFS, out_dir=str(tmp_path))
    assert path.endswith(".html") and pathlib.Path(path).exists()
    line = printed["line"]
    # loud one-line explanation of BOTH facts: no toolchain -> HTML written
    assert line.startswith("Worked steps: ") and ".html" in line
    assert "no LaTeX toolchain found" in line and ("pdflatex" in line or "tectonic" in line)


def test_html_fallback_when_compile_fails(tmp_path, monkeypatch):
    """A toolchain IS present but compilation raises: the message must be honest
    ('compilation failed'), never 'no toolchain found'."""
    A, ev, table = _events()
    monkeypatch.setattr(W, "have_latex", lambda: "tectonic")
    def boom(tex, out_pdf, engine):
        raise RuntimeError("tectonic exploded")
    monkeypatch.setattr(W, "_compile_pdf", boom)
    printed = {}
    monkeypatch.setattr("builtins.print",
                        lambda *a, **k: printed.setdefault("line", " ".join(map(str, a))))
    path = W.write_trace(ev, table, algebra=A, kind="HH", top=2, references=REFS, out_dir=str(tmp_path))
    assert path.endswith(".html") and pathlib.Path(path).exists()
    line = printed["line"]
    assert "compilation failed" in line
    assert "no LaTeX toolchain found" not in line


def test_output_dir_and_filename_contract(tmp_path, monkeypatch):
    A, ev, table = _events()
    monkeypatch.setattr(W, "have_latex", lambda: None)
    monkeypatch.setattr("builtins.print", lambda *a, **k: None)
    path = W.write_trace(ev, table, algebra=A, kind="HH", top=2, references=REFS, out_dir=str(tmp_path))
    p = pathlib.Path(path)
    assert p.parent == tmp_path
    assert p.name.startswith("HH_") and p.suffix == ".html"  # <kind>_<hash>.<ext>
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_writer.py -q`
Expected: FAIL — modules `render_latex`/`render_html`/`writer` do not exist.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/trace/render_latex.py`:

```python
"""LaTeX worked-steps renderer -> standalone article compiled to PDF (spec §3.8).
Matrices as pmatrix; resulting dims derived from events; References as
thebibliography. Float-free: all numbers come from event fields (ints/strings)."""
from quiverlab.trace.events import Dispatch, ResolutionTerm, RankStep
from quiverlab.trace.render_text import derive_dims


def _pmatrix(rs):
    if rs.elided or rs.matrix is None:
        return r"\text{%s}" % rs.note
    rows = r" \\ ".join(" & ".join(rs.matrix[i][j] for j in range(rs.ncols))
                        for i in range(rs.nrows))
    return r"\begin{pmatrix} %s \end{pmatrix}" % rows


def render_latex(events, title="", references=()):
    events = list(events)
    out = [r"\documentclass{article}", r"\usepackage{amsmath}",
           r"\begin{document}", r"\section*{Worked steps: %s}" % _tex_escape(title)]
    for e in events:
        if isinstance(e, Dispatch):
            out.append(r"\noindent\textbf{Chosen resolution:} %s\\" % _tex_escape(e.route))
            out.append(r"\textit{%s}\\" % _tex_escape(e.reason))
            out.append(r"defining relations: %d" % e.n_relations)
    terms = {e.degree: e for e in events if isinstance(e, ResolutionTerm)}
    ranks = {e.degree: e for e in events if isinstance(e, RankStep)}
    for n in sorted(terms):
        t = terms[n]
        out.append(r"\subsection*{Degree %d}" % n)
        out.append(r"Term with %d generators ($\dim C = %d$)." % (t.n_generators, t.collapsed_dim))
        if n in ranks:
            rs = ranks[n]
            sym = "d^{%d}" % n if rs.side == "cochain" else "b_{%d}" % n
            out.append(r"\[ %s = %s \qquad \operatorname{rank} = %d \]"
                       % (sym, _pmatrix(rs), rs.rank))
    dims = derive_dims(events)
    if dims:
        kind = "HH^" if any(getattr(e, "side", "") == "cochain"
                            for e in events if isinstance(e, RankStep)) else "HH_"
        cells = r",\quad ".join(r"%s{%d} = %d" % (kind, i, d) for i, d in enumerate(dims))
        out.append(r"\subsection*{Result}")
        out.append(r"\[ %s \]" % cells)
    if references:
        out.append(r"\begin{thebibliography}{9}")
        for key, entry in references:
            out.append(r"\bibitem{%s} %s" % (key, _tex_escape(entry)))
        out.append(r"\end{thebibliography}")
    out.append(r"\end{document}")
    return "\n".join(out) + "\n"


def _tex_escape(s):
    for a, b in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"),
                 ("_", r"\_"), ("#", r"\#"), ("$", r"\$")):
        s = s.replace(a, b)
    return s
```

`src/quiverlab/trace/render_html.py`:

```python
"""HTML worked-steps renderer -- NO JavaScript, NO external resources (Marco's
decision). LaTeX->PDF (pdflatex/tectonic) is the PRIMARY, typeset math output;
this HTML is the self-contained, offline fallback for when no LaTeX toolchain is on
PATH. Math is shown as readable *TeX source* (inside <pre><code>), NOT typeset --
there is no MathJax, no CDN <script>, and no external <link>, so the file renders
identically in any browser with the network off. Float-free: all numbers come from
event fields (ints/strings)."""
from quiverlab.trace.events import Dispatch, ResolutionTerm, RankStep
from quiverlab.trace.render_text import derive_dims

# Inline-only styling (no external stylesheet); keeps the file fully self-contained.
_STYLE = ("<style>body{font-family:sans-serif}"
          "pre{background:#f4f4f4;padding:6px;overflow-x:auto}</style>")


def _pmatrix(rs):
    if rs.elided or rs.matrix is None:
        return r"\text{%s}" % rs.note
    rows = r" \\ ".join(" & ".join(rs.matrix[i][j] for j in range(rs.ncols))
                        for i in range(rs.nrows))
    return r"\begin{pmatrix} %s \end{pmatrix}" % rows


def _math(expr):
    """Show the TeX SOURCE as text (no MathJax); escape HTML metachars only."""
    return "<pre><code>%s</code></pre>" % _esc(expr)


def render_html(events, title="", references=()):
    events = list(events)
    body = ["<!doctype html><html><head><meta charset='utf-8'>", _STYLE,
            "<title>Worked steps: %s</title></head><body>" % _esc(title),
            "<h1>Worked steps: %s</h1>" % _esc(title),
            "<p><i>Math is shown as TeX source (no JavaScript); compile the PDF with "
            "pdflatex/tectonic for typeset output.</i></p>"]
    for e in events:
        if isinstance(e, Dispatch):
            body.append("<p><b>Chosen resolution:</b> %s<br><i>%s</i><br>"
                        "defining relations: %d</p>" % (_esc(e.route), _esc(e.reason), e.n_relations))
    terms = {e.degree: e for e in events if isinstance(e, ResolutionTerm)}
    ranks = {e.degree: e for e in events if isinstance(e, RankStep)}
    for n in sorted(terms):
        t = terms[n]
        body.append("<h2>Degree %d</h2><p>Term with %d generators (dim C = %d).</p>"
                    % (n, t.n_generators, t.collapsed_dim))
        if n in ranks:
            rs = ranks[n]
            sym = "d^{%d}" % n if rs.side == "cochain" else "b_{%d}" % n
            body.append(_math(r"%s = %s \qquad \operatorname{rank} = %d"
                              % (sym, _pmatrix(rs), rs.rank)))
    dims = derive_dims(events)
    if dims:
        kind = "HH^" if any(getattr(e, "side", "") == "cochain"
                            for e in events if isinstance(e, RankStep)) else "HH_"
        cells = ",\\quad ".join(r"%s{%d} = %d" % (kind, i, d) for i, d in enumerate(dims))
        body.append("<h2>Result</h2>" + _math(cells))
    if references:
        body.append("<h2>References</h2><ol>")
        for key, entry in references:
            body.append("<li>[%s] %s</li>" % (_esc(key), _esc(entry)))
        body.append("</ol>")
    body.append("</body></html>")
    return "\n".join(body) + "\n"


def _esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
```

`src/quiverlab/trace/writer.py`:

```python
"""Renderer selection + output-path contract + the printed one-liner (spec §3.8).

Selection: pdflatex or tectonic on PATH -> compile LaTeX to PDF; otherwise (or if a
found toolchain fails to compile) write the self-contained no-JS HTML (TeX source)
and print a LOUD, HONEST one-liner distinguishing "no toolchain found" from
"compilation failed". Output:
./quiverlab_traces/HHc_<hash>.<ext> (cohomology) / HHh_<hash>.<ext> (homology)
(Plan 09 collects the newest *.pdf, else *.html, from this directory -- the glob is
extension-based, so the safe stem does not affect it)."""
import hashlib
import pathlib
import shutil
import subprocess
import tempfile

from quiverlab.trace.render_text import derive_dims
from quiverlab.trace.render_latex import render_latex
from quiverlab.trace.render_html import render_html

# Filesystem-safe filename stems for the caret-bearing kinds (no "^" in a filename).
_SAFE_STEM = {"HH^": "HHc", "HH_": "HHh"}


def have_latex():
    for engine in ("tectonic", "pdflatex"):
        if shutil.which(engine):
            return engine
    return None


def _hash(algebra, kind, top):
    h = hashlib.sha1(("%s|%s|%s" % (repr(algebra), kind, top)).encode("utf-8"))
    return h.hexdigest()[:4]


def _compile_pdf(tex, out_pdf, engine):
    """Compile `tex` to `out_pdf` with `engine`; return the page count (best effort)."""
    with tempfile.TemporaryDirectory() as d:
        src = pathlib.Path(d) / "trace.tex"
        src.write_text(tex)
        if engine == "tectonic":
            cmd = ["tectonic", "-o", d, str(src)]
        else:
            cmd = ["pdflatex", "-interaction=nonstopmode", "-output-directory", d, str(src)]
        subprocess.run(cmd, cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        built = pathlib.Path(d) / "trace.pdf"
        shutil.copyfile(built, out_pdf)
        data = built.read_bytes()
        return max(data.count(b"/Type /Page") - data.count(b"/Type /Pages"), 1)


def write_trace(events, table, algebra, kind, top, references=(), out_dir=None):
    out = pathlib.Path(out_dir) if out_dir is not None else (pathlib.Path.cwd() / "quiverlab_traces")
    out.mkdir(parents=True, exist_ok=True)
    stem = "%s_%s" % (_SAFE_STEM.get(kind, kind), _hash(algebra, kind, top))
    title = "%s of %s" % (kind, repr(algebra).splitlines()[0])
    engine = have_latex()
    html_note = "no LaTeX toolchain found -- install pdflatex or tectonic for a PDF"
    if engine is not None:
        pdf = out / (stem + ".pdf")
        try:
            pages = _compile_pdf(render_latex(events, title=title, references=references),
                                 str(pdf), engine)
            print("Worked steps: %s (%d pp)" % (_rel(pdf), pages))
            return str(pdf)
        except Exception:
            # a toolchain WAS found but the compile failed: say so honestly; never
            # claim "no toolchain found" when one is on PATH.
            html_note = "LaTeX compilation failed (%s); wrote HTML fallback" % engine
    html = out / (stem + ".html")
    html.write_text(render_html(events, title=title, references=references))
    print("Worked steps: %s (HTML, no JavaScript; %s)" % (_rel(html), html_note))
    return str(html)


def _rel(p):
    try:
        return str(p.relative_to(pathlib.Path.cwd()))
    except ValueError:
        return str(p)
```

- [ ] **Step 4: Run the focused suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_writer.py tests/test_no_floats.py -q`
Expected: pass (LaTeX has pmatrix + dims + thebibliography; HTML has NO `<script>`/MathJax and shows the TeX source; selection prefers PDF when `have_latex` is stubbed, falls back to HTML with the loud two-fact message; filename `HH_<hash>.<ext>` in the chosen dir).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/trace/render_latex.py src/quiverlab/trace/render_html.py src/quiverlab/trace/writer.py tests/trace/test_writer.py
git commit -m "feat(trace): LaTeX/HTML renderers, PDF-or-HTML selection, output path + one-liner

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 11: References — provenance map, `result.references`, References section end-to-end

**Files:**
- Create: `src/quiverlab/trace/provenance.py`, `tests/trace/test_references.py`
- Modify: `src/quiverlab/core/algebra.py` (call `write_trace` when verbose; **do NOT reassign `.references`** — it stays `self.citations()` as set in Task 8)

**Interfaces:**
- Consumes: `quiverlab.bibliography()` (Plan 06 registry), the recorded `Dispatch`/`LiftStep` events.
- Produces:
  - `ENGINE_REFERENCES: dict[str, tuple[str,...]]` mapping an engine `route` to Plan 06's lowercase REGISTRY keys (`"normalized bar complex" -> ("bar",)`, `"hanlab fast GF(p) rank" -> ("bar",)`, `"bardzell" -> ("bardzell",)`, `"chouhy-solotar" -> ("chouhy_solotar",)`).
  - `references_for(events) -> tuple[str,...]` — REGISTRY keys implied by the trace's `Dispatch` (engine); operation keys (cup/bracket) are `()` until the operation trace lands. Deduplicated, stable order. This is what a result's `.references` stores.
  - `resolve_references(keys) -> tuple[tuple[str,str],...]` — `(bibtex_key, formatted)` pairs read from Plan 06's `bibliography()` by iterating it into a `{entry.key: entry}` map (no subscripting; `.keys` is a tuple, iteration yields entry views with `.key/.formatted/.bibtex_key`); a registry key absent from the bibliography raises a loud `KeyError` (the Task-1 gate keeps this from firing).
  - `Algebra.hochschild_cohomology/homology`: leave `table.references = self.citations()` (Task 8; the merged family+engine union) **unchanged** — do NOT overwrite it engine-only. When verbose and no explicit `trace=`, call `write_trace(..., references=resolve_references(references_for(rec)))`: the References SECTION is resolved from the trace's ENGINE provenance keys (it MAY additionally list the family keys already on `table.references`).

- [ ] **Step 1: Write the failing test**

`tests/trace/test_references.py`:

```python
"""Provenance -> References section, resolved through Plan 06's bibliography()."""
import pathlib

import quiverlab
from quiverlab import truncated_polynomial, CC
from quiverlab.trace.recorder import Trace
from quiverlab.trace.provenance import (
    references_for, resolve_references, ENGINE_REFERENCES,
)
from quiverlab.trace.events import Dispatch


def test_bar_engine_maps_to_registry_key():
    tr = [Dispatch(route="normalized bar complex", reason="", n_relations=1)]
    assert references_for(tr) == ("bar",)                # REGISTRY key, not a BibTeX id
    assert ENGINE_REFERENCES["bardzell"] == ("bardzell",)


def test_resolve_uses_plan06_registry():
    pairs = resolve_references(("bar",))                 # registry key in
    assert len(pairs) == 1
    bibtex_key, entry = pairs[0]                          # (bibtex_key, formatted) out
    assert bibtex_key == "Hochschild1945"               # backed by Plan 06's `bar` entry
    assert "Hochschild" in entry


def test_resolve_unknown_key_raises_loudly():
    import pytest
    with pytest.raises(KeyError):
        resolve_references(("not_a_registry_key",))


def test_result_object_carries_the_citations_union():
    A = truncated_polynomial(2, field=CC)
    table = A.hochschild_cohomology(2, trace=Trace())
    # .references is the merged family+engine union (== A.citations()); this no-family
    # fixture makes that exactly ("bar",). Task 11 must NOT overwrite it engine-only.
    assert table.references == A.citations() == ("bar",)


def test_verbose_run_writes_html_with_references(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Force the HTML path so the citation is greppable (a compiled PDF may compress it);
    # PDF selection is covered separately in test_writer.py.
    from quiverlab.trace import writer as W
    monkeypatch.setattr(W, "have_latex", lambda: None)
    quiverlab.verbose = True
    try:
        A = truncated_polynomial(2, field=CC)
        A.hochschild_cohomology(2)   # no explicit trace -> verbose auto-writes a file
    finally:
        quiverlab.verbose = False
    out = tmp_path / "quiverlab_traces"
    assert out.is_dir()
    files = list(out.glob("HHc_*.html"))
    assert files, "no worked-steps file written"
    # the References section shows the bibtex id resolved from Plan 06's `bar` entry
    assert "Hochschild1945" in files[0].read_text()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_references.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'quiverlab.trace.provenance'`.

- [ ] **Step 3: Write minimal implementation**

`src/quiverlab/trace/provenance.py`:

```python
"""Trace provenance: which literature a computation used, resolved through
Plan 06's bibliography() registry (spec §3.9; Plan-07 owns the RENDERING, Plan 06
owns the registry data). ENGINE_REFERENCES maps an engine `route` to Plan 06's
lowercase REGISTRY keys (never BibTeX ids); the Task-1 freshness gate asserts every
such key resolves. bibliography() returns a Bibliography whose `.keys` is a tuple
and whose iteration yields entry views (.key / .formatted / .bibtex_key / ...);
there is no subscripting and no .keys() method, so we build a key->entry map by
iterating."""
import quiverlab

ENGINE_REFERENCES = {
    "normalized bar complex": ("bar",),
    "hanlab fast GF(p) rank": ("bar",),
    "bardzell": ("bardzell",),
    "chouhy-solotar": ("chouhy_solotar",),
}

# Cup/bracket (Tamarkin-Tsygan) operation references are wired when the cup/bracket
# trace lands (Plan 04+) against Plan 06's operation REGISTRY keys; empty until then,
# so a verbose run never asks the registry for a key it does not have.
_OPERATION_REFERENCES = ()


def references_for(events):
    """Plan-06 REGISTRY keys implied by a trace: the engine (from every Dispatch
    whose route names a known engine) plus operation keys when LiftSteps are present.
    These registry keys are what a result's `.references` stores (they resolve to
    formatted citations at render time via resolve_references)."""
    from quiverlab.trace.events import Dispatch, LiftStep
    keys = []
    for e in events:
        if isinstance(e, Dispatch):
            for k in ENGINE_REFERENCES.get(e.route, ()):
                if k not in keys:
                    keys.append(k)
    if any(isinstance(e, LiftStep) for e in events):
        for k in _OPERATION_REFERENCES:
            if k not in keys:
                keys.append(k)
    return tuple(keys)


def _entries_by_key():
    """Map registry key -> entry view by iterating bibliography() (no subscripting)."""
    return {e.key: e for e in quiverlab.bibliography()}


def resolve_references(keys):
    """(bibtex_key, formatted) pairs for a tuple of REGISTRY keys, read from Plan 06's
    bibliography(). A registry key absent from the bibliography raises loudly (never
    silently dropped) -- the Task-1 freshness gate keeps this from firing."""
    by_key = _entries_by_key()
    pairs = []
    for k in keys:
        if k not in by_key:
            raise KeyError(
                "citation registry has no key %r (Plan 06 bibliography drift; update "
                "ENGINE_REFERENCES / the freshness gate)" % (k,))
        e = by_key[k]
        pairs.append((getattr(e, "bibtex_key", k), e.formatted))
    return tuple(pairs)
```

In `src/quiverlab/core/algebra.py`, finish the `hochschild_cohomology` body from Task 8. Task 8 already set `table.references = self.citations()` (the FROZEN family+engine union) and left a `# NOTE: the verbose worked-steps FILE is written in Task 11` placeholder before `return table`. Task 11 **must NOT reassign `table.references`** — the merged `.references = self.citations()` contract is pinned by `tests/citations/test_result_references.py` and `tests/families/test_acceptance.py`, and overwriting it engine-only regresses them. Replace **only that NOTE line** with the verbose write, whose References SECTION is resolved from the ENGINE provenance keys of the trace (not from, and without changing, `table.references`):

```python
        if want and trace is None and rec is not None:
            from quiverlab.trace.provenance import references_for, resolve_references
            from quiverlab.trace.writer import write_trace
            # References SECTION = engine keys implied by the trace's Dispatch, resolved
            # through bibliography(); table.references stays self.citations() (untouched).
            write_trace(list(rec), table, algebra=self, kind="HH^", top=top,
                        references=resolve_references(references_for(rec)))
        return table
```

Mirror the same tail in `hochschild_homology` (with `kind="HH_"`); it too must leave `table.references = self.citations()` untouched.

- [ ] **Step 4: Run the focused suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_references.py tests/trace/test_freshness_gate.py tests/test_no_floats.py -q`
Expected: pass (`.references == ("bar",)` and `resolve_references(("bar",))` yields the `Hochschild1945` entry; a verbose run writes a file whose References section names Hochschild; the freshness gate still green — every provenance key resolves in `bibliography()`).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/trace/provenance.py tests/trace/test_references.py src/quiverlab/core/algebra.py
git commit -m "feat(trace): provenance map + result.references + References section via bibliography()

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 12: Renderer parity — same events → consistent claims across text/HTML/LaTeX

**Files:**
- Create: `tests/trace/test_parity.py`

**Interfaces:**
- Consumes: `render_text`, `render_html`, `render_latex`, `derive_dims`.
- Produces: no source — the cross-renderer binding lock.

- [ ] **Step 1: Write the failing test**

`tests/trace/test_parity.py`:

```python
"""Renderer parity (spec §8): the SAME events produce CONSISTENT claims in all
three renderers -- identical resulting dims and identical per-degree ranks. The
binding discipline holds regardless of output format."""
import re

from quiverlab import truncated_polynomial, CC
from quiverlab.trace.recorder import Trace
from quiverlab.trace.events import RankStep
from quiverlab.trace.render_text import render_text, derive_dims
from quiverlab.trace.render_html import render_html
from quiverlab.trace.render_latex import render_latex

# A generic (bibtex_key, formatted) fixture: parity is about the renderers agreeing,
# not about any specific citation, so this stays decoupled from Plan 06's registry.
REFS = (("Refkey2020", "A. Author, A Journal 1 (2020), 1-2."),)


def _events():
    A = truncated_polynomial(2, field=CC)
    tr = Trace()
    table = A.hochschild_cohomology(2, trace=tr)
    return list(tr), table


def _dims_in(s):
    return [int(m) for m in re.findall(r"HH[\^_]?\{?\d+\}?\s*=\s*(\d+)", s)]


def test_all_three_renderers_agree_on_dims():
    ev, table = _events()
    txt = render_text(ev, title="t", references=REFS)
    html = render_html(ev, title="t", references=REFS)
    tex = render_latex(ev, title="t", references=REFS)
    assert derive_dims(ev) == table.dims == [2, 1, 1]
    assert _dims_in(txt) == [2, 1, 1]
    assert _dims_in(html) == [2, 1, 1]
    assert _dims_in(tex) == [2, 1, 1]


def test_all_three_render_the_same_ranks():
    ev, table = _events()
    ranks = [e.rank for e in ev if isinstance(e, RankStep)]
    assert ranks == [0, 1, 0]
    for render in (render_text, render_html, render_latex):
        s = render(ev, title="t", references=REFS)
        # every recorded rank appears as a "rank = k" / "rank} = k" claim
        found = [int(m) for m in re.findall(r"rank[^\d=]*=\s*(\d+)", s)]
        assert found == ranks


def test_all_three_carry_the_reference():
    ev, table = _events()
    for render in (render_text, render_html, render_latex):
        assert "Refkey2020" in render(ev, title="t", references=REFS)
```

- [ ] **Step 2: Run test to verify it fails, then passes**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_parity.py -q`
Expected: PASS (all renderers already exist from Tasks 9–11; this task locks their agreement). If any renderer disagrees on dims/ranks, the bug is in that renderer — fix it, do not weaken the parity assertion.

- [ ] **Step 3: (no implementation — parity lock)**

If Step 2 fails, apply superpowers:systematic-debugging to the disagreeing renderer.

- [ ] **Step 4: Run the focused suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_parity.py tests/test_no_floats.py -q`
Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add tests/trace/test_parity.py
git commit -m "test(trace): renderer parity lock (text/HTML/LaTeX agree on dims and ranks)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 13: Acceptance — end-to-end, internals chapter, README, frozen statement, background full suite

**Files:**
- Create: `tests/trace/test_acceptance.py`
- Create: `docs/internals/11-viz-trace.md`
- Modify: `docs/internals/README.md` (add the chapter to the index)
- Modify: `README.md` (add a "Draw it / read the worked steps" section)

**Interfaces:**
- Consumes: the whole public surface — `A.draw`, `A.tikz`, `A.hochschild_cohomology(..., verbose=True)`, `quiverlab.verbose`.
- Produces: the Plan-07 exit criterion + the frozen statement for Plans 08/09.

- [ ] **Step 1: Write the failing test**

`tests/trace/test_acceptance.py`:

```python
"""Plan 07 acceptance (spec §3.7, §3.8, §3.9, D9): one algebra, drawn, TikZ'd, and
computed with a worked-steps document whose claims equal the computed values."""
import pathlib

import quiverlab
from quiverlab import Quiver, truncated_polynomial, CC


def test_end_to_end_draw_tikz_and_worked_steps(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    A = truncated_polynomial(2, field=CC)

    # draw + tikz share the exact layout
    fig = A.draw(file="A.svg")
    assert (tmp_path / "A.svg").exists() and fig.__class__.__name__ == "Figure"
    assert r"\begin{tikzpicture}" in A.tikz()

    # verbose default (D9) writes a worked-steps file; claims == computed
    quiverlab.verbose = True
    try:
        table = A.hochschild_cohomology(2)
    finally:
        quiverlab.verbose = False
    assert table.dims == [2, 1, 1]
    assert table.references == A.citations() == ("bar",)   # merged citations union (no family stamp here)
    out = tmp_path / "quiverlab_traces"
    files = list(out.glob("HHc_*"))                  # .pdf or .html per toolchain
    assert files, "verbose run must write a worked-steps file (D9)"


def test_per_call_verbose_false_writes_nothing(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    quiverlab.verbose = True
    try:
        truncated_polynomial(2, field=CC).hochschild_cohomology(2, verbose=False)
    finally:
        quiverlab.verbose = False
    assert not (tmp_path / "quiverlab_traces").exists()
```

- [ ] **Step 2: Run test to verify it passes**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest tests/trace/test_acceptance.py -q`
Expected: PASS if Tasks 1–12 are complete. If it fails, the failure localizes the missing piece; fix it before proceeding.

- [ ] **Step 3: Write the internals chapter**

`docs/internals/11-viz-trace.md`:

```markdown
# 11 — Visualization and worked-steps traces

## The mathematics

A quiver `Q` is a directed graph: vertices, and arrows (directed edges, possibly
loops or parallel). An algebra `kQ/I` is that graph plus a list of relations
(parallel-path linear combinations). To *look* at it we draw the graph and print
the relations beneath. To *learn from a computation* — a Hochschild cohomology,
say — we record the resolution chosen, each term, each differential as a matrix
over the working field, each rank, and the resulting dimensions, then lay them out
as a short worked example. Nothing here is new mathematics; it is the same
computation the engine runs, written down so a human can read it.

## How it is represented

**Layout.** `viz.layout.layout(quiver, relations)` returns a `LayoutData`: a dict
`positions` from each vertex to an exact `(x, y)` (`x` an integer column depth, `y`
a `Fraction` row), a tuple of `EdgeRoute`s (straight or parallel, with a `Fraction`
bend), a tuple of `LoopRoute`s (integer base angle), and the relation strings.
Coordinates are `int`/`Fraction` on purpose: they never touch algebra, but keeping
them exact means the float-ban gate covers viz with no exemption, and the layout is
golden-testable to the last coordinate. For the commuting square
`kQ/(a*b - c*d)` the positions are

    1 -> (0, 0)   2 -> (1, 1/2)   3 -> (1, -1/2)   4 -> (2, 0)

**Trace events.** A computation records a flat list of typed events (all plain
dataclasses): a `Dispatch` (which resolution, and why), one `ResolutionTerm` per
degree (its generator count), one `RankStep` per degree (the differential matrix
over the field, or a one-line elision note above 400 cells, plus its rank), and —
for the Chouhy–Solotar engine and the cup/bracket operations — `AmbiguityEvent`,
`DifferentialEvent`, and `LiftStep`. The list is capped at 5000 events so a deep
computation cannot blow up memory.

## How the computation runs

1. **Layering** (`viz.layout.layer`). Compute strongly-connected components
   (iterative Tarjan), condense them, and take the longest path on the condensation
   — so a cycle or loop collapses to one column and every vertex gets an integer
   depth. Vertices in a column are centered on half-integer rows.
2. **Drawing** (`viz.draw.draw_quiver`). Circles at `positions`, `FancyArrowPatch`
   for arrows (parallel bundles fan out via `ConnectionStyle.Arc3(rad=<Fraction>)`),
   `Arc` for loops (integer angle), the relations as a text block below. `tikz`
   emits the *identical* layout as `\node`/`\draw`, so a paper and a screen agree.
3. **Recording** (`trace.recorder.Trace`). `A.hochschild_cohomology(top,
   verbose=…, trace=…)` resolves verbosity (per-call overrides the global
   `quiverlab.verbose`, default `True`), records the engine-choice `Dispatch`, and
   runs the engine with the recorder; the bar engine appends a `ResolutionTerm` and
   a `RankStep` at each degree.
4. **Rendering** (`trace.writer.write_trace`). If `pdflatex` or `tectonic` is on
   `PATH`, compile the LaTeX rendering to `./quiverlab_traces/HHc_<hash>.pdf` and
   print `Worked steps: quiverlab_traces/HHc_<hash>.pdf (N pp)`; otherwise write
   the self-contained, JavaScript-free HTML rendering (math as TeX source) and print
   a loud one-line explanation (no toolchain found, or — if one was found — that
   compilation failed). The resulting
   dimensions in every rendering are *derived from the recorded ranks*
   (`HH^n = dim C^n − rank_n − rank_{n-1}`), so a golden test can assert the
   document's claims equal the engine's own `.dims`.
5. **References.** The engine's `route` maps (via `trace.provenance`) to Plan 06
   REGISTRY keys (e.g. `bar`); those resolve through Plan 06's `bibliography()`
   (`.keys` tuple + entry-view iteration exposing `.key/.formatted/.bibtex_key`)
   into `(bibtex_key, formatted)` lines in the document's **References** section.
   The result itself carries the merged `table.references = self.citations()` (the
   family + engine key union), which Plan 07 does not modify.

## A worked micro-example — HH*(k[x]/(x²)) over ℂ

`A = truncated_polynomial(2, field=CC)` has dimension 2, so the normalized bar
complex has `C^n = 2` for every `n`. Running `A.hochschild_cohomology(2)` records,
per degree, the 2×2 coboundary matrix and its rank:

    d^0 = [[0,0],[0,0]]  rank 0
    d^1 = [[0,0],[2,0]]  rank 1
    d^2 = [[0,0],[0,0]]  rank 0

so `HH^0 = 2−0−0 = 2`, `HH^1 = 2−1−0 = 1`, `HH^2 = 2−0−1 = 1`: dims `[2, 1, 1]`.
Over `GF(2)` the `2` becomes `0`, `d^1` has rank 0, and the dims are `[2, 2, 2]` —
the classic characteristic pathology, visible directly in the traced matrix. These
numbers were produced by running the code (they are the golden trace).

## Where to look in the code

| concept | file | function / class |
|---|---|---|
| layered layout, exact coords | `src/quiverlab/viz/layout.py` | `layout`, `layer`, `LayoutData` |
| matplotlib rendering | `src/quiverlab/viz/draw.py` | `draw_quiver`, `Algebra.draw` |
| TikZ rendering | `src/quiverlab/viz/tikz.py` | `tikz_quiver`, `Algebra.tikz` |
| event taxonomy | `src/quiverlab/trace/events.py` | `Dispatch`, `ResolutionTerm`, `RankStep`, `DifferentialEvent`, `LiftStep`, `AmbiguityEvent`, `ReductionStep` |
| recorder + elision | `src/quiverlab/trace/recorder.py` | `Trace`, `rankstep`, `resolve_verbose` |
| engine emission | `src/quiverlab/hochschild/bar.py` | `hochschild_cohomology_dims` |
| renderers | `src/quiverlab/trace/render_text.py` / `render_latex.py` / `render_html.py` | `render_text`, `render_latex`, `render_html`, `derive_dims` |
| selection + output path | `src/quiverlab/trace/writer.py` | `write_trace`, `have_latex` |
| provenance + References | `src/quiverlab/trace/provenance.py` | `references_for`, `resolve_references` |
```

Add to `docs/internals/README.md`'s chapter index a line: `11. Visualization and worked-steps traces (viz + trace)`.

- [ ] **Step 4: Update the README**

Append to `README.md` (before the license line):

```markdown
## Draw it, and read the worked steps

```python
from quiverlab import Quiver, CC

Q = Quiver(vertices=[1, 2, 3, 4],
           arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
A = Q.algebra(relations=["a*b - c*d"], field=CC)

A.draw(file="square.svg")     # matplotlib PNG/SVG: loops, parallels, relations below
print(A.tikz())               # same layout, paste-into-paper TikZ

A.hochschild_cohomology(2)    # writes quiverlab_traces/HHc_<hash>.pdf (or .html) and
                              # prints: Worked steps: quiverlab_traces/HHc_3f2a.pdf (N pp)
```

Worked-steps documents are on by default (`quiverlab.verbose = True`); every claim
in them is a golden-file-tested equality with the value the engine computed. Turn
them off per call (`A.hochschild_cohomology(2, verbose=False)`) or globally
(`quiverlab.verbose = False`). PDFs need `pdflatex` or `tectonic` on `PATH`;
otherwise a self-contained, JavaScript-free HTML document (math shown as TeX source)
is written with a one-line note.
```

- [ ] **Step 5: Run the full suite ONCE in the background, await, then commit**

Launch the full suite as a tracked background job (per the suite pattern):
`NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 /Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python -m pytest -q`
Await its completion. Expected: **all pass** — the full Plan-01…06 suite plus all of `tests/viz/` and `tests/trace/`, with `tests/test_no_floats.py` green (viz is float-free; no gate exemption). Then also confirm the pure-Python path once: `QUIVERLAB_NO_NUMBA=1 NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 …/.venv/bin/python -m pytest -q`. Only after the background full suite is green:

```bash
git add tests/trace/test_acceptance.py docs/internals/11-viz-trace.md docs/internals/README.md README.md
git commit -m "feat(viz,trace): Plan 07 acceptance — draw/tikz + worked-steps traces, internals ch.11

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Frozen statement for Plans 08 and 09

Plans 08 (release/docs) and 09 (webapp) consume the following **verbatim** (do not rename, move, or change types without a coordinated interface bump):

**Trace output-path contract.** Verbose computations write **one** worked-steps file to `./quiverlab_traces/` (relative to the current working directory), named `<stem>_<hash>.<ext>` where `stem` is a filesystem-safe token — `HHc` for cohomology, `HHh` for homology (the caret-bearing `HH^`/`HH_` kinds are mapped to these safe stems; §3.8's `HH_A_<hash>` is illustrative) — `hash` is a 4-hex short hash, and `ext` is `pdf` when a LaTeX toolchain (`pdflatex`/`tectonic`) is on `PATH`, otherwise `html`. Exactly one file per computation; `.pdf` is preferred, `.html` is the fallback. Plan 09's `_collect_trace` globs `quiverlab_traces/*.pdf` (then `*.html`) and copies the **newest by mtime** — so the directory, the extension precedence (pdf ≻ html), and "newest is the just-computed one" are the load-bearing guarantees; the exact filename (and its safe stem) is not. The child process must run with `cwd` = the job's artifact dir (Plan 09 already `os.chdir`es), so the directory lands where the collector looks.

**Printed one-liner (spec §3.8).** PDF: `Worked steps: quiverlab_traces/HHc_<hash>.pdf (N pp)`. HTML fallback, **no toolchain**: `Worked steps: quiverlab_traces/HHc_<hash>.html (HTML, no JavaScript; no LaTeX toolchain found -- install pdflatex or tectonic for a PDF)`. HTML fallback, **compile failed** (toolchain present but errored): `Worked steps: quiverlab_traces/HHc_<hash>.html (HTML, no JavaScript; LaTeX compilation failed (<engine>); wrote HTML fallback)`. Always the loud, honest explanation — never "no toolchain" when one was found. The HTML is self-contained and JavaScript-free (math as TeX source); LaTeX→PDF is the primary typeset output.

**Verbose flag.** `quiverlab.verbose: bool` (module attribute, default `True`, D9). Plan 09 sets it per request (`ql.verbose = bool(req.artifacts.pdf)`) and restores it. Per-call `verbose=` on `hochschild_cohomology`/`hochschild_homology` overrides the global; passing an explicit `trace=` list populates events **without** writing a file.

**Event taxonomy (single import surface `quiverlab.trace.events`).** Seven plain dataclasses: `Dispatch(route, reason, n_relations)`, `ReductionStep(word, rule_lead, before, after)`, `AmbiguityEvent(degree, chain_words)`, `ResolutionTerm(degree, n_generators, collapsed_dim)`, `DifferentialEvent(degree, chain, terms)`, `LiftStep(degree, kind, detail)`, `RankStep(degree, side, nrows, ncols, rank, field, matrix, elided, note)`. The first six are re-exported from their home modules (`groebner.events`, `resolutions_cs.trace`); `RankStep` is defined in `trace.events`. `Dispatch` labels both the construction route (`monomial`/`groebner`) and the engine choice (`route` = engine name).

**Renderer selection.** `pdflatex` or `tectonic` on `PATH` ⇒ LaTeX→PDF (primary math output); else a self-contained, JavaScript-free HTML showing TeX source (loud note). `quiverlab.trace.render_text`/`render_latex`/`render_html` each take `(events, title=, references=)`; `derive_dims(events)` is the binding computation shared by all three (resulting dims are always derived from the recorded ranks, never echoed).

**Result references.** `HHTable` carries a post-hoc `.references: tuple[str,...]` of Plan 06 **registry** keys (lowercase, e.g. `("bar",)`); `quiverlab.trace.provenance.resolve_references` turns them into `(bibtex_key, formatted)` pairs by **iterating** Plan 06's `bibliography()` — a `Bibliography` with a `.keys` **tuple** and entry-view iteration exposing `.key/.formatted/.bibtex_key` (no `.keys()` method, no subscripting) — which the renderers list in the References section. The engine→registry-key map is `quiverlab.trace.provenance.ENGINE_REFERENCES`.

**Viz surface.** `A.draw(file=None) -> Figure` (PNG/SVG by extension; built via `matplotlib.figure.Figure` + `FigureCanvasAgg`, so it never mutates the global backend — no `matplotlib.use()`), `A.tikz() -> str` (same coordinates), `quiverlab.viz.layout.layout(quiver, relations=()) -> LayoutData` (exact `int`/`Fraction`; axis limits snapped to ints via `math.floor/ceil` because matplotlib≥3.11 rejects `Fraction` bounds). `matplotlib>=3.7` is a **hard** dependency (in `[project].dependencies`, not an extra) — **Plan 08's packaging task must keep it there**. `print(A)` shows dimension (line 0, unchanged), basis, vertices, arrows, and relations.
```
