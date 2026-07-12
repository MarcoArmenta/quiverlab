# quiverlab Plan 02 — hanlab Engine Port Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** quiverlab gains hanlab's battle-tested F_p compute engine — fast exact linear algebra (pure/sparse/numba, equality-gated), bar homology/cohomology with fast ranks, bimodule/twisted coefficients, the Tamarkin–Tsygan calculus (cup, cap, Gerstenhaber bracket), cyclic homology, minimal/periodic/Bardzell resolutions with guards and checkpointed deepening, and Cartan/Coxeter/Nakayama machinery — wired behind Plan 01's public API with automatic dispatch and bar-oracle cross-checks.

**Architecture:** hanlab modules are copied (never moved) from the read-only bank into `src/quiverlab/engine/` as an internal package, keeping their module names, with three mechanical transformations: package-relative imports, removal of `__main__` demo blocks, and float-literal eradication (the AST gate is non-negotiable). An adapter converts Plan-01 `Algebra` objects over `PrimeField(p)` into engine algebras; `hochschild_cohomology/homology` gain an `engine=` parameter with `"auto"` dispatch (fast for GF(p), pure bar otherwise) and a cross-check battery pins engine ≡ bar on the Plan-01 anchor zoo. New public invariants (Cartan/Coxeter, all fields, from quiver provenance; Nakayama/Frobenius/symmetric over GF(p) via the engine) land with loud errors where a field is not yet supported.

**Tech Stack:** Python ≥ 3.10; hard deps become `numpy>=1.21`, `sympy>=1.12`; optional extra `fast = ["numba>=0.64"]`; pytest.

## Global Constraints

- **The bank is READ-ONLY.** Source of truth: `/Users/marco/Desktop/HomologicalNetworks/HomologicalAlgebra/HansConjecture/` (below: `$BANK`). Files are COPIED out with `cp`; never write, move, or delete anything under `$BANK`. Both hanlab and quiverlab are MIT © 2026 Marco Armenta.
- **Attribution header** (exact text, first lines of every ported file, above the module docstring):
  ```python
  # Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
  # github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.
  # Mechanical changes only: package-relative imports, __main__ blocks removed,
  # float literals eradicated (quiverlab AST gate), env guard renamed.
  ```
- **Import rewrite table** (apply to every ported module AND test): `from X import ...` / `import X` where `X ∈ {hh_engine, scan2, scan3, coxeter, coxeter2, coxeter_spectrum, bimodule, tt_calculus, cyclic, resolutions, resolutions_minimal, resolutions_periodic, resolutions_bardzell, deepen, linalg_fast, _kernels}` → `from quiverlab.engine.X import ...` / `from quiverlab.engine import X`. No `sys.path` hacks survive the port (hanlab's `__init__` hack and `tests/conftest.py` are NOT copied).
- **Strip every `if __name__ == "__main__":` block** from ported modules (dev scratch; also removes coxeter.py's float-bearing demo print at bank line ~438).
- **Float-literal eradication** — the Task-1 AST gate scans `src/` and must stay green. Known instances and their exact rewrites are itemized in the tasks below; for ANY other float/complex literal or `float()` call the gate flags, apply the same style of integer-exact rewrite and record it in the task report. Floats inside docstrings/comments are fine (the gate ignores them). Wall-clock values returned by `time.time()` are permitted (no literal, no `float()` call; they never touch algebra).
- **Env guard rename:** `HANLAB_NO_NUMBA` → `QUIVERLAB_NO_NUMBA` (single read site: `_kernels.py` bank line 28; update the docstrings that mention it and any test that sets it). `USE_KERNELS` keeps its name.
- **Engine is internal.** Nothing from `quiverlab.engine` is added to the top-level `__all__`. Public surface grows only via the Task-13 methods.
- **Resource discipline:** export `NUMBA_NUM_THREADS=2` and `OMP_NUM_THREADS=2` in every test command of this plan. No parallel processes.
- **Venv:** ALWAYS `/Users/marco/Desktop/quiverlab/.venv/bin/python` (system python is 3.8). Task 1 installs `[dev,fast]` so numba parity tests genuinely run; Task 14 re-runs the whole suite with `QUIVERLAB_NO_NUMBA=1` (pure path).
- **Ported tests** go to `tests/engine/` (same filenames), with imports rewritten per the table. A bank test file is ported **iff** after rewriting it imports only ported modules; the expected inclusion list is in Task 14's checklist and each port task names its files. Excluded bank tests (they import labdb/open_zoo/resolutions_cs/scan/scan_open/cluster): `test_labdb, test_open_zoo_broaden, test_cs_resolution, test_cs_induced_action, test_cluster_runner, test_scan_drivers`.
- **Engine defaults stay hanlab's** (e.g. p = 32003 as the char-0 proxy prime; primes {2,3,5} in scans). Do not "improve" ported numerics.
- Run tests with `... .venv/bin/python -m pytest -q` from the repo root; green at every commit. Conventional commits; append your harness's standard co-author trailer.
- **Not ported in Plan 02** (deferred): `scan.py`, `scan_open.py`, `open_zoo.py`, `labdb.py`, `reduction_algebra.py`, `resolutions_cs.py`, `cluster/` (Plans 04/06), and from `coxeter_spectrum.py` the two numeric functions `spectral_radius` / `mahler_measure` (exact sympy reimplementations arrive in Plan 05). JSON artifacts (`coxeter_results.json`, `coxeter2_nakayama.json`, `coxeter_spectrum_report.json`, `open_zoo_catalog*.json`) are copied ONLY if a ported test reads them (check with grep before copying; record which).

---

### Task 1: Branch, dependencies, engine skeleton

**Files:**
- Modify: `pyproject.toml`
- Create: `src/quiverlab/engine/__init__.py`
- Test: existing suite must stay green (82 passed)

**Interfaces:**
- Consumes: nothing new.
- Produces: `quiverlab.engine` package; numpy as a hard dep; `[fast]` extra; branch `plan-02-hanlab-port`.

- [ ] **Step 1: Create the branch**

Run: `cd /Users/marco/Desktop/quiverlab && git checkout -b plan-02-hanlab-port`
Expected: `Switched to a new branch 'plan-02-hanlab-port'`.

- [ ] **Step 2: Update pyproject dependencies**

In `pyproject.toml`, change the `dependencies` line and add the `fast` extra so the relevant sections read exactly:

```toml
dependencies = ["numpy>=1.21", "sympy>=1.12"]

[project.optional-dependencies]
dev = ["pytest>=8"]
fast = ["numba>=0.64"]
```

- [ ] **Step 3: Create the engine package**

`src/quiverlab/engine/__init__.py`:

```python
"""quiverlab.engine: the F_p compute engine, ported from hanlab.

Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.

INTERNAL API. The engine computes over prime fields F_p with exact int64
arithmetic (structure constants held unreduced, reduced mod p at rank time;
p = 32003 is the char-0 proxy). Public quiverlab entry points dispatch here
automatically for algebras over GF(p); everything else uses the pure Plan-01
paths. Set QUIVERLAB_NO_NUMBA=1 to force the pure-Python kernels.

Modules keep their hanlab development names (hh_engine, scan3, coxeter, ...).
"""
```

- [ ] **Step 4: Reinstall with extras and verify baseline**

Run: `cd /Users/marco/Desktop/quiverlab && .venv/bin/python -m pip install -q -e '.[dev,fast]' && NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q`
Expected: `82 passed` (numpy/numba install cleanly; nothing else changes).

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml src/quiverlab/engine/__init__.py
git commit -m "chore: engine package skeleton, numpy hard dep, fast extra"
```

---

### Task 2: Port the F_p linear-algebra stack (_kernels, linalg_fast)

**Files:**
- Create: `src/quiverlab/engine/_kernels.py`, `src/quiverlab/engine/linalg_fast.py` (copies from `$BANK/hanlab/`), `tests/engine/__init__.py` (empty), `tests/engine/test_kernels.py`, `tests/engine/test_linalg_fast.py` (copies from `$BANK/tests/`)
- Test: the two ported test files

**Interfaces:**
- Consumes: numpy, optional numba.
- Produces: `quiverlab.engine._kernels` (numba kernels + pure twins, `USE_KERNELS`, `QUIVERLAB_NO_NUMBA` guard) and `quiverlab.engine.linalg_fast` (`sparse_rank_mod_p`, `rank_mod_p_auto`) — exactly hanlab's semantics.

- [ ] **Step 1: Copy the four files**

```bash
BANK=/Users/marco/Desktop/HomologicalNetworks/HomologicalAlgebra/HansConjecture
cp "$BANK/hanlab/_kernels.py" src/quiverlab/engine/_kernels.py
cp "$BANK/hanlab/linalg_fast.py" src/quiverlab/engine/linalg_fast.py
cp "$BANK/tests/test_kernels.py" tests/engine/test_kernels.py
cp "$BANK/tests/test_linalg_fast.py" tests/engine/test_linalg_fast.py
```

- [ ] **Step 2: Apply the mechanical transformations**

1. Prepend the attribution header (Global Constraints) to both source files.
2. Import rewrites per the table (e.g. `import _kernels` → `from quiverlab.engine import _kernels`; `from hh_engine import rank_mod_p` → `from quiverlab.engine.hh_engine import rank_mod_p` — note: `linalg_fast` imports `rank_mod_p` from `hh_engine`, which does not exist yet; see Step 3).
3. Env rename in `_kernels.py` (bank line ~28): `os.environ.get("HANLAB_NO_NUMBA", ...)` → `os.environ.get("QUIVERLAB_NO_NUMBA", ...)`; update every docstring/comment mention in both files and both tests; rewrite any test that sets `HANLAB_NO_NUMBA` to set `QUIVERLAB_NO_NUMBA`.
4. **Float-literal rewrite** in `linalg_fast.py` (bank line ~80): the signature default `max_density=0.05` and its comparison must become integer-exact. Replace the parameter with `max_density_permille=50` and rewrite the density test from the form `nnz / size <= max_density` to `1000 * nnz <= max_density_permille * size` (adjust to the file's actual variable names — the inequality must be algebraically identical). Update the docstring accordingly. If any caller in later-ported modules passes `max_density=...`, rewrite that call site the same way (grep when porting).
5. Strip `__main__` blocks if present.

- [ ] **Step 3: Break the forward dependency for this task only**

`linalg_fast.rank_mod_p_auto` falls back to `hh_engine.rank_mod_p`, which is ported in Task 3. To keep Task 2 self-contained and tested, move that import INSIDE the function body (lazy import), exactly:

```python
def rank_mod_p_auto(M, p, sparse_min_size=200_000, max_density_permille=50):
    from quiverlab.engine.hh_engine import rank_mod_p
    ...
```

(keep the rest of the function verbatim). If `tests/engine/test_linalg_fast.py` imports `hh_engine` symbols directly, defer THOSE test functions to Task 3 by moving them into a clearly-marked block at the bottom of the file with a comment `# enabled in Task 3 (hh_engine)` and a module-level `pytest.importorskip("quiverlab.engine.hh_engine")` — do NOT delete them.

- [ ] **Step 4: Run the ported tests, then the full suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/engine/test_kernels.py tests/engine/test_linalg_fast.py -q` then `... -m pytest -q`
Expected: ported tests pass (kernel/pure parity across primes {2,3,5,32003}); full suite green; AST gate green (no float literals survived).

- [ ] **Step 5: Commit**

```bash
git add src/quiverlab/engine tests/engine
git commit -m "feat(engine): port F_p kernel stack (_kernels, linalg_fast) with QUIVERLAB_NO_NUMBA guard"
```

---

### Task 3: Port hh_engine + scan2 (engine core and algebra builders)

**Files:**
- Create: `src/quiverlab/engine/hh_engine.py`, `src/quiverlab/engine/scan2.py`, `tests/engine/test_engine_validation.py`, `tests/engine/test_multivertex_engine.py`, `tests/engine/test_scan2_builders.py`
- Modify: `tests/engine/test_linalg_fast.py` (re-enable the deferred block from Task 2, if any)

**Interfaces:**
- Consumes: `engine._kernels`, `engine.linalg_fast`.
- Produces: `engine.hh_engine.Algebra` (int64 structure-constant algebra: the engine currency), `rank_mod_p`, `cn_basis`, `differential_matrix`, `hochschild_homology_dims`, `check_associative`, plus hanlab's closed-form builders (`truncated_polynomial`, `two_gen_local`, ...); `engine.scan2` builders (`tensor_product`, `trivial_extension`, `triangular_extension`, module builders) used as fixtures by many bank tests.

- [ ] **Step 1: Copy and transform**

```bash
cp "$BANK/hanlab/hh_engine.py" src/quiverlab/engine/hh_engine.py
cp "$BANK/hanlab/scan2.py"     src/quiverlab/engine/scan2.py
cp "$BANK/tests/test_engine_validation.py"  tests/engine/
cp "$BANK/tests/test_multivertex_engine.py" tests/engine/
cp "$BANK/tests/test_scan2_builders.py"     tests/engine/
```

Apply the mechanical transformations (attribution, import rewrites, `__main__` strip) to all five. Audit `scan2.py` for float literals with the gate's own logic before running (`.venv/bin/python -c "import ast,sys; [print(v) for v in __import__('ast').walk(ast.parse(open('src/quiverlab/engine/scan2.py').read())) if isinstance(v, ast.Constant) and isinstance(v.value, (float, complex))]"`) and rewrite integer-exactly anything found, recording it in your report. Same audit for `hh_engine.py`.

- [ ] **Step 2: Restore the Task-2 deferred tests**

If Task 2 deferred `hh_engine`-dependent test functions in `tests/engine/test_linalg_fast.py`, remove the `importorskip` and the marker comment now.

- [ ] **Step 3: Run ported tests, then full suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/engine -q` then `... -m pytest -q`
Expected: all green (engine validation cross-checks the int64 bar complex against known values; multivertex covers the unit-adapted trick; builders build).

- [ ] **Step 4: Commit**

```bash
git add src/quiverlab/engine tests/engine
git commit -m "feat(engine): port hh_engine core and scan2 builders"
```

---

### Task 4: Port scan3 (fast-rank cohomology, complexity, quantum CIs)

**Files:**
- Create: `src/quiverlab/engine/scan3.py`, `tests/engine/test_cohomology_fast_rank.py`, `tests/engine/test_complexity.py`

**Interfaces:**
- Consumes: `engine.hh_engine`, `engine.linalg_fast`.
- Produces: `engine.scan3.hochschild_cohomology_dims`, `cochain_basis`, `coboundary_matrix`, `complexity_of`, `quantum_ci`.

- [ ] **Step 1: Copy and transform**

```bash
cp "$BANK/hanlab/scan3.py" src/quiverlab/engine/scan3.py
cp "$BANK/tests/test_cohomology_fast_rank.py" tests/engine/
cp "$BANK/tests/test_complexity.py"           tests/engine/
```

Mechanical transformations. **Known wrinkle:** `scan3.py` contains an inline (function-body) import of an excluded module (the recon grep matched it on a non-top-level line). Find it (`grep -n "scan\|labdb\|open_zoo" src/quiverlab/engine/scan3.py`), and: if it sits in a demo/driver function or `__main__` block, delete that function/block; if a genuinely-needed helper imports it, report NEEDS_CONTEXT with the exact lines — do not guess. Same treatment for `test_complexity.py`'s excluded-module import if present (the recon flagged it): port only the test functions whose imports are within the ported set; leave a `# not ported: depends on <module> (Plan 06)` comment listing dropped test names.

- [ ] **Step 2: Run ported tests, then full suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/engine -q` then `... -m pytest -q`
Expected: green.

- [ ] **Step 3: Commit**

```bash
git add src/quiverlab/engine tests/engine
git commit -m "feat(engine): port scan3 cohomology, complexity, quantum CIs"
```

---

### Task 5: Port bimodule (coefficients, dual, twisted)

**Files:**
- Create: `src/quiverlab/engine/bimodule.py`, `tests/engine/test_bimodule_coefficients.py`, `tests/engine/test_twisted_bimodule.py`

**Interfaces:**
- Consumes: `engine.hh_engine`.
- Produces: `engine.bimodule.{hochschild_homology_with_coefficients, hochschild_cohomology_with_coefficients, regular_bimodule, dual_bimodule, twisted_bimodule, cn_basis_coeff, differential_matrix_coeff}`.

- [ ] **Step 1: Copy, transform, resolve fixture imports**

```bash
cp "$BANK/hanlab/bimodule.py" src/quiverlab/engine/bimodule.py
cp "$BANK/tests/test_bimodule_coefficients.py" tests/engine/
cp "$BANK/tests/test_twisted_bimodule.py"      tests/engine/
```

Mechanical transformations. `test_bimodule_coefficients.py` imports scan2 builders as fixtures (recon-flagged) — those now resolve to `quiverlab.engine.scan2` (ported in Task 3); rewrite accordingly. Any residual import of a non-ported module: drop only the affected test functions with the `# not ported:` comment convention.

- [ ] **Step 2: Run, then full suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/engine -q` then `... -m pytest -q`
Expected: green.

- [ ] **Step 3: Commit**

```bash
git add src/quiverlab/engine tests/engine
git commit -m "feat(engine): port bimodule coefficients (regular, dual, twisted)"
```

---

### Task 6: Port coxeter (F_p linalg + induced actions + Nakayama)

**Files:**
- Create: `src/quiverlab/engine/coxeter.py`, `tests/engine/test_automorphism.py`, `tests/engine/test_nakayama_automorphism.py`, `tests/engine/test_linalg_modp.py`, `tests/engine/test_basis_transport.py` (iff import-closed, see Step 1)

**Interfaces:**
- Consumes: `engine.hh_engine`, `engine.scan3`.
- Produces: `engine.coxeter.{rref_mod_p, nullspace_mod_p, solve_mod_p, inverse_mod_p, quotient_induced, sigma_chain_matrix, induced_on_HH_homology, induced_on_HH_cohomology, is_identity, nakayama_automorphism, is_frobenius, frobenius_form}`.

- [ ] **Step 1: Copy and transform**

```bash
cp "$BANK/hanlab/coxeter.py" src/quiverlab/engine/coxeter.py
cp "$BANK/tests/test_automorphism.py"          tests/engine/
cp "$BANK/tests/test_nakayama_automorphism.py" tests/engine/
cp "$BANK/tests/test_linalg_modp.py"           tests/engine/
```

Mechanical transformations. The bank file's `__main__`/demo block contains the float-bearing diagnostic print (bank line ~438) — the strip rule removes it; verify no float literal survives. `test_nakayama_automorphism.py` was recon-flagged for an excluded-module import: apply the fixture-resolution rule (scan2 now exists; anything else → drop affected functions with the comment convention). Check `test_basis_transport.py` imports (`grep "^from\|^import" $BANK/tests/test_basis_transport.py`): port it iff closed over ported modules.

- [ ] **Step 2: Run, then full suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/engine -q` then `... -m pytest -q`
Expected: green.

- [ ] **Step 3: Commit**

```bash
git add src/quiverlab/engine tests/engine
git commit -m "feat(engine): port coxeter (mod-p linalg, induced HH actions, Nakayama)"
```

---

### Task 7: Port tt_calculus (cup, cap, Gerstenhaber bracket)

**Files:**
- Create: `src/quiverlab/engine/tt_calculus.py`, `tests/engine/test_tt_calculus.py`, `tests/engine/test_gerstenhaber.py`

**Interfaces:**
- Consumes: `engine.hh_engine` (`cn_basis`, `differential_matrix`), `engine.scan3` (`cochain_basis`, `coboundary_matrix`), `engine.coxeter`.
- Produces: `engine.tt_calculus.{cup_product_matrix, cap_product_matrix, gerstenhaber_bracket_matrix}` + cochain helpers — **the operations no other system ships.**

- [ ] **Step 1: Copy and transform**

```bash
cp "$BANK/hanlab/tt_calculus.py" src/quiverlab/engine/tt_calculus.py
cp "$BANK/tests/test_tt_calculus.py"  tests/engine/
cp "$BANK/tests/test_gerstenhaber.py" tests/engine/
```

Mechanical transformations (the arXiv float in the module docstring is prose — leave it).

- [ ] **Step 2: Run, then full suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/engine -q` then `... -m pytest -q`
Expected: green — including the graded-commutativity, Leibniz, and graded-Jacobi law tests that make this port trustworthy.

- [ ] **Step 3: Commit**

```bash
git add src/quiverlab/engine tests/engine
git commit -m "feat(engine): port Tamarkin-Tsygan calculus (cup, cap, Gerstenhaber bracket)"
```

---

### Task 8: Port cyclic homology

**Files:**
- Create: `src/quiverlab/engine/cyclic.py`, `tests/engine/test_cyclic_homology.py`

**Interfaces:**
- Consumes: `engine.hh_engine`.
- Produces: `engine.cyclic.{connes_B_matrix, cyclic_homology_dims, check_mixed_complex}`.

- [ ] **Step 1: Copy, transform, run, commit**

```bash
cp "$BANK/hanlab/cyclic.py" src/quiverlab/engine/cyclic.py
cp "$BANK/tests/test_cyclic_homology.py" tests/engine/
```

Mechanical transformations. Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/engine -q` then full suite. Expected: green. Commit:

```bash
git add src/quiverlab/engine tests/engine
git commit -m "feat(engine): port cyclic homology (Connes B, mixed complex)"
```

---

### Task 9: Port the resolution suite (protocol, minimal, periodic, deepen)

**Files:**
- Create: `src/quiverlab/engine/resolutions.py`, `src/quiverlab/engine/resolutions_minimal.py`, `src/quiverlab/engine/resolutions_periodic.py`, `src/quiverlab/engine/deepen.py`, and in `tests/engine/`: `test_resolution_protocol.py`, `test_resolution_contract.py`, `test_minimal_resolution.py`, `test_minimal_memory_guard.py`, `test_periodic_resolution.py`, `test_periodic_symmetric_family.py`, `test_deepen.py`, plus `test_truncated_resolution.py`, `test_error_paths.py`, `test_recheck.py` iff import-closed (check each with grep; recon flagged the latter two for excluded imports — apply the fixture-resolution/drop rules)

**Interfaces:**
- Consumes: `engine.hh_engine`, `engine._kernels`.
- Produces: `engine.resolutions.Resolution` (protocol), `engine.resolutions_minimal.{minimal_resolution, minimal_homology_dims, hochschild_dimension}` (+ the `_init/_advance_resolution` stepper, memory/term guards), `engine.resolutions_periodic` (periodicity detection), `engine.deepen.deepen` (checkpointed resumable driver).

- [ ] **Step 1: Copy and transform**

```bash
cp "$BANK/hanlab/resolutions.py"          src/quiverlab/engine/resolutions.py
cp "$BANK/hanlab/resolutions_minimal.py"  src/quiverlab/engine/resolutions_minimal.py
cp "$BANK/hanlab/resolutions_periodic.py" src/quiverlab/engine/resolutions_periodic.py
cp "$BANK/hanlab/deepen.py"               src/quiverlab/engine/deepen.py
cp "$BANK"/tests/test_resolution_protocol.py "$BANK"/tests/test_resolution_contract.py \
   "$BANK"/tests/test_minimal_resolution.py "$BANK"/tests/test_minimal_memory_guard.py \
   "$BANK"/tests/test_periodic_resolution.py "$BANK"/tests/test_periodic_symmetric_family.py \
   "$BANK"/tests/test_deepen.py tests/engine/
```

Mechanical transformations, plus these **exact float-literal rewrites in `deepen.py`**:
- bank line ~28: `_TIME_PRED_FACTOR = 2.0` → `_TIME_PRED_FACTOR = 2`
- bank line ~137: `per_degree[-1].get("secs", 0.0) if per_degree else 0.0` → `per_degree[-1].get("secs", 0) if per_degree else 0`
(Runtime values stay floats from `time.time()` — permitted; only literals go.)

**`resolutions_periodic.py`** was recon-flagged for an inline excluded-module import: locate it; if it lives in a demo/driver function or `__main__`, delete that block; otherwise NEEDS_CONTEXT with the lines. `deepen.py` checkpoint paths must not be hardcoded to bank locations — grep for absolute paths; if found, parameterize to an argument with the current behavior as default under the caller's CWD, and note it in the report.

- [ ] **Step 2: Run, then full suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/engine -q` then `... -m pytest -q`
Expected: green (memory-guard and checkpoint tests included).

- [ ] **Step 3: Commit**

```bash
git add src/quiverlab/engine tests/engine
git commit -m "feat(engine): port resolution suite (protocol, minimal, periodic, deepen)"
```

---

### Task 10: Port Bardzell (monomial minimal resolution) + QPA reference battery

**Files:**
- Create: `src/quiverlab/engine/resolutions_bardzell.py`, `tests/engine/test_bardzell_resolution.py`, `tests/engine/test_qpa_reference_validation.py` (iff import-closed — this is the bank's literature/QPA cross-check battery and SHOULD come; if it imports an excluded module, port the closed subset and report what was dropped)

**Interfaces:**
- Consumes: `engine.resolutions.Resolution`.
- Produces: `engine.resolutions_bardzell.{BardzellResolution, MonomialPresentation}` — deep exact homology for monomial algebras (bank record: degree 1702).

- [ ] **Step 1: Copy, transform, run, commit**

```bash
cp "$BANK/hanlab/resolutions_bardzell.py" src/quiverlab/engine/resolutions_bardzell.py
cp "$BANK/tests/test_bardzell_resolution.py" tests/engine/
cp "$BANK/tests/test_qpa_reference_validation.py" tests/engine/   # iff import-closed
```

Mechanical transformations. Run tests/engine then full suite: green. Commit:

```bash
git add src/quiverlab/engine tests/engine
git commit -m "feat(engine): port Bardzell resolution and QPA reference battery"
```

---

### Task 11: Port coxeter2 + the exact subset of coxeter_spectrum

**Files:**
- Create: `src/quiverlab/engine/coxeter2.py`, `src/quiverlab/engine/coxeter_spectrum.py`, and in `tests/engine/`: `test_dynkin_zoo.py`, `test_theoremB_convention.py`, `test_theoremB_multivertex.py`, `test_coxeter_spectrum.py` (exact-function tests only)

**Interfaces:**
- Consumes: `engine.hh_engine`, `engine.scan3`, `engine.coxeter`, `engine.bimodule`, sympy.
- Produces: `engine.coxeter2.{cartan_from_raw, coxeter_polynomial_from_cartan, coxeter_element_order, charpoly_of_induced, linear_path_algebra, quiver_path_algebra, dynkin_quiver, cyclic_nakayama}`; `engine.coxeter_spectrum` (exact functions only: `is_cyclotomic_product`, `cartan_of_quiver`, `trivial_extension_cartan`, `star_quiver`, `nakayama_charpoly_hh`, `column_degeneration`, `dual_column_action`).

- [ ] **Step 1: Copy and transform**

```bash
cp "$BANK/hanlab/coxeter2.py"         src/quiverlab/engine/coxeter2.py
cp "$BANK/hanlab/coxeter_spectrum.py" src/quiverlab/engine/coxeter_spectrum.py
cp "$BANK"/tests/test_dynkin_zoo.py "$BANK"/tests/test_theoremB_convention.py \
   "$BANK"/tests/test_theoremB_multivertex.py "$BANK"/tests/test_coxeter_spectrum.py tests/engine/
```

Mechanical transformations, plus the **numeric-function excision in `coxeter_spectrum.py`**:
1. DELETE the functions `spectral_radius` and `mahler_measure` entirely (bank lines ~105-125, the `1.0` literals — exact sympy reimplementations are Plan 05's job). Leave a module-docstring note: `spectral_radius/mahler_measure: not ported; exact reimplementations arrive with the invariants plan (Plan 05).`
2. Bank line ~340: `det=int(round(float(sp.Matrix(CT.tolist()).det())))` — if the containing function is ported, rewrite exactly to `det=int(sp.Matrix(CT.tolist()).det())` (the determinant of an integer matrix is an exact sympy Integer; `round(float(...))` was a legacy cast). If the containing function is a report/driver that depends on deleted functions, delete it too and note it.
3. In `test_coxeter_spectrum.py`, delete only the test functions that call the deleted functions, with the `# not ported: exact reimplementation in Plan 05` comment convention. Everything else stays.

- [ ] **Step 2: Run, then full suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/engine -q` then `... -m pytest -q`
Expected: green; AST gate green (the deleted numeric functions were the last float carriers).

- [ ] **Step 3: Commit**

```bash
git add src/quiverlab/engine tests/engine
git commit -m "feat(engine): port coxeter2 and the exact coxeter_spectrum subset"
```

---

### Task 12: The adapter and engine dispatch

**Files:**
- Create: `src/quiverlab/engine/adapter.py`, `tests/engine/test_adapter_dispatch.py`
- Modify: `src/quiverlab/core/algebra.py` (`hochschild_cohomology`/`hochschild_homology` grow `engine="auto"`), `src/quiverlab/hochschild/table.py` (no change expected — `engine` label already a constructor arg)

**Interfaces:**
- Consumes: Plan-01 `Algebra` (domain, T, unit, unit_adapted), `PrimeField`, engine `hh_engine.Algebra`, `scan3.hochschild_cohomology_dims`, `hh_engine.hochschild_homology_dims`, `HHTable`, `FieldError`.
- Produces: `to_engine(A) -> engine Algebra` (GF(p) only, loud otherwise); `A.hochschild_cohomology(top, engine="auto"|"bar"|"fast")` and same for homology; `HHTable.engine` reports which engine ran.

- [ ] **Step 1: Write the failing tests**

`tests/engine/test_adapter_dispatch.py`:

```python
import pytest
from quiverlab import CC, GF, Quiver, truncated_polynomial
from quiverlab.errors import FieldError


def _zoo(field):
    Q1 = Quiver(vertices=[1], arrows={"x": (1, 1)})
    Q2 = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)})
    Q3 = Quiver(vertices=[1, 2], arrows={})
    return [
        Q1.algebra(relations=["x^2"], field=field),
        Q1.algebra(relations=["x^4"], field=field),
        Q2.algebra(relations=["a*b"], field=field),
        Q2.algebra(relations=[], field=field),
        Q3.algebra(relations=[], field=field),
    ]


@pytest.mark.parametrize("p", [2, 3, 5])
def test_engine_matches_bar_cohomology_and_homology(p):
    for A in _zoo(GF(p)):
        slow_c = A.hochschild_cohomology(4, engine="bar").dims
        fast_c = A.hochschild_cohomology(4, engine="fast").dims
        assert fast_c == slow_c, f"HH^ mismatch over GF({p}) on {A!r}"
        slow_h = A.hochschild_homology(4, engine="bar").dims
        fast_h = A.hochschild_homology(4, engine="fast").dims
        assert fast_h == slow_h, f"HH_ mismatch over GF({p}) on {A!r}"


def test_auto_dispatch_picks_fast_for_gfp_and_bar_for_cc():
    A = truncated_polynomial(2, field=GF(2))
    B = truncated_polynomial(2, field=CC)
    assert "fast" in A.hochschild_cohomology(2).engine
    assert "bar" in B.hochschild_cohomology(2).engine
    assert A.hochschild_cohomology(2).dims == [2, 2, 2]
    assert B.hochschild_cohomology(2).dims == [2, 1, 1]


def test_fast_engine_refuses_non_prime_fields_loudly():
    B = truncated_polynomial(2, field=CC)
    C = truncated_polynomial(2, field=GF(4))
    for A in (B, C):
        with pytest.raises(FieldError):
            A.hochschild_cohomology(2, engine="fast")


def test_gf4_still_works_via_bar_auto():
    A = truncated_polynomial(2, field=GF(4))
    assert A.hochschild_cohomology(2).dims == [2, 2, 2]
    assert "bar" in A.hochschild_cohomology(2).engine
```

- [ ] **Step 2: Run to verify failure**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/engine/test_adapter_dispatch.py -q`
Expected: FAIL — `hochschild_cohomology() got an unexpected keyword argument 'engine'`.

- [ ] **Step 3: Implement the adapter**

`src/quiverlab/engine/adapter.py`:

```python
"""Bridge between quiverlab's domain-generic Algebra and the hanlab engine.

The engine computes over prime fields F_p with numpy int64 structure
constants. `to_engine` converts a quiverlab Algebra over GF(p) (PrimeField);
every other domain is refused loudly (spec: the fast engine is a GF(p)
accelerator, the pure bar path serves all fields)."""
import numpy as np

from quiverlab.errors import FieldError
from quiverlab.fields.primefield import PrimeField


def to_engine(A):
    dom = A.domain
    if not isinstance(dom, PrimeField):
        raise FieldError(
            f"the fast engine computes over prime fields only; this algebra is over {dom.name}",
            hint="use engine='bar' (any field), or construct the algebra over GF(p)",
        )
    from quiverlab.engine.hh_engine import Algebra as EngineAlgebra

    m = A.dim
    T = np.zeros((m, m, m), dtype=np.int64)
    for i in range(m):
        for j in range(m):
            vec = A.T[i][j]
            for t in range(m):
                T[i, j, t] = int(vec[t])
    unit = np.array([int(c) for c in A.unit], dtype=np.int64)
    return EngineAlgebra(m, T, unit)


def engine_cohomology_dims(A, top):
    from quiverlab.engine.scan3 import hochschild_cohomology_dims
    E = to_engine(A.unit_adapted())
    return hochschild_cohomology_dims(E, top, A.domain.p)


def engine_homology_dims(A, top):
    from quiverlab.engine.hh_engine import hochschild_homology_dims
    E = to_engine(A.unit_adapted())
    return hochschild_homology_dims(E, top, A.domain.p)
```

**Signature note (resolve at implementation time, do not guess):** open the ported `hh_engine.py` and `scan3.py` and match the REAL constructor and dims-function signatures — hanlab's `Algebra` may take `(m, T, unit, name=...)` or keyword forms, and the dims functions may take `N`/`primes=[...]`/`p` in a different shape (they may also return per-prime dicts rather than lists — normalize to a plain `list[int]` of dimensions for degrees `0..top` inside these two wrappers). The Step-1 equality tests against the pure bar are the correctness gate; make them pass by adjusting ONLY these wrappers, never the ported engine.

- [ ] **Step 4: Wire the dispatch into Algebra**

In `src/quiverlab/core/algebra.py`, replace the two Hochschild methods with:

```python
    def hochschild_cohomology(self, top, max_cells=4_000_000, engine="auto"):
        """Dimensions of HH^0..HH^top, exact. engine: 'auto' (fast over GF(p),
        bar otherwise), 'bar' (pure, any field), 'fast' (GF(p) only, loud otherwise)."""
        from quiverlab.hochschild.bar import hochschild_cohomology_dims
        from quiverlab.hochschild.table import HHTable

        if engine not in ("auto", "bar", "fast"):
            from quiverlab.errors import QuiverlabError
            raise QuiverlabError(f"unknown engine {engine!r}",
                                 hint="choose 'auto', 'bar', or 'fast'")
        use_fast = engine == "fast" or (
            engine == "auto" and type(self.domain).__name__ == "PrimeField"
        )
        if use_fast:
            from quiverlab.engine.adapter import engine_cohomology_dims
            dims = engine_cohomology_dims(self, top)
            return HHTable(dims, "HH^", repr(self).splitlines()[0],
                           engine="hanlab engine (F_p fast rank)")
        return hochschild_cohomology_dims(self, top, max_cells=max_cells)

    def hochschild_homology(self, top, max_cells=4_000_000, engine="auto"):
        """Dimensions of HH_0..HH_top, exact. Same engine semantics as cohomology."""
        from quiverlab.hochschild.bar import hochschild_homology_dims
        from quiverlab.hochschild.table import HHTable

        if engine not in ("auto", "bar", "fast"):
            from quiverlab.errors import QuiverlabError
            raise QuiverlabError(f"unknown engine {engine!r}",
                                 hint="choose 'auto', 'bar', or 'fast'")
        use_fast = engine == "fast" or (
            engine == "auto" and type(self.domain).__name__ == "PrimeField"
        )
        if use_fast:
            from quiverlab.engine.adapter import engine_homology_dims
            dims = engine_homology_dims(self, top)
            return HHTable(dims, "HH_", repr(self).splitlines()[0],
                           engine="hanlab engine (F_p fast rank)")
        return hochschild_homology_dims(self, top, max_cells=max_cells)
```

(Use `isinstance` with a lazy `PrimeField` import instead of the `type(...).__name__` comparison if circular imports permit — prefer `from quiverlab.fields.primefield import PrimeField` at method level and `isinstance(self.domain, PrimeField)`; the `__name__` fallback is only for an import-cycle emergency, and if you use it, say so in your report.)

Also confirm the existing pure paths still label their `HHTable.engine` as the bar complex (Plan-01 default `"normalized bar complex"`) — the dispatch tests distinguish engines by substring `"fast"` vs `"bar"`.

- [ ] **Step 5: Run the new tests, then the full suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/engine/test_adapter_dispatch.py -q` then `... -m pytest -q`
Expected: all green — the 5-algebra × {2,3,5} × {HH^, HH_} equality battery is the port's central correctness gate.

- [ ] **Step 6: Commit**

```bash
git add src/quiverlab/engine/adapter.py src/quiverlab/core/algebra.py tests/engine/test_adapter_dispatch.py
git commit -m "feat: engine adapter and auto dispatch for GF(p) Hochschild"
```

---

### Task 13: Public invariants (Cartan/Coxeter for all fields; engine-backed extras for GF(p))

**Files:**
- Create: `src/quiverlab/invariants/__init__.py`, `src/quiverlab/invariants/cartan.py`, `tests/invariants/__init__.py` (empty), `tests/invariants/test_cartan.py`, `tests/invariants/test_engine_backed.py`
- Modify: `src/quiverlab/core/algebra.py` (add methods), `src/quiverlab/__init__.py` (no new exports needed — methods only)

**Interfaces:**
- Consumes: `Algebra` provenance (`.quiver`, `.basis_labels`), sympy, `engine.coxeter` (`nakayama_automorphism`, `is_frobenius`), `engine.cyclic`, adapter.
- Produces: `A.cartan_matrix() -> list[list[int]]`, `A.coxeter_matrix() -> list[list[int]]`, `A.coxeter_polynomial() -> sympy.Poly` (any field, from provenance); `A.cyclic_homology(top) -> HHTable`, `A.nakayama_automorphism()`, `A.is_frobenius() -> bool`, `A.is_symmetric() -> bool` (GF(p) via engine; loud `FieldError` naming the phase plan otherwise).

- [ ] **Step 1: Write the failing tests**

`tests/invariants/test_cartan.py`:

```python
import sympy
import pytest
from quiverlab import CC, GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.errors import QuiverlabError


def test_cartan_kA2():
    A = linear_path_algebra(2)
    # C[i][j] = dim e_i A e_j = number of basis paths v_i -> v_j
    assert A.cartan_matrix() == [[1, 1], [0, 1]]


def test_cartan_field_independent():
    for field in (CC, GF(2), GF(5)):
        Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)})
        A = Q.algebra(relations=["a*b"], field=field)
        assert A.cartan_matrix() == [[1, 1, 0], [0, 1, 1], [0, 0, 1]]


def test_cartan_dual_numbers():
    assert truncated_polynomial(2).cartan_matrix() == [[2]]


def test_coxeter_polynomial_A2():
    # Phi = -C^{-T} C for kA_2 has characteristic polynomial t^2 + t + 1 (Coxeter number 3)
    A = linear_path_algebra(2)
    t = sympy.Symbol("t")
    assert sympy.expand(A.coxeter_polynomial().as_expr() - (t**2 + t + 1)) == 0


def test_cartan_without_provenance_fails_loudly():
    # a bare structure-constant algebra has no path basis: the Cartan matrix must
    # refuse loudly with an actionable hint, never guess
    from quiverlab import Algebra
    T = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]
    A = Algebra.from_structure_constants(T, unit=[1, 0], field=CC)
    with pytest.raises(QuiverlabError):
        A.cartan_matrix()
```

`tests/invariants/test_engine_backed.py`:

```python
import pytest
from quiverlab import CC, GF, truncated_polynomial
from quiverlab.errors import FieldError


def test_cyclic_homology_gfp():
    A = truncated_polynomial(2, field=GF(3))
    t = A.cyclic_homology(4)
    assert t.kind == "HC_"
    assert len(t.dims) == 5 and all(isinstance(d, int) for d in t.dims)


def test_cyclic_homology_cc_loud():
    with pytest.raises(FieldError):
        truncated_polynomial(2, field=CC).cyclic_homology(4)


def test_symmetric_dual_numbers_gfp():
    A = truncated_polynomial(2, field=GF(5))
    assert A.is_frobenius() is True
    assert A.is_symmetric() is True


def test_frobenius_cc_loud():
    with pytest.raises(FieldError):
        truncated_polynomial(2, field=CC).is_frobenius()
```

- [ ] **Step 2: Run to verify failure**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/invariants -q`
Expected: FAIL — missing methods.

- [ ] **Step 3: Implement the provenance-based Cartan/Coxeter module**

`src/quiverlab/invariants/cartan.py`:

```python
"""Cartan and Coxeter data from the quiver presentation (any field).

For a monomial bound quiver algebra with basis the irreducible paths,
C[i][j] = dim e_i A e_j = #(basis paths from vertex i to vertex j) -- an
integer matrix independent of the ground field. The Coxeter matrix is
Phi = -C^{-T} C (requires C invertible over Q, e.g. finite global dimension);
the Coxeter polynomial is charpoly(Phi)."""
import sympy

from quiverlab.errors import QuiverlabError


def cartan_matrix(A):
    Q = A.quiver
    if Q is None or A.basis_labels is None:
        raise QuiverlabError(
            "Cartan matrix needs the quiver presentation",
            hint="construct the algebra via Quiver.algebra(...); "
                 "structure-constant algebras carry no path basis",
        )
    verts = list(Q.vertices)
    vindex = {v: k for k, v in enumerate(verts)}
    n = len(verts)
    C = [[0] * n for _ in range(n)]
    for label in A.basis_labels:
        if label.startswith("e_"):
            # trivial path at vertex v: label 'e_<v>'; recover v by matching reprs
            v = next(w for w in verts if f"e_{w}" == label)
            C[vindex[v]][vindex[v]] += 1
        else:
            word = tuple(label.split("*"))
            C[vindex[Q.word_source(word)]][vindex[Q.word_target(word)]] += 1
    return C


def coxeter_matrix(A):
    C = sympy.Matrix(cartan_matrix(A))
    if C.det() == 0:
        raise QuiverlabError(
            "Cartan matrix is singular: the Coxeter matrix -C^{-T} C is undefined",
            hint="this happens e.g. for infinite global dimension with |det C| != 1",
        )
    Phi = -C.inv().T * C
    if any(x.q != 1 for x in Phi):  # entries should be integers for det = +-1
        # keep exact rationals if det != +-1; return exact sympy entries as ints when possible
        return [[sympy.nsimplify(Phi[i, j]) for j in range(Phi.cols)] for i in range(Phi.rows)]
    return [[int(Phi[i, j]) for j in range(Phi.cols)] for i in range(Phi.rows)]


def coxeter_polynomial(A):
    C = sympy.Matrix(cartan_matrix(A))
    if C.det() == 0:
        raise QuiverlabError(
            "Cartan matrix is singular: no Coxeter polynomial",
            hint="see coxeter_matrix",
        )
    Phi = -C.inv().T * C
    t = sympy.Symbol("t")
    return sympy.Poly(Phi.charpoly(t).as_expr(), t)
```

`src/quiverlab/invariants/__init__.py`:

```python
from quiverlab.invariants.cartan import cartan_matrix, coxeter_matrix, coxeter_polynomial  # noqa: F401
```

- [ ] **Step 4: Add the Algebra methods**

Append to the `Algebra` class in `src/quiverlab/core/algebra.py`:

```python
    def cartan_matrix(self):
        """Integer Cartan matrix from the quiver presentation (any field)."""
        from quiverlab.invariants.cartan import cartan_matrix
        return cartan_matrix(self)

    def coxeter_matrix(self):
        """Coxeter matrix -C^{-T} C (exact; loud if the Cartan matrix is singular)."""
        from quiverlab.invariants.cartan import coxeter_matrix
        return coxeter_matrix(self)

    def coxeter_polynomial(self):
        """Characteristic polynomial of the Coxeter matrix, as an exact sympy Poly."""
        from quiverlab.invariants.cartan import coxeter_polynomial
        return coxeter_polynomial(self)

    def _require_prime_field(self, what):
        from quiverlab.errors import FieldError
        from quiverlab.fields.primefield import PrimeField
        if not isinstance(self.domain, PrimeField):
            raise FieldError(
                f"{what} is available over GF(p) today (fast engine); "
                f"this algebra is over {self.domain.name}",
                hint="construct the algebra over a prime field, or wait for the "
                     "later phase that generalizes this invariant",
            )

    def cyclic_homology(self, top):
        """Dimensions of HC_0..HC_top (Connes mixed complex; GF(p) via the engine)."""
        self._require_prime_field("cyclic homology")
        from quiverlab.engine.adapter import to_engine
        from quiverlab.engine.cyclic import cyclic_homology_dims
        from quiverlab.hochschild.table import HHTable
        dims = cyclic_homology_dims(to_engine(self.unit_adapted()), top, self.domain.p)
        return HHTable(dims, "HC_", repr(self).splitlines()[0],
                       engine="hanlab engine (F_p fast rank)")

    def nakayama_automorphism(self):
        """Nakayama automorphism matrix (Frobenius algebras over GF(p) via the engine)."""
        self._require_prime_field("the Nakayama automorphism")
        from quiverlab.engine.adapter import to_engine
        from quiverlab.engine.coxeter import nakayama_automorphism
        return nakayama_automorphism(to_engine(self.unit_adapted()), self.domain.p)

    def is_frobenius(self):
        """Is the algebra Frobenius? (GF(p) via the engine.)"""
        self._require_prime_field("the Frobenius test")
        from quiverlab.engine.adapter import to_engine
        from quiverlab.engine.coxeter import is_frobenius
        return bool(is_frobenius(to_engine(self.unit_adapted()), self.domain.p))

    def is_symmetric(self):
        """Is the algebra symmetric? (Frobenius with identity Nakayama automorphism; GF(p).)"""
        self._require_prime_field("the symmetry test")
        if not self.is_frobenius():
            return False
        from quiverlab.engine.adapter import to_engine
        from quiverlab.engine.coxeter import is_identity, nakayama_automorphism
        E = to_engine(self.unit_adapted())
        return bool(is_identity(nakayama_automorphism(E, self.domain.p), self.domain.p))
```

**Signature note (resolve at implementation time):** match the real ported signatures of `cyclic_homology_dims`, `nakayama_automorphism`, `is_frobenius`, `is_identity` (argument order/keywords for the prime; return shapes). Normalize `cyclic_homology_dims`'s return to `list[int]` for degrees `0..top` inside the method. The Step-1 tests define the contract; adjust ONLY the wrapper calls, never the engine.

- [ ] **Step 5: Run the new tests, then the full suite**

Run: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/invariants -q` then `... -m pytest -q`
Expected: all green.

- [ ] **Step 6: Commit**

```bash
git add src/quiverlab/invariants tests/invariants src/quiverlab/core/algebra.py
git commit -m "feat: Cartan/Coxeter invariants (all fields) and engine-backed GF(p) extras"
```

---

### Task 14: Acceptance — deep homology demo, pure-path sweep, docs

**Files:**
- Create: `tests/engine/test_acceptance_plan02.py`
- Modify: `README.md` (Status section), `docs/plans/ROADMAP.md` (mark Plan 02 delivered)

**Interfaces:**
- Consumes: everything above.
- Produces: the Plan-02 exit criteria, executable.

- [ ] **Step 1: Write the acceptance test**

`tests/engine/test_acceptance_plan02.py`:

```python
"""Plan 02 exit criteria: the engine reaches depths the bar oracle cannot,
agrees with the bar oracle where both run, and the whole suite works on the
pure path (QUIVERLAB_NO_NUMBA=1 sweep is run by the task checklist, not here)."""
from quiverlab import GF, Quiver


def test_deep_monomial_homology_beyond_bar_reach():
    # k[x]/(x^2) over GF(2): HH_n = 2 for all n (known closed form).
    # Bardzell/minimal engine reaches degree 30 in seconds; the bar oracle's
    # C_30 would have 2 * 1^30 cells here (trivial), so use a 3-dim example
    # where bar at degree 30 is impossible: kA_3 with rad^2 = 0 (m = 5).
    from quiverlab.engine.adapter import to_engine
    from quiverlab.engine.resolutions_bardzell import BardzellResolution, MonomialPresentation

    # A_3 with both arrows and relation a*b: monomial presentation
    pres = MonomialPresentation(
        vertices=[1, 2, 3],
        arrows=[("a", 1, 2), ("b", 2, 3)],
        relations=[("a", "b")],
    )
    res = BardzellResolution(pres)
    # dims of the resolution terms must exist out to degree 30 without error;
    # rad-square-zero of A_3: the resolution is finite (gl.dim <= 2), so deep
    # terms vanish -- exactly the kind of certified structural fact the engine
    # provides and the bar complex never could at this degree.
    dims = [res.term_dim(n) for n in range(31)]
    assert dims[0] > 0 and dims[-1] == 0


def test_engine_agrees_with_bar_on_truncated_x4():
    # k[x]/(x^4) over GF(5), engine vs pure bar through the public API route.
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    A = Q.algebra(relations=["x^4"], field=GF(5))
    assert (A.hochschild_homology(5, engine="fast").dims
            == A.hochschild_homology(5, engine="bar").dims)
```

**Signature note:** `MonomialPresentation`/`BardzellResolution` argument shapes and the term-dimension accessor name must be matched to the ported code (bank convention per the inventory: `arrows=(id, src, tgt)`, `relations=path-tuples`; the accessor may be `term_dim`, `dims`, or similar). Adjust the test to the real API — the CONTRACT is: build the A_3 rad-square-zero presentation, walk the resolution to degree 30, see finite gl.dim reflected as vanishing terms.

- [ ] **Step 2: Run the acceptance test, the full suite, and the pure-path sweep**

Run, in order:
1. `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest tests/engine/test_acceptance_plan02.py -q` — expected: pass.
2. `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q` — expected: full suite green (Plan-01's 82 + all ported + new; record the final count).
3. `QUIVERLAB_NO_NUMBA=1 NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q` — expected: same count green on the pure-Python kernel path.

- [ ] **Step 3: Update the docs**

README.md Status section: replace the "Foundations phase" paragraph with a Plan-02 status naming the new capabilities (fast GF(p) engine; cup/cap/Gerstenhaber bracket at engine level; cyclic homology, Cartan/Coxeter, Nakayama/Frobenius/symmetric over GF(p); deep monomial resolutions) and stating honestly that classy `A.cup(u, v)` awaits the cohomology-classes machinery of a later phase. ROADMAP.md: mark Plan 02 delivered with the final test count.

- [ ] **Step 4: Commit**

```bash
git add tests/engine/test_acceptance_plan02.py README.md docs/plans/ROADMAP.md
git commit -m "feat: Plan 02 acceptance - deep engine, dual-path sweep, docs"
```

---

## Ported-test inclusion checklist (expected; verify per-file at port time)

Included: test_kernels, test_linalg_fast, test_linalg_modp, test_engine_validation, test_multivertex_engine, test_scan2_builders, test_cohomology_fast_rank, test_complexity (closed subset), test_bimodule_coefficients, test_twisted_bimodule, test_automorphism, test_nakayama_automorphism (closed subset), test_basis_transport (iff closed), test_tt_calculus, test_gerstenhaber, test_cyclic_homology, test_resolution_protocol, test_resolution_contract, test_minimal_resolution, test_minimal_memory_guard, test_periodic_resolution, test_periodic_symmetric_family (closed subset), test_deepen, test_truncated_resolution (iff closed), test_error_paths (iff closed), test_recheck (iff closed), test_bardzell_resolution, test_qpa_reference_validation (iff closed), test_dynkin_zoo, test_theoremB_convention, test_theoremB_multivertex, test_coxeter_spectrum (exact subset).
Excluded (import labdb/open_zoo/resolutions_cs/scan/scan_open/cluster): test_labdb, test_open_zoo_broaden, test_cs_resolution, test_cs_induced_action, test_cluster_runner, test_scan_drivers.

## Plan self-review notes (done at writing time)

- Spec coverage of this phase vs ROADMAP row 02: kernel stack (T2), minimal+guards+deepen (T9), Bardzell (T10), TT calculus (T7), cyclic (T8), bimodule coefficients (T5), Cartan/Coxeter/Nakayama (T6, T11, T13), hanlab tests travel along (every port task + checklist). Deliberately deferred and named: CS resolution (Plan 04), zoo/batch (Plan 06), exact spectral_radius/mahler_measure (Plan 05), classy cup/bracket API (Plan 04/05).
- Port tasks specify exact source paths, exact transformations (including every recon-identified float literal with its integer-exact rewrite), and gate every step with the ported tests plus the full suite; new code (adapter, dispatch, invariants) is given in full, with explicitly-marked signature-resolution points whose correctness is pinned by equality tests against the pure bar oracle rather than guessed signatures.
- Type consistency: `engine=` parameter values, `HHTable(dims, kind, algebra_repr, engine=...)` (matches Plan-01 constructor), `to_engine`/`engine_*_dims` names used identically in Tasks 12–14.
