# quiverlab — project memory

Exact computation with finite-dimensional algebras `kQ/I` over CC / GF(p) / GF(p^n):
certified finiteness, Hochschild (co)homology, TT-calculus, resolutions (bar /
minimal / Bardzell / Chouhy–Solotar), modules & invariants. **Exact only — floats
fail loudly by design.**

## Commands (always from repo root)

- Python is **always** `/Users/marco/Desktop/HomologicalNetworks/quiverlab/.venv/bin/python`.
  System python is 3.8 — never use it.
- Tests: `NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python -m pytest -q -m <marker>`
  - Markers (`pyproject.toml`): `fast` (CI OS×Py matrix), `deep` (heavy
    engine/resolution/CS suites, ~19 min, one Linux cell), `slow` (opt-in;
    implies `deep`), `qpa` (needs `[qpa]` extra).
  - Test buckets are **auto-assigned by directory** in `tests/conftest.py`:
    `tests/{engine,resolutions_cs,modules,families,batch}/` → **deep**; everything
    else → fast. So deeper-engine tests run under `-m deep`, not `-m fast`.
- Pure-Python kernel path: `QUIVERLAB_NO_NUMBA=1` (or numba absent). Engines must
  agree **exactly** on both the numba and pure paths — parity is gated.
- Strict docs build: `.venv/bin/mkdocs build --strict` (exit 0). Executes tutorial
  notebooks and pip-builds the wheel — takes several minutes.
- Extras: `[fast]`=numba, `[qpa]`, `[docs]`, `[dev]`.

## Hard conventions

- **No floats in `src/`.** AST gate `tests/test_no_floats.py` scans `src/` only.
  All algebra is exact (sympy / int mod p). Non-algebraic glue (`docs/gui/`,
  `webapp/`) is exempt.
- **Path composition is left-to-right:** `a*b` = first `a` then `b`, requiring
  `target(a) = source(b)` (Assem–Simson–Skowroński). Internal currency is the
  unit-adapted structure-constant `Algebra`.
- **Engine internals stay internal.** Public surface is `import quiverlab`; never
  reach into `quiverlab.engine.*` from app/GUI code.
- **Read-only bank** `HomologicalAlgebra/HansConjecture/` is the port source — copy
  with attribution, never modify.
- Conventional commits; green tests at every commit. Commit/push only when asked;
  branch first if on the default branch (`main`).

## Layout

`src/quiverlab/`: `core/` (Quiver, Algebra, dispatch), `fields/`, `groebner/`,
`hochschild/`, `engine/` (fast GF(p) stack + resolutions), `resolutions_cs/` (CS,
top-level), `modules/`, `invariants/`, `families/`, `batch/`, `citations/`,
`trace/`, `viz/`. Plans+roadmap in `docs/plans/`; algebraist "under the hood"
chapters in `docs/internals/`. GUI (Pyodide, Plans 10–11) in `docs/gui/`; Plan 09 is the *planned* server tier (not built).

## Deeper engines — the current focus

Public HH dispatch is `src/quiverlab/core/algebra.py`:
`Algebra.hochschild_cohomology(top, max_cells=4_000_000, engine="auto", auto_cs=False, ...)`
and `hochschild_homology(...)`. **Valid `engine=` is exactly `{"auto","bar","fast","cs"}`.**
`auto` routes to CS only when `auto_cs=True` and the tip is non-monomial admissible;
`fast`/auto-over-GF(p) uses the **bar basis** (accelerates rank only) via
`engine/adapter.py` → `scan3`/`hh_engine`.

**Bardzell / minimal / periodic are engine-internal** — NOT reachable through the
public `engine=`; they are selected via the `resolution=` kwarg of
`engine.hh_engine.hochschild_homology_dims` (used by tests, the periodic wrappers,
and CS). Over int64 GF(p).

- **Bardzell** (monomial `kQ/I`) — `engine/resolutions_bardzell.py`:
  `MonomialPresentation` (factories `truncated_polynomial(a)`, `cyclic_nakayama(n,ell)`,
  `local_radsq(g)`; combinatorics `associated_paths`, `left_decomposition` = Bardzell/Anick
  chains) and `class BardzellResolution(Resolution)`.
- **Minimal `A^e` resolution** (any f.d. algebra over GF(p), iterated syzygies) —
  `engine/resolutions_minimal.py`: `minimal_resolution(A,N,p,...)`,
  `minimal_homology_dims(...)`, `hochschild_dimension(...)`, `class AeEngine`.
- **Periodic** — `engine/resolutions_periodic.py`: thin family wrappers
  (`QuantumCIResolution`→CS, `CyclicNakayamaResolution`→Bardzell).
- **Resolution contract / bar** — `engine/resolutions.py`: `Resolution` ABC,
  `BarResolution`, `TruncatedPolynomialResolution` (closed-form `k[x]/(x^a)`).
- **Chouhy–Solotar** — top-level package `src/quiverlab/resolutions_cs/` (runs over
  any exact `Domain`, not just GF(p)):
  `resolution.py::ChouhySolotarResolution` (differential `delta_terms`/`d_terms`/
  `matrix`, `assert_dd_zero`, `assert_order_condition`), `ambiguities.py::SSequence`,
  `homology.py` (`cs_cohomology_dims`, `cs_homology_dims`, `cs_hh_basis`),
  `comparison.py::Comparison` (CS↔bar class transport), `engine_facade.py::CSResolution`,
  `build.py::reduction_system_of(A)`. Reach it via `engine="cs"` or
  `import quiverlab.resolutions_cs`.
- **Disambiguation:** `src/quiverlab/modules/resolution.py` is a *separate* minimal
  resolution for **right modules** (module Ext, Plan 05) — not HH. Don't conflate.

**Known deeper-engine stretch item (next CS task):** the CS differential raises
`NotImplementedError` for **non-quadratic non-monomial tips** (tip length ≥ 3 with
nonzero tail, `v_n ≠ u_0`) — see `resolutions_cs/resolution.py` `_require_in_scope`/
`_d_general`. Lifting it needs a `right_decomposition` upgrade (right block
factorization so the even `f_n` first term is exact for all tips).

**Tests & oracles:** `tests/engine/` (Bardzell/minimal vs normalized-bar oracle,
deep-depth past the bar blow-up, pure/numba parity, memory guards) and
`tests/resolutions_cs/` (validation battery: `test_battery_bar.py` CS≡bar degreewise,
`test_battery_bardzell.py` CS≡Bardzell on monomial, `test_battery_bank_oracle.py`
byte-level bank closed-form, `test_battery_literature.py`, plus `d∘d=0`/order,
comparison, homology, dispatch). Oracles: normalized bar complex, primes
`(32003, 2, 3, 5)`, the hanlab bank closed-form, literature values.

**Design docs:** `docs/plans/2026-07-18-plan-02-hanlab-port.md` (minimal/Bardzell),
`docs/plans/2026-07-18-plan-04-chouhy-solotar.md` (CS), `docs/plans/ROADMAP.md`.
Internals: `docs/internals/05-resolutions.md`, `docs/internals/09-chouhy-solotar.md`.

## Status (2026-07-22)

Plans 01–11 delivered and merged to `main`. GUI (Plans 10–11 — Pyodide quiver
canvas + live wait estimates) is live on the docs landing page. Deeper-engine work
is the current focus.
