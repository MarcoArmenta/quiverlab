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
  `minimal_homology_dims(...)`, `minimal_cohomology_dims(...)` (Plan 16 Hom-collapse:
  coh acts `a·w·b` vs homology's `b·w·a`; corner blocks are the SWAPPED tag
  `e_v A e_w`), `hochschild_dimension(...)`, `class AeEngine`.
- **Periodic** — `engine/resolutions_periodic.py`: thin family wrappers
  (`QuantumCIResolution`→CS, `CyclicNakayamaResolution`→Bardzell).
- **Resolution contract / bar** — `engine/resolutions.py`: `Resolution` ABC,
  `BarResolution`, `TruncatedPolynomialResolution` (closed-form `k[x]/(x^a)`).
- **Chouhy–Solotar** — top-level package `src/quiverlab/resolutions_cs/` (runs over
  any exact `Domain`, not just GF(p)):
  `resolution.py::ChouhySolotarResolution` (differential `delta_terms`/`d_terms`/
  `matrix`, `assert_dd_zero`, `assert_order_condition`; differentials are
  **canonical/byte-reproducible** since Plan 17 — the correction solve is reduced
  mod its nullspace via `fields.linalg.reduce_mod_nullspace`), `ambiguities.py::SSequence`,
  `homology.py` (`cs_cohomology_dims`, `cs_homology_dims`, `cs_hh_basis`),
  `comparison.py::Comparison` (CS↔bar class transport), `engine_facade.py::CSResolution`,
  `build.py::reduction_system_of(A)`. Reach it via `engine="cs"` or
  `import quiverlab.resolutions_cs`.
- **Disambiguation:** `src/quiverlab/modules/resolution.py` is a *separate* minimal
  resolution for **right modules** (module Ext, Plan 05) — not HH. Don't conflate.

**Plan 12 (2026-07-22, delivered):** the former non-quadratic non-monomial
`NotImplementedError` is **lifted** — `_require_in_scope` is gone; every admissible
presentation computes, certified per instance (d²=0 + order gate + bar window; only
refusal left is an inconsistent correction solve). Two Plan-12 facts to keep in mind:
(1) block decompositions cut at **first reducibility** — the witness tip may straddle
the block boundary (the old exact-pair condition silently missed chains, e.g. tips
`{xx,yy,xyx}` lost `xyxx`/`xyxyx` and Bardzell HH was wrong from degree 2);
(2) the odd (2-term) differential's first term uses
`MonomialPresentation.right_decomposition` (CS §4 `f_n` even), equal to `u_0` only for
quadratic/palindromic tips. Batteries: `tests/engine/test_bardzell_straddle.py`,
`tests/resolutions_cs/test_battery_straddle.py`.

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

## Status (2026-07-23)

Plans 01–11 delivered and merged to `main`. GUI (Plans 10–11 — Pyodide quiver
canvas + live wait estimates) is live on the docs landing page. Deeper-engine work
is the current focus, driven by the standing queue
`docs/plans/DEEPER-ENGINES-BACKLOG.md` (on "continue": take the topmost unchecked
item). Delivered so far: Plan 12 (straddling ambiguities + right decomposition + CS
scope lift), Plan 13 (minimal A^e engine multi-vertex — corner-typed projective
terms; was silently wrong on any multi-vertex input), Plan 14 (comparison map
rebuilt as a homotopy lift — the closed-form block map was not a chain map for any
tip of length ≥ 3, even monomial; cup/bracket transport now serves every admissible
presentation in the bar window, and cap transport exists,
`Comparison.cap_of_cs_classes`), Plan 15 (`engine/deepen.py` corner-mode
checkpoints — multi-vertex algebras now deepen/resume; the payload gains only the
per-degree corner `tags`, everything else is rebuilt from `(A, prime)`; cross-mode
ckpt_dir reuse refuses loudly), Plan 16 (`minimal_cohomology_dims` — Hom-collapse
HH^• on the same minimal/corner resolution for any f.d. algebra; validated vs bar
coh over 4 primes and degreewise vs CS coh to depth 8 — the second deep oracle),
Plan 17 (CS canonicalization — the correction γ is reduced to the free-variables-zero
coset representative, CS differentials byte-reproducible by construction; the 7
former `xfail(strict=False)` byte pins are strict, gated by an adversarial-solver
test that shifts the solve by a nullspace vector and demands identical bytes),
Plan 18 (standing-zoo diversity — zoo records may carry `vertices`/`arrows`
(multi-vertex reduction systems), five diversity records added (2 straddling-
monomial + 3 multi-vertex incl. `line_abc_cde`), diversity gates in
`test_zoo.py`, batch specs carry quiver data, and `_analyze_open` no longer
unit-adapts multi-vertex algebras — unit-adaptation is a local-only convention
that trips the Plan-13 radical guard), Plan 19 (field generality — the five
engine-backed invariants (`cyclic_homology`, `complexity`, `is_frobenius`,
`nakayama_automorphism`, `is_symmetric`) compute over EVERY exact Domain:
GF(p) keeps the engine byte-unchanged; off GF(p): generic (b,B) mixed complex
(`hochschild/cyclic.py`, no quiver needed), relative-Tor Betti complex
(`invariants/betti.py`, H_n = engine rks[n], GF(p)-parity-gated), socle
criterion + inner-ν symmetry (`invariants/frobenius.py`, self-certifying form,
loud when inconclusive); `_require_prime_field` is gone — the only refusal
left is structure-constants algebras off GF(p) for path-basis-needing
invariants, honestly worded), Plan 09 (the server tier — `webapp/`, backlog
Tier-1 item 7: FastAPI `create_app` + spawn-child resource-capped worker fleet
over a single SQLite/WAL queue; two-tier compute (instant sync / queued jobs) +
email magic-link big-job tier; bilingual EN/ES pages with vendored KaTeX,
`/literature`, feedback+admin; deploy assets (Dockerfile/compose/Caddy/
PROVISIONING). All algebra delegated to `import quiverlab` — no user-code exec,
`quiverlab.engine.*` never imported, base wheel stays lean. Executed
subagent-driven with adversarial critics; whole-branch review's four cross-layer
majors fixed. Branch `plan-09-web`, UNMERGED), Plan 20 (native deep-degree CS
**cup** past the bar 2n+1 window — backlog Tier-2 item 1: a comparison-lifted
diagonal `Δ: P → P ⊗_A P` in `resolutions_cs/diagonal.py` (double-PELT ambient,
Koszul-signed tensor differential `d_p⊗1 + (−1)^p·1⊗d_q`), built degreewise by a
per-degree lift-solve `d^{P⊗P}·Δ_n(σ)=Δ_{n−1}(d_nσ)` — the structural clone of
`_d_general`'s correction solve (`solve`+`reduce_mod_nullspace`, canonical; same
loud `NotImplementedError` scope edge, no contracting homotopy). Sign-free `a·w·b`
cup collapse in `cup.py`; `Comparison.cup_of_cs_classes(engine="auto"/"native"/
"transport")` routes native past-window, transported in-window (byte-unchanged).
Leibniz is the sign arbiter; the transported cup anchors it in-window. Cap = Plan
21 (same Δ); bracket stays transported/window-bounded by design. Branch
`plan-20-native-cs-cup`, UNMERGED), Plan 23 (module-theoretic surface, engine
slice — backlog Tier 1b items 1/2/4: the opposite algebra `A^op` as a
first-class Algebra (reversed quiver + transposed structure constants,
`modules/opposite.py`, cached involution), the duality `D` on modules
(`modules/duality.py::dualize`; `D(P_v^{op})≅injective(A,v)` reconciles the
implicit `I_v=D(Ae_v)`), the corner-transpose `Tr`, and the **AR translates**
`τM=D(Tr M)` / `τ⁻M=Tr(DM)` (`M.tau()`, `M.tau_minus()`, `M.dualize()`,
`M.transpose()`, `Algebra.opposite()`). Injective resolutions +
`M.injective_dimension()=pd_{A^op}(DM)` dual to Plan 05 on the same engine
(`modules/injective.py`). Exact module iso (`modules/hom.py::is_isomorphic`,
invertible-hom certificate; char-0 generic-rank/Noether–Deuring; loud beyond
budget). Oracles: theory (hereditary Coxeter `dim τM=Φ⁻ᵀ dim M`, kA_n tables,
self-injective ⇒ inj.dim∈{0,∞}, Nakayama orbits, `max_v inj.dim S_v=gl.dim`)
AND QPA `-m qpa` (`DTr`/`TrD`+`IsomorphicModules`, `ProjectiveResolution`,
`DualOfModule`, `InjDimensionOfModule`; `modules/qpa_module.py::graded_form`
translates our modules to QPA's row convention) across the zoo incl. Plan-18
`line_abc_cde`. Tier 1b item 3 (no-code module GUI/webapp) SPECIFIED in the
plan doc, NEXT slice. Branch `plan-23-module-surface`, UNMERGED), Plan 24
(**left modules alongside right, right as default** — backlog Tier 1b: the
module surface takes `side="right"|"left"`; right stays the default and
byte-unchanged (right-module repr byte-identical, whole existing suite green).
A left A-module IS a right A^op-module, so the Plan-23 op+D engine makes the
side a presentation-only tag on `Module` (`self.side`, `self.base_algebra`,
`self.with_side`); ALL algorithms read only `(algebra, action)` and are blind
to the side — no math forked. `Algebra.simple/projective/injective/module`
gain `side=` (left routes through `A^op` + re-tag); `D` and `Tr` are now
**side-aware** (exchange left↔right over the SAME base algebra — the classical
contravariant form), so `M.dualize()` of a right A-module is now a LEFT
A-module (representation byte-identical, only the tag/repr move — `injective.py`,
`τ`, `τ⁻` numerically unchanged). `is_isomorphic`/`hom`/`ext` refuse loudly
across sides/algebras (`modules/hom.py::_assert_comparable`, before the
dim-vector fast-paths). Two Plan-23 tests spelling `I_v=D(Ae_v)` via
`A.opposite().projective(v).dualize()` moved to the honest
`A.projective(v, side="left").dualize()`. Oracles: ASS2006 (left/right
duality D, AR translates); QPA `-m qpa` crosschecks left τ/τ⁻/inj.dim/
resolutions by FEEDING QPA THE OPPOSITE ALGEBRA (right-module native).
`tests/modules/test_left_modules.py` + `tests/qpa/test_left_modules_qpa.py`.
Branch `plan-24-left-modules`, UNMERGED).
