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
chapters in `docs/internals/`. GUI (Pyodide, Plans 10–11 + the Plan-26 no-code module panel) in `docs/gui/`;
the server tier (Plan 09) + result cache (Plan 25) in `webapp/`.

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
  `comparison.py::Comparison` (CS↔bar class transport; `cup_of_cs_classes` /
  `cap_of_cs_classes` route native past-window (Plans 20/21), transported in-window),
  `diagonal.py` (lifted diagonal Δ), `cup.py` (native `a·w·b` cup, Plan 20),
  `cap.py` (native `b·w·a` cap, Plan 21), `engine_facade.py::CSResolution`,
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
majors fixed. Branch `plan-09-web`, merged 2026-07-25), Plan 20 (native deep-degree CS
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
`plan-20-native-cs-cup`, merged 2026-07-25), Plan 22 (verification transparency, backlog
Tier 1a — `docs/verification.md` documents HOW everything is tested: the two
oracle classes (theory/literature pins on constructed examples; cross-engine +
QPA/GAP agreement wherever QPA implements the feature), a subsystem→oracles→tests
table, the marker/bucket scheme with audited counts (1377 tests: fast 504 / deep
868 / qpa 5), the CI matrix, and an honest-scope section naming where QPA cannot
compare and the theory oracle that covers each gap; wired into the mkdocs nav +
a README section; standing rule: every future plan adds its new oracles to the
verification page as part of acceptance. Docs+audit only, no `src/` change.
Branch `plan-22-verification-transparency`, merged 2026-07-25), Plan 21 (native deep-degree CS **cap** past the
bar window — backlog Tier-2: the homology-side sign-free `b·w·a` collapse of the SAME
Plan-20 diagonal Δ (`resolutions_cs/cap.py::native_cap`,
`(f∩z)(ρ)=Σ_Δ coeff·b_c·x·b_a·f(τ)·b_mid`; f eats the degree-p first factor τ, ρ
survives — the covariance flip of the Plan-16 note, cohomology `a·w·b` vs homology
`b·w·a` on swapped corner tags). `Comparison.cap_of_cs_classes(engine="auto"/"native"/
"transport")` routes native past-window, transported in-window (byte-unchanged). The
sign convention is ARBITRATED, not assumed — the in-window native ≡ transported cap
anchor (non-commutative quantum CI distinguishes `b·w·a` from `a·w·b`), exact unit cap
`1∩z=z`, exact cap-Leibniz `b(f∩z)=(-1)^{p+1}(δf∩z)+(-1)^p(f∩bz)`, and module identity
`(z∩f)∩g ~ z∩(f∪g)` via the native cup all hold simultaneously. Edges: n=p→C_0, p>n
raises. Folded-in Plan-20 review backlog cleared (QCI skip-guards now FAIL not skip,
element-wise CS-basis bridge, session-scoped Δ_4 fixture, `native_cup` chain-cache
routing, constructed inconsistent-lift refusal). Branch `plan-21-native-cs-cap`,
merged 2026-07-25 — cap oracles added to `docs/verification.md` at merge), Plan 23 (module-theoretic surface, engine
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
plan doc, NEXT slice. Branch `plan-23-module-surface`, merged 2026-07-25), Plan 24
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
Branch `plan-24-left-modules`, merged 2026-07-25). Plan 25 (webapp result cache — backlog Tier-3
item: never recompute a known example. `webapp/server/cache.py` canonicalizer
(`canonical_key` = sha256 of sorted-keys JSON over the versioned request + library
version; dict order irrelevant, tuple/list round-trip collides, version bump
invalidates) + a new `result_cache` table in the single SQLite/WAL store — NOT a
jobs column — that maps the canonical key → a finished job and "pins" that job
against the retention purge so its artifacts survive to back replays. All three
tiers check the cache FIRST (`/api/compute`, `/api/jobs`, `/api/jobs/big`); a hit
replays instantly across users with zero recompute. Big-job crux: the cache check
is the first statement in `submit_big`, BEFORE `big_jobs_enabled`, so a cached big
example is served with NO token/email/`email_hash` — email verification gates the
COST of computing, not access to the mathematics — even when SMTP is off. Worker
records on success only (failures may be transient); `cache_put` is idempotent
(benign identical-request race, in-flight dedup skipped as it needs a jobs-schema
change); `sweep_cache_once` (version purge + LRU size cap) runs beside `sweep_once`.
Rows are math-only (no email/ip/token). Config `QLWEB_CACHE_ENABLED`/
`QLWEB_CACHE_MAX_ENTRIES`. Branch `plan-25-webapp-result-cache`, merged 2026-07-25),
Plan 26 (no-code module input GUI + webapp — backlog Tier 1b item 3, the last:
every representation theorist specifies a module with zero code and reads off the
classical module invariants. Webapp request **schema v2** gains a `module` block
(explicit `{dims, maps, side}` — one exact-entry matrix per arrow — OR a
zero-typing `{builtin: {kind: simple|projective|injective, vertex, side}}`) plus
`ext_target` (the N in Ext), guarded to schema 2, canonicalizing through the
Plan-25 `canonical_key` (side default explicit so an omitted `side` and `"right"`
collide; a `ComputeRequest.model_dump` override drops an absent module block so
every existing family/quiver request keeps its byte-identical key). Module compute
kinds — dimension_vector, rad_top_soc, ext, tau/tau_minus, projective/injective
_resolution, projective/injective_dimension — served by ALL THREE tiers with
references/citations (GSZ2001 `minimal_resolution`/`module_ext`, ASS2006
`assem_book`); a relation-violating module raises the library's loud error as a
clean typed 4xx (never a 500). Matrix entries are exact DATA (int/`"1/2"` strings,
never eval, floats refused); per-arrow BLOCK matrices `dim[t]×dim[s]` expand to
the full vertex-ordered action `A.module` consumes, using the representation
quiver's directions — right (A) and left (A^op) alike. `estimator.sizing_dim`
sizes on the larger of algebra/module dim so big modules route off instant like
oversized families. Two runners carry the SAME dispatch: `webapp/server/runner.py`
(server) and `docs/gui/runner.py` (Pyodide client). GUI: a "Module (no code)"
panel on the canvas — per-vertex dimension picker, per-arrow matrix grid (dims
follow source/target live), a right/left side toggle, S(v)/P(v)/I(v) pick-lists;
results render like existing blocks (MathJax + citations); interface-freshness
pins the module surface. Branch `plan-26-no-code-modules`, merged 2026-07-25), Plan 27 (Yoneda Ext-algebra + Koszulity, Tier 2 — the Ext quiver equals Q with
left-to-right products so E(A) ≅ (A^!)^op, arbitrated by worked anchors + the
QPA ext-quiver pin; Ext dims = graded Betti numbers off Plan-05 minimal
resolutions (minimality ⇒ δ≡0, cocycles = summand projections), Yoneda product
by canonical chain-map lifting (reduce_mod_nullspace, byte-reproducible),
degreewise minimal generators/relations over k^{Q_0};
`Algebra.ext_algebra(top)` → `modules/ext_algebra.py::YonedaPresentation`
(honest certification; `as_algebra()` iff gl.dim exact-finite); Koszulity
three-valued in `modules/koszul.py` — G-quadratic certifier (Priddy PBW),
quadratic dual on Q^op, Fröberg falsifier; citations priddy/froberg_koszul/
polishchuk_positselski; 43 live QPA tests (`ExtAlgebraGenerators`,
`ExtOverAlgebra`, `IsQuadraticIdeal`; QPA has NO IsKoszul — honest scope on the
verification page). Branch `plan-27-ext-algebra`, merged 2026-07-25), Plan 28
(container HPC batch tier + offline laptop app, Tier 1c — the `quiverlab-hpc`
CLI lives in the WHEEL (`src/quiverlab/hpc/`: stdlib-validated spec core
promoted from the webapp runner, which now delegates with byte-stable results
AND Plan-25 cache keys pinned by frozen goldens; verbs run/render/sample-config/
estimate/gui/version/selftest; `[hpc]` extra = pyyaml; exit 75 = clean
checkpoint stop via the finally-wired `engine/deepen.py`, sbatch templates
requeue on it); container/ Dockerfile+quiverlab.def (python:3.12-slim by
digest, tectonic pre-warmed for offline PDF via the extended `quiverlab.trace`
ladder PDF→HTML→txt, result_schema envelope), GHCR OCI primary + SIF release
asset via .github/workflows/container.yml; slurm/ templates auto-detect
apptainer→singularity→venv (drac-local emulator = the venv path);
`quiverlab-hpc gui` serves the Plan-09 webapp fully OFFLINE on localhost
(embedded worker, big-jobs/SMTP off, seeded Plan-25 result cache from
webapp/precomputed/manifest.yaml — 5 placeholder examples, Marco to curate —
WAL-checkpointed at build), estimator gains exact memory estimates (bilingual),
GUI + job pages export the cluster config YAML; `quiverlab/hpc/resources.py`
detects cores/RAM (cgroup+SLURM-aware) and reports GPUs as detected-but-UNUSED
(exact CPU engines — sbatch/docs say request cores+RAM, never GPUs). Branch
`plan-28-hpc-container`, merged 2026-07-25). Also merged 2026-07-25: the
literature-oracle deep research (4 author-cluster reports in
`docs/plans/2026-07-25-literature-oracles-deep-research.md` + a Tier-1a
battery backlog item headlined by a QPA-verified `is_symmetric` bug on
multi-vertex symmetric Nakayama and TrivialExtension — DO FIRST), and the two
METAGOALS (README top): (1) no code required, (2) any computation in
representation theory — with the C1–C8 coverage program in
`docs/plans/ROADMAP.md` (C1 categorical glue: morphisms + Krull–Schmidt; C4
τ-tilting engine = τ-tilting/torsion/stability/silting in one, white space
even in QPA) grounded in `2026-07-25-metagoal-coverage-deep-research.md`.
 Plan 29 (literature-oracle batteries + fixes, Tier 1a — 79 new oracle tests
across invariants/engine/resolutions_cs/families with ZERO
literature-vs-engine mismatches except the DOCUMENTED CRS-2004 Example-2.20
discrepancy (bar oracle gives dim HH¹=1, not the paper's 0 — the verified
value is pinned, honest-scope entry on the verification page); the
`is_symmetric` bug FIXED (the GF(p) ν==identity branch was
sufficient-not-necessary — silent False on Brauer stars kZ_n/J^L with
n|(L−1); replaced by the Skowroński–Yamagata nondegenerate-trace-form
certifier, loud on small fields) + `is_weakly_symmetric` + QPA symmetry
crosschecks; **module Tor**: `modules/tor.py::tor_dims(A, M, N, top)` right⊗left
with `Algebra.tor`, anchored by `dim Tor_n(M,N) = dim Ext^n(M, DN)` on every
test + resolve-either-side balance (QPA has NO native Tor — probed live; the
QPA bridge computes Ext-via-NthSyzygy dimension shifting); 17 citation keys
added; TrivialExtension symmetry now refuses LOUDLY pending its double-quiver
presentation (backlogged, xfail fences auto-flip — delivered by Plan 31). Branch
`plan-29-literature-batteries`, merged 2026-07-25), Plan 30 (Krull–Schmidt +
GUI/report completeness, Tier 1d — Marco's feedback: `modules/decompose.py`
`decompose`/`is_indecomposable` (exact Fitting splitting on Domain-native
minimal polynomials; Dickson/CIW trace-rank locality certificate rigorous for
char 0 or char > dim M; LOUD char-p refusals; QPA
`DecomposeModuleWithMultiplicities` crosschecks); webapp/GUI: `tor` +
`decompose` compute kinds (schema v2, cache keys byte-stable), the
second-argument no-code matrix editor for Ext/Tor targets, τ/τ⁻ blocks
certify indecomposability (honest omission when uncertifiable), resolution
tables render term | `P_1^{2} ⊕ P_3` (dropped #summands + dim-vector columns
per Marco); the exhaustive worked-steps bundle: module trace events with
rendered differential matrices (stated elision past the size threshold), the
"projectives and injectives of A" Loewy section, `trace.tex` downloadable
everywhere (webapp whitelist + job page, `quiverlab-hpc render --format tex`,
GUI button), kA₂ replay golden (every differential of the S₁ resolution
verbatim). Branch `plan-30-decompose-reporting`, merged 2026-07-25 after
Plan 29 — the skip-guarded tor wiring auto-activated in the merged tree).
Plan 31 (TrivialExtension double-quiver presentation, Tier 1a —
`TrivialExtension(A)` of a presented `A` over QQ/GF(p) now returns a genuine
`kQ_T/I_T`-presented Algebra (via `Quiver.algebra`): quiver = Q_A plus one
arrow dual to each corner-homogeneous basis element of the bimodule socle
`soc_{A^e}(A)` (direction reversed), relations extracted from the ⋉ structure
by a length-lex kernel enumeration, certified per instance by
`dim kQ_T/I_T = 2·dim A` (loud `QuiverlabError` otherwise); presentation-less
bases keep the unchanged ⋉ structure-constants build as
`_trivial_extension_structure_constants` (honest refusals preserved, and that
build doubles as the internal iso-invariance oracle). `is_symmetric`/
`is_weakly_symmetric`/`is_frobenius`/`is_selfinjective` are now True on every
`T(A)` via the unchanged Plan-29 trace-form certifier, and the four
`xfail(strict=False)` fences in `tests/invariants/test_symmetric_regression.py`
are real asserts. QPA 1.37's native `TrivialExtensionOfQuiverAlgebra` is a
construction oracle (dim + arrow count + `IsSymmetricAlgebra`/
`IsWeaklySymmetricAlgebra`/`IsSelfinjectiveAlgebra`). Oracles: `T(kA_n) ≅
kZ_n/J^{n+1}` (n=2,3,4), `T(k[x]/(x^a)) = k⟨x,y⟩/(x^a,y²,xy−yx)` (plain
commutator in every characteristic; a=2 reproduces the `k[x,y]/(x²,y²)`
`HH_• = [4,4,5,6]` pin), Cartan identity `C_T = C_A + C_Aᵀ`, presented ≡ ⋉
iso-invariance (bar-HH degreewise), and CS ≡ bar on presented `T` (CS refused
the old structure-constant build). No Fernández–Platzeck citation by design
(metadata not BibTeX-verifiable; the build is per-instance-certified +
QPA-oracled). New tests: `tests/families/test_trivial_extension_presented.py`
(deep), `tests/qpa/test_trivial_extension_qpa.py` (qpa). Branch
`plan-31-trivial-extension`).
Plan 32 (oracle-class markers, 2026-07-26, branch `plan-32-oracle-markers`,
Marco's four-part test taxonomy as ORTHOGONAL pytest markers — v0.1.0 release
gate): every test classifiable and each class runnable standalone —
`oracle_literature` (670: frozen literature/theory values), `oracle_crossengine`
(396: two independent implementations agree), `oracle_selfcert` (604: internal
certificates — d∘d=0, order, dimension/canonicality), the existing `qpa` marker
IS the fourth class (112, never double-marked); unmarked = contract &
infrastructure; overlap allowed (union 1227). Audited: the verification-page
class table == live collection, gated by `tests/release/test_oracle_classes.py`
(3 deep tests — the only count delta; the sweep itself is bucket-byte-identical,
verified 804/1180/112 excluding the gate). 89 files marked; 12 edge rulings in
the plan doc (`2026-07-26-plan-32-oracle-markers.md`).
Plan 34 (homework-grade reports + honest artifacts, 2026-07-26, branch
`plan-34-homework-reports`, /spawn /devils: 4 workers + 4 paired adversarial
critics + 3-agent fix round, RELEASE-GATING — Marco: "I don't think it is
ready"): rad/top/soc returned/rendered as FULL representations ({dims, maps}
in the Plan-26 input schema via the shared `module_blocks` serializer, dim
column dropped, GF(p^n) entries domain-notated + display_only); Ext/Tor
batteries with certified interior modules (kA5 intervals — kA4 has exactly
one; Euler-form asymmetric arbiter, Tor≡Ext∘D + balance + live QPA anchors,
symmetric ext_target schema guard); the artifact contract PDF+HTML+JSON all
cached (HTML = complete record, no elision, scrollable; trace.json = exact
versioned event stream, served/promoted/downloadable everywhere; PDF = the
page-bounded homework document — event-computed MaxMatrixCols preamble +
shrink-only resizebox 11–25 cols, stated elision >25 pointing at HTML/JSON;
recorder keeps FULL matrices, 250k-cell memory backstop only); homework-depth
content (definitions + matrices + pivot/rank lines + justifications for
rad/top/soc/resolutions/Ext/Tor/tau/tau-minus/decompose, side-aware left
narration, loud drift gates vs engine dims, correct rad=MJ attribution);
renderers REFUSE foreign stream objects loudly (ALL_EVENTS gate), Result
footers segmented per computation, self-map labels carry the module's name;
worker-recorded why-no-PDF reasons; MathML converter with a real
escaped-source fallback + leak-guard test. PLUS the Marco-2026-07-26 DISPATCH
AMENDMENT: engine="auto" no longer dies at the bar/fast depth wall — it falls
back to Chouhy-Solotar for quiver-presented algebras exactly where
DepthLimitError fired (recorded in the dispatch trace, in-window results
byte-unchanged, explicit engines keep honest walls, presentation-less refuses;
the stale "later phases" bar hint rewritten; Pillar-4 pin updated;
`tests/hochschild/test_auto_cs_fallback_p34.py`).
Suite recounted post-Plan-34: 2295 tests (fast 930 / deep 1243 / qpa 122).
Plan 33 (nontrivial literature examples at scale, Tier 1a — Marco: the
test/paper examples were too small (kA₂/kA₃/single loops); scale them to the
books/literature and push the quantum ones deeper. SCALE, not new theorems —
Plan 29 already pins the small directions. Scale batteries: generalized quantum
CI `k⟨x,y⟩/(x^a,y^b,yx−q·xy)` for `(a,b)∈{(2,4),(3,4),(4,4),(2,5),(5,5)}`
(Bergh–Erdmann coh `[2,2,1,0,…]` past deg 8, hom `[a+b−1,…]`; char-0 branch
only, `qci_hh_oracle`); preprojective Π(A₄/A₅/D₄/D₅) dims 20/35/28/60,
self-injective, Loewy = h−1 (Coxeter 5/6/6/8) — HH values xeng/QPA-only;
Bardzell depth kZ₂₀/J¹¹ dim 220 to degree 300 (context-managed
`sys.setrecursionlimit` raise around the Bardzell walk in
`engine/resolutions_bardzell.py`, pure/numba parity unaffected; deg-300
regression); symmetric Brauer stars kZ₄/J⁹, kZ₅/J¹¹; Taft Λ₅/Λ₆ HH + cyclic;
canonical C(2,2,2,2,2) `HH²=t−3=2` (first ≥2 case); Boolean B₃ incidence dim 27
(nerve vanishing `HH^{≥1}=0`); presented T(kD₄)/T(kA₅)/T(kA₆); wild m-Kronecker
`HH¹=m²−1`, Coxeter `t²−(m²−2)t+1`; exterior Λ(k³)/Λ(k⁴) Koszul via
`modules.koszul.g_quadratic_certificate`. Builders (C2 src): `QuantumCI(q,a,b)`
(byte-identical `QuantumCI(q)`), preprojective auto degree-bound table, the
Bardzell recursion guard. Papers (C3): JSC worked-examples rebuilt around the
research top-10 with every number recomputed by replayable
`paper-jsc/computations/` scripts + the representation-theory-first interior
pass; JOSS Research-impact folds the QCI-(a,b), m-Kronecker, and Bardzell-depth
one-liners in band. Honest-scope labels binding (verification page):
preprojective and exterior HH values are cross-engine-only (no published table);
Π(D₅) HH-at-scale, Λ(kⁿ≥4) depth, decompose ≳ dim 50, dim-30+ non-monomial HH
past ~degree 10, and the Π(E₆)/Π(D₆) builds are deferred to the SUBMISSION
step-4 cluster list. Oracle-class markers per Plan 32. Branch
`plan-33-nontrivial-examples`.)

**Marco's report-completeness pass (2026-07-29, branch
`marco-report-completeness`)** — feedback on the containerized macOS desktop app
(his `exmple-a.html` = the draw page after Compute, `example-b.html` = the
downloaded worked-steps report). Twelve items, all delivered:

*Presentation, both GUI renderers (`docs/gui/gui.js`, vendored byte-identical to
`webapp/static/gui/gui.js`, and `webapp/static/app.js`)* — the Plan-34 scroll box
is GONE (`mathScroll` → `mathFit` + a post-typeset `fitMath` shrink-to-fit pass;
matrices show COMPLETE, an over-wide one scales down, the page body still never
scrolls sideways); an arrow acting as the EXACT zero map is named in one line
instead of printing a zero block (`matIsZero`); τ/τ⁻ render the translate's FULL
per-arrow matrices (`block.repr`) and, when the request names a second module N,
N's translate too (`block.targets`); projective/injective resolutions render
their differentials, and a differential equal to an earlier one prints
`d_3 = d_2` instead of repeating the matrix.

*Block data (`quiverlab.hpc.spec` + its Pyodide twin `docs/gui/runner.py`, kept
shape-identical)* — `projective_dimension`/`injective_dimension` gained `latex`
(their absence is what typeset the two literal "undefined"s), and an UNRESOLVED
probe now states the certified lower bound `\operatorname{pd} M > 32` + a `note`,
never a bare `\infty` the engine did not prove; `tau`/`tau_minus` gained
`targets` (the Ext/Tor argument's translate, honest error entry on a loud
refusal; omitted entirely when no target is named, so single-module blocks stay
byte-identical). One frozen golden re-frozen (`module_left_a2`), documented in
`tests/webapp/test_runner_delegation.py`.

*The report (`quiverlab.trace`)* — new `results_html.py` renders EVERY computed
result block into a "Computed results" section (`render_html(..., results=)`,
threaded through `writer.write_trace` and `spec.run`), so the saved HTML is the
whole session, not just one traced computation; a request with NO traceable
computation now still writes the bundle (`_write_results_only`); `.ql-eq`'s
`overflow-x:auto` is gone (wide equations shrink via the integer `_fit_pct`
rule); the Result is a degree TABLE (`_dims_table`) like the GUI's Ext/Tor
tables; repeated differentials are referenced via `_MatrixEcho` (an ELIDED
matrix is never matched — its body was not recorded); and `ResolutionTerm` gained
an optional `corners` field (the CS generators' `(o, t)` pairs) so the worked
resolution steps NAME the term, `C_n = ⊕ P(v,w)`, `P(v,w) = A e_v ⊗ e_w A` — a
term without recorded corners (bar over structure constants) claims nothing.

New EN/ES i18n keys for every string this added to the webapp module blocks.
Tests: `tests/trace/test_report_completeness_m0729.py` (14, `oracle_selfcert`),
`tests/webapp/test_module_blocks_m0729.py` (13, cross-runner contract, unmarked
per the Plan-32 extras-gated ruling), plus additions to
`tests/webapp/test_display_only_p34.py`; `tests/trace/_result_table.py` reads the
Result table back out of rendered HTML for the renderer tests. Verification page
updated (new self-cert oracles + audited counts 633 / 1322, suite 2406).

SECOND PASS (same day, Marco on the regenerated report): a new **"The modules"**
report section describes every module the computation was about -- `M`, and `N`
when a second module was named -- with its Loewy series, top/socle and the exact
matrix of every arrow (`trace/modules.py::module_description`, threaded as
`render_html(..., modules=[(label, Module)])` from both runners; each module is
guarded individually so an undescribable one is skipped, never fatal). Krull-Schmidt
summands go through the new shared serializer
`modules/qpa_module.py::summand_blocks`: a summand isomorphic to a STANDARD
indecomposable is NAMED `S_v`/`P_v`/`I_v` with its matrices omitted (via the new
`modules/hom.py::identify_standard` -- dimension-vector prefilter, then the exact
`is_isomorphic` certificate; undecidable leaves it unnamed and shown in full,
never guessed), every other summand carries its full action; both GUI renderers
mirror the naming. Headings now say WHAT they hold: the bottom "Result" section is
`Hochschild homology`/`Hochschild cohomology` (`_hh_heading`) and is SKIPPED
entirely when the Computed results already carry that table (`_results_carry_hh`),
so the same numbers are never printed twice; the module-steps Ext/Tor footer is
headed `Ext`/`Tor`.

THIRD PASS + CI FIX (same day): every displayed matrix is now an INDEXED GRID --
`render_html.matrix_grid` / `_event_grid` (report) and `matrixGrid` (both GUI
renderers): an HTML table with a header row of column indices, a header column of
row indices and a light-grey rule, 1-based, entries verbatim + escaped; a
zero-dimensional matrix stays the symbol `0`. `pmatrix` remains only as the TeX
SOURCE form. Tests read matrices back out of the page via the new
`tests/trace/_matrix_grid.py` helper, so they assert ENTRIES, not a presentation.
CI FIX (a REAL pre-existing bug, not caused by this work): the whole **Windows**
fast matrix had been failing on `UnicodeDecodeError` because `Path.write_text`
defaults to the LOCALE codec -- cp1252 on Windows -- so the report's em dashes were
written as single cp1252 bytes and every utf-8 reader blew up (invisible on
macOS/Linux, whose locale codec IS utf-8). Every artifact write in
`trace/writer.py`, `hpc/spec.py`, `webapp/worker/worker.py`, plus the
locale-dependent `open()`s in `engine/deepen.py` and `engine/scan2.py`, now pass
`encoding="utf-8"` explicitly, gated by `tests/trace/test_artifact_encoding.py`
(an AST scan of the shipping tree + a live round-trip under a forced cp1252
locale).

KNOWN GAP (not in Marco's list, flagged not fixed): the INSTANT tier deletes its
artifact dir by design (`webapp/server/instant.py`), so a computation served
instantly produces no report at all. Only queued/cached jobs expose
`trace_steps.html`.
