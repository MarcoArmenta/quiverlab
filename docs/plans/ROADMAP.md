# quiverlab v1 — Implementation Roadmap

Spec: `docs/specs/2026-07-12-quiverlab-design.md` (approved 2026-07-12).
Each plan below produces working, testable software on its own and ends with all
tests green. Plans are written when their phase starts (not speculatively); this
roadmap fixes scope and interfaces between them.

| # | Plan | Delivers | Spec sections |
|---|------|----------|---------------|
| 01 | **Foundations** (`2026-07-12-plan-01-foundations.md`, written) | Package skeleton; exact fields CC/GF(p)/GF(p^n) + generic exact linear algebra; Quiver + relation parser; structure-constant `Algebra` (unit-adapted); monomial `kQ/I` with certified finiteness; normalized bar complex → HH^n/HH_n dims; two starter builders; float-ban AST gate | §3.2, §3.3 (monomial), §4, §5 (components 1, 2, 4; bar of 5), §7 |
| 02 | **hanlab engine port** (`2026-07-12-plan-02-hanlab-port.md`, DELIVERED — 591 tests green on both the numba and pure-Python `QUIVERLAB_NO_NUMBA=1` paths) | Fast GF(p) kernel stack (pure/sparse/numba, equality-gated) behind the Domain interface; `minimal` A^e resolution + guards + checkpointed deepen; `bardzell`; Tamarkin–Tsygan calculus (cup/cap/bracket) + cyclic homology + bimodule coefficients; Cartan/Coxeter/Nakayama; hanlab test suite travels along. Source is read from the read-only bank `HomologicalAlgebra/HansConjecture/` (copy with attribution; never modify the bank). Deferred and named: CS resolution + memory-guard tests (Plan 04); zoo/labdb/periodic-symmetric-family (Plan 06); classy `A.cup`/bracket API and exact spectral invariants (Plan 04/05) | §5 (1, 5, 6, 8 partial), §8 ring 1 |
| 03 | Gröbner + general kQ/I | Noncommutative overlap completion, degree bound, admissibility certificate; general relations lower to `Algebra`; reduction systems emitted for CS | §5 (component 3), §3.3 (general) |
| 04 | Chouhy–Solotar + operation transport | The general CS resolution; comparison morphisms/homotopy liftings; CS validation battery (CS≡Bardzell on monomial; CS≡bar/minimal degreewise; literature oracle) | §6, §8 battery |
| 05 | Modules + invariants | `module_ext` generalized (any vertex set, any domain); simples/projectives/injectives, Hom/End, Ext^n; gl.dim, remaining invariants; `sweep` | §5 (7, 8), §3.5–§3.6, §3.9 |
| 06 | Families + batch | Full §3.4 catalog, `zoo()`, `families()`; `quiverlab.batch` (labdb lift) | §5 (9), §3.4 |
| 07 | Viz + trace | `draw()` (matplotlib hard dep added here) + `tikz()`; the worked-steps trace subsystem (PDF/HTML/text renderers, `verbose=True` default, eliding rules, golden-file tests) | §5 (10, 11), §3.7, §3.8, D9 |
| 08 | QPA extra + release | `[qpa]` extra (passagemath-gap), crosscheck oracle + CI job; GitHub Actions matrix; docs site + CI-run tutorials; PyPI packaging; JOSS paper draft | §5 (12), §8 ring 3, §9–§11 |

**Deeper-engine plans (post-v1, 2026-07-22 →):**

| # | Plan | Delivers |
|---|------|----------|
| 12 | **Straddling ambiguities & the right decomposition** (`2026-07-22-plan-12-ambiguity-blocks.md`, DELIVERED) | Fixes the latent CS §3 block-combinatorics bug (exact-pair condition missed straddling-overlap chains — repro `k⟨x,y⟩/(xx,yy,xyx)`, Bardzell HH wrong from degree 2 vs bar/minimal oracles); adds `right_decomposition`; corrects the Bardzell and CS odd (2-term) differentials to CS §4 `f_n` (first term from the right factorization); lifts the CS non-quadratic non-monomial `NotImplementedError` (the Plan-04 stretch item) with a new battery: straddle-monomial, QCI(3,2) with CS §7.2 φ-formula pins, cubic-tip-with-tail — all gated by d²=0 + order condition + live bar agreement |
| 13 | **Minimal A^e engine: multi-vertex support** (`2026-07-22-plan-13-minimal-multivertex.md`, DELIVERED) | Fixes the silent zero-resolution bug on multi-vertex input (local-only radical formula swallowed the kernel); builds the corner-typed minimal **projective** resolution `P_n = ⊕ A^e·(e_v⊗e_w)` (local algebras keep the kernel-accelerated free path bit-for-bit); nilpotent-closure guard refuses non-path-type bases loudly; validated vs bar on kA_2 / commutative square / kZ_3·rad² and by corner Betti = Bardzell chain counts `6,5,2,1,0` on kQ/(abc,cde) (independent re-derivation of Plan 12's straddle chain); `complexity` now exact for multi-vertex; `deepen` stays local-only (loud refusal) |
| 14 | **CS operations beyond the quadratic window** (`2026-07-23-plan-14-cs-operations.md`, DELIVERED) | Rebuilds the CS↔bar comparison map as a degreewise homotopy lift `Φ_n = h∘Φ_{n-1}∘d_n` — the closed-form block map silently failed the chain-map identity for ANY tip of length ≥ 3 (third uniform-zoo latent bug) and non-monomial n ≥ 3 refused outright; cup/bracket transport now serves every admissible presentation inside the bar window, plus the previously missing cap transport (`PhiHom` covariant collapse, `cap_of_cs_classes`), gated by chain-map/roundtrip/cup-route/unit-cap/module identities on the Plan-12 algebras |
| 15 | **Corner-mode checkpoints for `deepen`** (`2026-07-23-plan-15-deepen-corner-checkpoint.md`, DELIVERED) | Lifts the Plan-13 boundary: the checkpointed resumable driver now serves multi-vertex algebras (cluster-scale scans). Corner data is deterministic from `(A, prime)` and rebuilt on resume; the payload persists only the extra per-degree corner `tags`; HH finalization gains the corner contracted-complex branch (mirroring `minimal_homology_dims`); cross-mode ckpt_dir reuse refuses loudly (`QuiverlabError`); local path bit-for-bit unchanged. Validated vs the batch corner engine on CN(3,2) over 4 primes, kA_2 termination, non-monomial `kQ/(ab−cd)`, memory-wall/finalize-only/resume parity |
| 16 | **HH cohomology from the minimal/corner resolution** (`2026-07-23-plan-16-minimal-cohomology.md`, DELIVERED) | `minimal_cohomology_dims`: `Hom_{A^e}(−, A)` on the SAME minimal projective resolution — deep HH^• for **any** f.d. algebra over GF(p), the engine's second deep oracle. The coh collapse acts `a·w·b` (homology is `b·w·a`); corner blocks are the swapped tag `e_v A e_w` (load-bearing: kA₂'s coh corner is 1-dim where the homology corner is 0). Validated vs bar coh over 4 primes (local + multi-vertex zoos), Happel/Künneth `[1,0,0]` pins, degreewise vs `cs_cohomology_dims` to depth 8 on the GF(3) quantum CI, truncation-prefix semantics |
| 17 | **CS canonicalization** (`2026-07-23-plan-17-cs-canonicalization.md`, DELIVERED) | The order-condition correction γ (unique only mod the solve nullspace; nullity grows with degree on the quantum CI) is reduced to the free-variables-zero coset representative (`fields.linalg.reduce_mod_nullspace`, wired in `_d_general`) — CS differentials are **byte-reproducible by construction**. The 7 former `xfail(strict=False)` pins (2×3 bank byte batteries + the paper d₂ coefficient pin) are strict plain tests; an adversarial-solver gate shifts the solve by a nullspace vector and demands byte-identical differentials plus d²=0/order. No-op on prior outputs (Plan-04 stretch item E2) |
| 18 | **Standing-zoo diversity audit** (`2026-07-23-plan-18-zoo-diversity.md`, DELIVERED) | The standing zoo gains the two shapes whose absence hid both 2026-07-22 bugs: zoo records may carry `vertices`/`arrows` (multi-vertex reduction systems — previously inexpressible), five curated records added (straddle_xx_yy_xyx / straddle_xx_yy_yxy: the catalog's FIRST monomial records, mixed tips {2,3}; cn_3_2, comm_square, line_abc_cde), each live-certified (dim + minimal≡bar over 2 primes; line quiver by the Bardzell chain-count pin 6,5,2,1,0). Diversity gates pin ≥2 straddling-monomial + ≥3 multi-vertex records. Batch surface follows: specs carry quiver data (legacy specs byte-unchanged) and `_analyze_open` stops unit-adapting multi-vertex algebras (would trip the Plan-13 radical guard). New debt filed: open-zone scan cohomology (stale None since Plan 16) |
| 19 | **Field generality of engine-backed invariants** (`2026-07-23-plan-19-field-generality.md`, DELIVERED) | The five engine-backed invariants compute over every exact Domain (QQ, CC, QQi, GF(p^n)); GF(p) keeps the engine byte-unchanged. `cyclic_homology`: generic (b,B) mixed complex on the normalized bar basis (`hochschild/cyclic.py`) — any unital algebra, no quiver needed; gated by GF(p) engine parity, exact mixed-complex identities, and an independent char-0 second model (Connes λ-complex, Loday 2.1.5). `complexity`: relative-Tor (Cibils) Betti complex (`invariants/betti.py`) — H_n = the minimal resolution's rks[n] over every field, GF(p)-parity-gated incl. multi-vertex + straddling records. Frobenius family (`invariants/frobenius.py`): conclusive socle-permutation criterion (Skowroński–Yamagata), verified socle-dual form, ν = G⁻¹Gᵀ self-certified by its defining identity; `is_symmetric` upgrades to the definitional "ν inner" (Schwartz–Zippel sweep, loud when inconclusive — never a silent wrong answer). `_require_prime_field` retired; the residual refusal (structure-constants algebras off GF(p), path-basis-needing invariants) is honestly worded — no "later phase" promise survives in `src/` (gated by a grep test) |
| 20 | **Native deep-degree CS cup** (`2026-07-24-plan-20-native-cs-cup.md`, DELIVERED — branch `plan-20-native-cs-cup`, merged 2026-07-25) | The Hochschild cup PAST the bar 2n+1 window, computed natively on the CS small model — no bar object ever built. A comparison-lifted diagonal `Δ: P → P ⊗_A P` (`resolutions_cs/diagonal.py`) modelled as a **double-PELT** with the Koszul-signed tensor differential `d_p⊗1 + (−1)^p·1⊗d_q` (sign on the second summand only), built degreewise by a **per-degree lift-solve** `d^{P⊗P}·Δ_n(σ) = Δ_{n−1}(d_n σ)` — the structural clone of `_d_general`'s correction solve (`solve` + `reduce_mod_nullspace`, canonical/byte-reproducible; an inconsistent solve raises the identical loud `NotImplementedError`, no contracting homotopy needed). The sign-free `a·w·b` cup collapse with two interior values lives in `cup.py`; `Comparison.cup_of_cs_classes(engine="auto"/"native"/"transport")` routes native past the window (auto), transported in-window (byte-unchanged, still its refusal at `"transport"`). Oracles: **Leibniz** (the sign arbiter, exact GF(p)), the in-window transported cup **anchor** (native ≡ transport mod coboundary), a bumped-window transport **bridge**, `(d^{P⊗P})²=0` + chain-map + graded-commutativity/associativity, byte-reproducibility, one QQ domain-generic smoke, and the first genuine multi-vertex diagonal (`comm_square`). The **cap** is Plan 21 (same Δ); the **bracket** stays transported/window-bounded by design (needs the CS brace/circle machinery) |
| 21 | **Native deep-degree CS cap** (`2026-07-25-plan-21-native-cs-cap.md`, DELIVERED — branch `plan-21-native-cs-cap`, merged 2026-07-25) | The Hochschild cap PAST the bar window, computed natively by reading the SAME Plan-20 lifted diagonal Δ the HOMOLOGY way — no new construction. The covariance flip is the whole idea (as in the minimal engine's Plan-16 note: cohomology acts `a·w·b`, HOMOLOGY acts `b·w·a` on SWAPPED corner tags): `(f∩z)(ρ) = Σ_Δ coeff·b_c·x·b_a·f(τ)·b_mid` (`resolutions_cs/cap.py`), sign-free — all Koszul sign already lives in Δ; f eats the degree-p first factor τ, ρ survives, read against ρ's `hom` corner. `Comparison.cap_of_cs_classes(engine="auto"/"native"/"transport")` routes native past-window (auto), transported in-window (byte-unchanged, still its refusal at `"transport"`). The sign convention is ARBITRATED, not assumed: the in-window native ≡ transported cap **anchor** (the non-commutative quantum CI distinguishes `b·w·a` from `a·w·b`), the exact **unit cap** `1∩z=z`, the exact **cap-Leibniz** `b(f∩z)=(-1)^{p+1}(δf∩z)+(-1)^p(f∩bz)`, and the **module identity** `(z∩f)∩g ~ z∩(f∪g)` via the NATIVE cup all hold simultaneously. Degree edges: `n=p → C_0`; `p>n` raises. Bumped-window bridge, multi-vertex (comm_square + cn_3_2), QQ domain-generic smoke. Folded-in Plan-20 review backlog: hardened QCI skip-guards (dims regression now FAILs), element-wise CS-basis bridge, session-scoped Δ_4 fixture, `native_cup` chain-cache routing, constructed inconsistent-lift refusal. The **bracket** stays transported/window-bounded by design |
| 09 | **Server tier** (`2026-07-18-plan-09-web.md` + 2026-07-24 interface-drift amendment, DELIVERED — branch `plan-09-web`, merged 2026-07-25) | The public no-code web GUI `webapp/`: FastAPI `create_app` (server-rendered bilingual EN/ES pages + JSON `/api/`) with vendored self-hosted KaTeX; a spawn-child, resource-capped (`setrlimit` CPU+AS on Linux) worker fleet consuming a single SQLite/WAL jobs table that is the only queue+state store (atomic `BEGIN IMMEDIATE` claim, startup requeue of stranded `running` rows, graceful stop, retention sweep). Two compute tiers — instant synchronous (hard wall-net → converts to queued) and queued jobs with a permalink — plus the spec-§17 email magic-link big-job tier (HMAC single-use expiring token, per-email-hash caps, SMTP relay, emails never logged/stored-past-notice). Library-introspected family catalog + versioned request schema (family AND quiver kinds); all mathematics delegated to `import quiverlab` (no `eval`/`exec`/user-code; `quiverlab.engine.*` never imported); strict CSP + genericized errors + salted-hash IPs + ULID-gated traversal-proof downloads + library-sourced citations (`/literature`). Feedback (honeypot + per-IP cap) + header-token constant-time admin view. Deploy: Dockerfile (`.[web,fast]`, non-root, healthcheck), docker-compose (required prod secrets, `stop_grace_period`), Caddy TLS, PROVISIONING runbook. Base `pip install quiverlab` stays lean (webapp excluded from the wheel; web deps `[web]`-gated). Executed subagent-driven with paired adversarial critics; the whole-branch review's four cross-layer majors (async-error genericization, header admin token, last-hop XFF trust, lean base wheel) fixed |
| 22 | **Verification transparency** (`2026-07-25-plan-22-verification-transparency.md`, DELIVERED — branch `plan-22-verification-transparency`, merged 2026-07-25) | Displays in the repo that every shipped feature is unit tested and says HOW — both against QPA and against theory from the literature. Audit-first (1377 tests: fast 504 / deep 868 / qpa 5; every fact traced to an actual test file). `docs/verification.md` (wired into the mkdocs nav after "Under the hood"): the **two oracle classes** in Marco's framing — (1) theory/literature pins on constructed examples (Happel hereditary vanishing; BGMS/Bergh–Erdmann quantum CI; classical `k[x]/(x^n)` and Künneth commutative-CI values; the read-only hanlab bank byte-level closed forms; self-certifying `d∘d=0`/order/Leibniz/`λ`/`ν` identities; Connes λ-complex + relative-Tor Betti second models; Bardzell chain-count pins), and (2) cross-engine + external agreement (bar ≡ minimal ≡ Bardzell ≡ CS degreewise over `{32003,2,3,5}`; numba vs pure parity gated by the twice-run deep leg; live GAP/QPA recomputation of `HH^n = Ext^n_{A^e}(A,A)` and module self-Ext). Subsystem → oracles → test-file table; the marker/bucket scheme; the CI matrix; and an **honest-scope** section that names exactly where QPA cannot compare (cup/cap/bracket, cyclic homology, the CS resolution, deep degrees, Frobenius/Nakayama, CC/`GF(p^n)`, distinct-module Ext) and the theory oracle that covers each instead. Concise README "How quiverlab is verified" section. No `src/` change |
| 25 | **Webapp result cache** (`2026-07-25-plan-25-webapp-result-cache.md`, DELIVERED — branch `plan-25-webapp-result-cache`, merged 2026-07-25) | Never recompute a known example. A `result_cache` table (new, inside the single SQLite/WAL store — not a jobs column) keyed by the CANONICALIZED request + library version (`webapp/server/cache.py`: sorted-keys JSON so dict order is irrelevant and a tuple/list round-trip collides; version bump invalidates naturally). All three tiers check the cache FIRST — a hit replays the finished result instantly, across users, with zero recompute and not even a rate-limit charge. Big-job crux: the cache check is the first statement in `submit_big`, BEFORE `big_jobs_enabled`, so a cached big example is served with NO token minted, NO email sent, NO `email_hash` stored (email verification gates the COST of computing, not access to the mathematics) — even when SMTP is off. The worker records on success only; a cache row "pins" its finished job against the ordinary retention purge (artifacts survive to back replays) until it is both past the cutoff AND LRU/version-evicted (`sweep_cache_once` beside `sweep_once`). Math-only rows (no email/ip/token). Config `QLWEB_CACHE_ENABLED`/`QLWEB_CACHE_MAX_ENTRIES`. In-flight dedup documented benign and skipped (needs a jobs-schema change); `cache_put` idempotent so concurrent identical completions never crash |
| 23 | **Module-theoretic surface — engine slice** (`2026-07-25-plan-23-module-surface.md`, DELIVERED — branch `plan-23-module-surface`, merged 2026-07-25; Tier 1b items 1/2/4) | The classical module invariants, no code required. The opposite algebra `A^op` as a first-class `Algebra` (reversed quiver, transposed structure constants, reversed labels/relations, index-preserving cached involution — `modules/opposite.py`); the duality `D = Hom_k(−,k)` on arbitrary modules (`modules/duality.py::dualize`, transpose every action + reverse every label; `D(P_v^{op}) ≅ injective(A,v)` reconciles the previously-implicit `I_v = D(Ae_v)`); the corner-transpose `Tr M = coker(d₁*)` from the minimal presentation; and the **Auslander–Reiten translates** `τM = D(Tr M)`, `τ⁻M = Tr(DM)` (`M.tau()`, `M.tau_minus()`). **Injective resolutions** + `injective_dimension(M) = pd_{A^op}(DM)` as the dual of Plan 05 on the same op+D engine (`modules/injective.py`). Exact module isomorphism (`modules/hom.py::is_isomorphic`, invertible-hom certificate; generic-rank/Noether–Deuring over char-0; loud when a big GF(pⁿ) exceeds budget). Oracles: theory/literature (hereditary Coxeter transformation `dim τM = Φ⁻ᵀ dim M`, kA₂/kA₃ AR tables, self-injective ⇒ inj.dim ∈ {0,∞}, Nakayama τ-orbits, `max_v inj.dim S_v = gl.dim`, `(A^op)^op = A`, `D∘D = id`) AND **QPA/GAP** `-m qpa` (`DTr`/`TrD` dim-vectors + `IsomorphicModules` on a translated module, `ProjectiveResolution`/`DualOfModule` resolution terms, `InjDimensionOfModule`) across the zoo incl. the Plan-18 multi-vertex `line_abc_cde`. Deferred: Tier 1b item 3 (no-code module input in GUI/webapp) — specified in the plan doc, NEXT slice |
| 24 | **Left modules alongside right, right as default** (`2026-07-25-plan-24-left-modules.md`, DELIVERED — branch `plan-24-left-modules`, merged 2026-07-25; Tier 1b) | The module surface takes `side="right"|"left"`, right the (byte-unchanged) default. A left A-module IS a right A^op-module, so the Plan-23 op+D engine makes the side a **presentation-only tag** on `Module` (`self.side`, `self.base_algebra`, `self.with_side`) — every algorithm (rad/top/soc, Hom/End/Ext, `is_isomorphic`, projective AND injective resolutions, τ/τ⁻, dimensions) reads only `(algebra, action)` and is blind to the side, so left construction routes through `A^op` with **no mathematics forked**. `Algebra.simple/projective/injective/module` gain `side=`; `D` and `Tr` are now **side-aware**, exchanging left↔right over the SAME base algebra (their classical contravariant form), so `M.dualize()` of a right A-module is now a LEFT A-module (representation byte-identical, only the tag/repr move; `injective.py`, τ, τ⁻ numerically unchanged). `is_isomorphic`/`hom`/`ext` refuse loudly across sides/algebras (`_assert_comparable`, before the dim-vector fast-paths — a left `S_v` and a right `S_v` with equal dim vector still refuse). The honest asymmetry: on kA₂ right `P(1)=e_1A` has dimvec `{1:1,2:1}`, left `P(1)=Ae_1` has `{1:1,2:0}`. Two Plan-23 tests spelling `I_v=D(Ae_v)` via `A.opposite().projective(v).dualize()` moved to the honest `A.projective(v, side="left").dualize()` (documented). Oracles: ASS2006 (left/right duality D, AR translates) AND **QPA** `-m qpa` crosschecking left τ/τ⁻/inj.dim/resolutions by FEEDING QPA THE OPPOSITE ALGEBRA (right-module native) |
| 26 | **No-code module input (GUI + webapp)** (`2026-07-25-plan-26-no-code-modules.md`, DELIVERED — branch `plan-26-no-code-modules`, merged 2026-07-25; Tier 1b item 3) | Every representation theorist specifies a module with zero code and reads off the classical module invariants. Webapp request **schema v2** adds a `module` block — explicit `{dims, maps, side}` (one exact-entry matrix per arrow) OR a zero-typing `{builtin: {kind: simple|projective|injective, vertex, side}}` — plus `ext_target` (the N in Ext), guarded to schema 2 and canonicalizing through Plan-25 `canonical_key` (side default explicit ⇒ omitted `side` and `"right"` collide; a non-module request's key stays byte-identical via a `model_dump` that drops the absent block). New module compute kinds — dimension vector, rad/top/soc, Ext^n(M,N), τ/τ⁻, projective & injective resolution, projective & injective dimension — served by ALL THREE tiers (instant/queued/big-cached), each block carrying references/citations (GSZ2001/ASS2006). Matrix entries are exact DATA (ints or strings like `"1/2"`, never evaluated, never floats); a relation-violating module raises the library's loud error surfaced as a clean typed **4xx**, never a 500. Per-arrow BLOCK matrices (dim[target]×dim[source]) expand to the full vertex-ordered action `A.module` consumes, correct for right (A) and left (A^op) alike; module-aware tier sizing (`estimator.sizing_dim`) routes big modules off the instant tier like oversized families. GUI: a "Module (no code)" panel on the Pyodide canvas — per-vertex dimension picker, per-arrow matrix grid whose dimensions follow the arrow's source/target dims, a right/left side toggle, and S(v)/P(v)/I(v) pick-lists; the client runner (`docs/gui/runner.py`) carries the twin dispatch (results render like existing blocks, MathJax + citations). Interface-freshness pins the module surface |

**Planned next (2026-07-25, post-merge):** Plans 21–26 are all DELIVERED AND
MERGED (rows above): the native cup **and** cap past the bar window (the
Gerstenhaber bracket stays transported/window-bounded by design — it needs the
CS brace/circle machinery), the verification-transparency surface
(`docs/verification.md` + the citation rule), the complete module-theoretic
surface (τ/τ⁻, injective resolutions/dimension, left+right modules with right
the byte-unchanged default, no-code module input in the GUI and webapp — the
Tier-1b vision delivered end-to-end), and the webapp result cache (a cached big
example needs no email). Next on "continue"
(`docs/plans/DEEPER-ENGINES-BACKLOG.md`, topmost unchecked first): the Tier-2
remainder — Ext-algebra / Yoneda-ring presentations, HH cohomology ring
structure + support varieties, BV structure, periodicity certificates,
single-degree HH mode, Han's-conjecture batch campaigns, open-zone scan
cohomology (stale `None` since Plan 16), A∞ (Kadeishvili), performance (numba
corner kernels; GF(p^n) fast engine), deeper-engine GUI surfacing, native
AR-quiver — plus the Tier-3 webapp polish list.

Standing constraints for every plan: exact arithmetic only, floats fail loudly
(AST gate enforces); read-only banks are never modified; ≤2 parallel agents during
execution; path composition is left-to-right (`a*b` = first `a`, then `b`,
requiring `target(a) = source(b)` — Assem–Simson–Skowroński convention); internal
currency is the unit-adapted structure-constant `Algebra`.

Interface freeze between plans: later plans consume `Domain`
(`coerce/add/sub/neg/mul/inv/is_zero/eq/characteristic`), `linalg.rank/nullspace/rref/solve`,
`Quiver`, `Relation`, `Algebra` (`.domain, .dim, .T, .unit, .multiply, .unit_adapted,
.basis_labels`), and `HHTable` exactly as defined in Plan 01.

**Standing constraint added 2026-07-18 (Marco):** `docs/internals/` — "under the
hood for algebraists" — explains at low coding level how each object is
represented and how each computation runs (chapters for Plans 01–02 exist; the
format is fixed there). Every subsequent plan's acceptance task adds/updates the
chapters for what that plan introduces (Plan 03: groebner + ReductionSystem;
Plan 04: CS resolution + comparison maps; Plan 05: modules — module
representation and module resolutions; Plan 06: families/batch; Plan 07:
viz/trace; Plan 09: webapp). Plans 03/04/09 were committed before this rule:
their executors receive the obligation at dispatch. Chapters ship on the Plan 08
docs site as an "Under the hood" section.

**Standing constraint added 2026-07-18 (Marco):** citations subsystem (spec §3.9).
Plan 06 ships the registry core + `references.bib` + family refs; Plan 07 renders
the References section in worked-steps PDFs; Plan 09 surfaces `result.references`
on result pages and ships `/literature` + the `literature` feedback category
(web spec §16). Plan writers for 06/07 inherit this at plan-writing time.

**Standing constraint added 2026-07-25 (Marco):** verification transparency
(Plan 22). `docs/verification.md` is the single living record of HOW quiverlab is
verified — the two oracle classes (theory/literature pins on constructed examples;
cross-engine + QPA/GAP agreement wherever QPA implements the feature), a
subsystem → oracles → test-file table, the marker/bucket scheme, the CI matrix,
and an honest-scope section. **Every subsequent plan's acceptance task adds its new
oracles to this page** (the same obligation as the internals chapters): when a plan
ships a new engine, invariant, or operation, it extends the verification tables
with the oracle that guards it and the test file that runs it, and — where QPA
cannot compare — names the theory oracle that covers that ground. The page must
never overclaim: if a subsystem lacks an oracle, it says so. **Every
literature/theory oracle carries its citation** (Marco, 2026-07-25 follow-up):
cite the literature we test against at the precision the repo can verify — author,
year, venue, and a theorem/example/proposition number *only when it is actually
recorded* (test/docstring/plan/bank), never guessed — reusing a `references.bib`
registry key (`src/quiverlab/citations/`) where one exists; the page keeps a
References section rendered consistently with `quiverlab.bibliography(...)`.

---

# The two metagoals and the coverage program (Marco, 2026-07-25 →)

quiverlab is measured against two standing metagoals (stated at the top of the
README):

1. **No code required.** Every computation must be reachable with zero code:
   GUI canvas + no-code module panel, config-file/container workflows, rendered
   results/PDF reports. Every coverage-program plan below carries a "GUI story"
   as part of acceptance — a feature is not done until a representation
   theorist can click it. (UX reference points: Enomoto's FD Applet and
   Geuenich's String Applet.)
2. **Any computation done in representation theory** (of finite-dimensional
   algebras). Whatever a representation theorist computes in a paper should be
   computable here, exactly, with certified, oracle-tested results — or the
   page says honestly why not yet.

Evidence base: `2026-07-25-metagoal-coverage-deep-research.md` (three
deep-research reports: classical-curriculum audit; modern research practice;
software-systems inventory incl. a QPA `.gd`-level feature diff). Key
structural findings: (i) quiverlab already owns the expensive kernels — the
top-tier gaps are **composition layers** (morphisms, decomposition), not new
math engines; (ii) τ-tilting / torsion lattices / stability walls / 2-term
silting are **one engine** (AIR mutation), and no maintained general tool
exists anywhere — QPA included; (iii) the gentle/string world and the AAG /
surface derived invariants are white space outside two web applets.

## Coverage program (phases; each becomes one or more numbered plans)

**C1 — The categorical glue (P1, unlocks everything).** Module homomorphisms
as first-class objects (Hom BASES, not dims; kernel/image/cokernel as Modules;
composition; split-mono/epi; short-exact-sequence + is-split utilities;
projective cover / injective envelope as maps); direct sums with
inclusions/projections + `is_direct_summand`; `End(M)` as an algebra;
**Krull–Schmidt decomposition** into indecomposables with multiplicities +
`is_indecomposable` (End-local): GF(p)/GF(p^n) via idempotent lifting through
the radical (MeatAxe-style, Las-Vegas-exact), char-0 via exact idempotent
lifting, loud beyond budget; radical/socle SERIES + Loewy-diagram display;
composition-factor listing. Internals exist (`radtopsoc.submodule/quotient`,
`hom_space` matrices) — largely a surfacing job. Oracles: QPA
(`HomOverAlgebra`, `DecomposeModule`, `AlmostSplitSequence` inputs), theory.

**C2 — Forms, roots, recognition (P1, cheap off the existing Cartan).**
Euler bilinear form `⟨d,e⟩ = d C^{-T} e^T` + Tits quadratic form + definiteness
⇒ Dynkin/Euclidean/wild classification for hereditary; **type DETECTION**
(recognize A/D/E/affine from a quiver); Kac/Gabriel root combinatorics —
indecomposable dimension vectors = positive roots (hereditary), real/imaginary
roots, reflections; recognizer batch (`is_hereditary`, `is_gentle`,
`is_string`, `is_special_biserial`, `is_nakayama`, `is_radical_square_zero`,
`is_gorenstein`, `is_semisimple`, `is_weakly_symmetric`, `is_connected`,
`is_basic`, blocks); quiver graph predicates. Oracles: QPA unit forms +
Gabriel/Kac literature pins.

**C3 — Auslander–Reiten theory completed (P1).** Almost-split sequence with
its MIDDLE TERM `0→τM→E→M→0`; irreducible maps / rad(M,N)/rad²; **AR-quiver
knitting** (semi-decides rep-finiteness; loud budget cap); Nakayama functor
`ν(M)` named; stable Hom. Needs C1. Oracles: QPA `AlmostSplitSequence` /
`PredecessorsOfModule`; ASS/ARS worked examples.

**C4 — The τ-tilting engine (P1, the modern flagship; one engine = four
areas).** `is_tau_rigid` (= `Hom(M, τM) = 0`, composable TODAY); support
τ-tilting pairs via mutation BFS from `(A, 0)` (semi-algorithm: complete iff
τ-tilting finite — loud cap otherwise); g-vectors from minimal projective
presentations; exchange graph + Hasse poset; torsion-class lattice with
brick-labelled arrows; bricks/semibricks ↔ wide subcategories; 2-term silting
relabeling; King θ-stability + the **wall-and-chamber picture drawn live in
the GUI for n = 2, 3** (the killer no-code demo); maximal green sequences.
Needs C1 (decomposition) for summand handling. Oracles: theory (AIR/DIJ/DIRRT
identities — e.g. #sτ-tilt = #torsion classes = #chambers), FD Applet /
feisele-blueprint spot values; QPA cannot compare (state on the verification
page).

**C5 — The gentle/string subsystem (P1 for the class; complete combinatorial
answers).** String & band module classification (Butler–Ringel), string-module
syzygies/τ (SBStrips-style), the special-biserial AR quiver, **AAG derived
invariant** + surface-model data (Opper–Plamondon–Schroll) — white space, ties
to our HH (Redondo–Román oracles already researched); Gorenstein-projectives +
singularity category for gentle (Kalck, combinatorial); Brauer graph algebra
constructor family. Oracles: literature pins + SBStrips/String Applet spot
values + our own bar/CS engines.

**C6 — Homological-dimensions family (P1/P2).** Finitistic dimension (findim),
dominant dimension, Igusa–Todorov φ/ψ, Gorenstein dimension + `is_gorenstein`,
delooping level (Gélinas), Ω/τ-periodicity tests. Composes on existing
syzygy/Ext engines. Oracles: QPA (domdim/Gorenstein/periodicity) + literature.

**C7 — Tilting & new-algebra constructions (P2).** `is_tilting`/cotilting +
Bongartz completion + complement mutation; minimal left/right approximations;
**`End(M)` as a quiver-with-relations algebra** (tilted algebras); Gabriel
quiver recovery / basic-ization of structure-constant algebras (unlocks
End(M)-as-algebra and consuming kG-block basic algebras); one-point
(co)extensions `A[M]`; repetitive algebra; Jacobian-algebra constructor from
(Q, W). Oracles: QPA tilting/approximations; ASS VI examples.

**C8 — Geometry, derived fingerprints, complexes (P2/P3).** Kac canonical
decomposition (Schofield / Derksen–Weyman); orbit dimensions + Voigt
(`dim Ext¹(M,M)` rigidity); degeneration/hom order for finite type; the
"derived fingerprint" compare-two-algebras panel (HH/HC/Cartan/center — all
already computed); user-facing complex/chain-map API (cones, truncations,
homology of complexes); Hall numbers over GF(q) (P3). Oracles: QPA
degeneration + literature (Zwara/Bongartz).

**Explicit out-of-scope (stated, with reasons):** species/valued quivers
(different foundation: division-algebra bimodules — recorded, not planned);
DWZ potential mutation / cluster categories (the Jacobian-algebra CONSTRUCTOR
is in scope above; the mutation engine is a separate field — Sage covers its
combinatorial shadow); DT/BPS/moduli-counting; group-algebra MeatAxe at scale
(we CONSUME block basic algebras, we don't derive them); coalgebras/comodules
(Simson vols 2–3).

**Standing rules for every coverage plan:** oracle-first acceptance per Plan 22
(verification page + citations); a GUI story per feature (metagoal 1); honest
semi-decision contracts ("complete: N" vs "budget hit — uncertified") wherever
termination is equivalent to finiteness (knitting, τ-tilting BFS).
