# Plan 33 — nontrivial literature examples at scale (batteries + papers)

- **Date:** 2026-07-26 · **Branch:** `plan-33-nontrivial-examples`
- **Driver:** Marco's feedback (saved to memory): test/paper examples are too
  small (kA₂/kA₃/single loops); the quantum ones are good — compute higher
  degrees; add nontrivial examples from the books and the literature. Also:
  positioning is representation-theory-first (already merged; the JSC interior
  pass rides this plan).
- **Research (two agents, 2026-07-26):** Plan 29 already pins most directions
  at small size — **Plan 33's job is SCALE**, not new theorems. Every value
  below was machine-verified during research or carries the stated citation;
  the honest-scope labels are BINDING (nothing labeled xeng gets a literature
  pin).

## The catalog (contract for the batteries)

Tier 1 (verified live in research; implement all):
- **Quantum CI generalized**: `k⟨x,y⟩/(x^a, y^b, yx − q·xy)` for
  (a,b) ∈ {(2,4),(3,4),(4,4),(2,5),(5,5)} over CC — Bergh–Erdmann coh
  `[2,2,1,0,…]` pushed to degree ≥ 8, homology `[a+b−1, a+b−2, a+b−2, …]`
  (`qci_hh_oracle`, Thm 3.1/3.2; char-0 branch only — the char-p branches need
  infinite fields we lack, documented). Verified in research: QCI(2,4) coh
  `[2,2,1,0,0,0,0,0,0]` to deg 8, hom `[5,4,4,4,4,4,4,4,4]`; QCI(4,4) coh to
  deg 7. GF(p) small-prime reshapes stay OUT (root-of-unity trap, Plan-29
  documented). Deep engine: minimal A^e internal for dim ≤ ~30 (deg 80 in 7s
  for dim 4), CS otherwise.
- **Preprojective Π(A₄), Π(A₅), Π(D₄), Π(D₅)** (dims 20/35/28/60): structural
  pins — dim (closed form dim Π(Aₙ) = C(n+2,3), verified A₂..A₅; the Dₙ form
  n(n−1)(2n−1)/3 matched D₄/D₅ only → pin builder values, promote the closed
  form only if a D₆ point is added), `is_selfinjective`, Loewy = h−1
  (Coxeter numbers 5/6/6/8) (`preprojective`, `assem_book`; Erdmann–Snashall
  for self-inj/LL). HH values: **xeng/QPA only** (no published table
  consulted). HH depth ≤ 6–7 and only where <60s (see envelope). Nakayama
  permutation: NOT pinned (convention check unresolved — flagged).
- **m-Kronecker m=3,4** (dims 5/6): `HH^• = [1, m²−1, 0, 0]` (Happel,
  `happel_question`; verified 8/15), Coxeter `t² − (m²−2)t + 1`
  (`lenzing_delapena_spectral` key if that is the registry name — VERIFY the
  key exists; else the registry's spectral key used by test_spectral.py).

Tier 2/3 (implement all that fit the runtime budget; cut from the bottom):
- **Canonical C(2,2,2,2,2)** (dim 19, 7 vertices): `dim HH² = t−3 = 2` — the
  first ≥2 case (same source as the live t=3,4 battery); Happel-trace
  cross-link.
- **Self-injective Nakayama at scale**: Taft Λₙ = kZ_n/J^n, n=5,6 (dims
  25/36): `HH_• = [n, n−1, n−1, …]` + HC even/odd = n/(n−1)
  (`taillefer_taft` — VERIFY the registry key name; Plan 29 added Taft keys).
  Symmetric Brauer stars kZ₄/J⁹ (36), kZ₅/J¹¹ (55): symmetry booleans +
  Bardzell depth (see envelope) (`skowronski_yamagata`; Rickard
  derived-invariance framing in docstring only, no numeric pin from it).
- **Exterior Λ(k³), Λ(k⁴)** (dims 8/16): Koszul via
  `modules.koszul.g_quadratic_certificate` (NOT `ext_algebra().koszul` —
  >120s trap), self-injective, Loewy n+1 (`priddy`, `froberg_koszul`).
  HH numeric values: **xeng only, char-qualified** (no published table).
- **Incidence of Boolean B₃** (dim 27): has 0̂,1̂ ⇒ order complex contractible
  ⇒ `HH^{≥1} = 0` at any depth (Gerstenhaber–Schack / Cibils; the registry
  keys the live diamond battery uses) — and higher HH is ~free to verify deep.
- **Presented trivial extensions at scale** (Plan 31): T(kD₄), T(kA₅),
  T(kA₆) (dim 42): certificate + 4 symmetry booleans + `C_T = C_A + C_Aᵀ` +
  `HH¹ ≠ 0` (`cmrs_split`) — all milliseconds; QPA crosscheck where cheap.
- **Commutative CI tensor**: `k[x]/(x³) ⊗ k[x]/(x³)` (dim 9),
  `Λ triple = k[x,y,z]/(x²,y²,z²)` (dim 8): symmetric ⇒ HHⁿ = HH_n;
  Künneth cross-link (`tensor_product`, `quantum_ci` for q=1 growth contrast).
- **Zoo multi-vertex deep**: 3-cycle, square, `line_abc_cde` to degree 24
  (verified ~free): extend existing depth pins.
- **Triangular-string A₅** (`redondo_roman_2014`, Ex 3): `HH^• =
  [1,5,0,0,0,2,0]` — vanish-then-revive at degree 5; anchor the revival on the
  minimal engine like the live A₃ test. slow bucket if >60s.

DEFERRED (recorded, not built): (D,A)-stacked Ex 1.2 and Cassidy non-Koszul
witness (Plan-27 feeders; build risk), Π(E₆)/Π(D₆) (AdmissibilityError at
tested bounds; certification cost grows — cluster), incidence ≅ S²,
Toupie figure-only example (never pin), Redondo–Román 2018 cup sets.

## Feasibility envelope (binding for bucket placement; timings are this-Mac upper bounds)

- Bardzell (monomial, internal route via `engine.hh_engine` +
  `MonomialPresentation.cyclic_nakayama`): deg 300 at dim 150–220 in <10s —
  THE depth showcase. **Python recursion cap**: default limit stops it
  (deg ~160 at dim 84; deg ~20 for k[x]/(x^100)); `sys.setrecursionlimit(50000)`
  lifts cleanly.
- Minimal A^e internal: QCI dim 4 → deg 80 in ~7s; **hangs at dim ≳ 30**
  (Π(A₅) unkillable-slow) — never use it there.
- `engine="cs"`: dim ≤ ~20 comfortable (T(comm-square) dim 18: deg 16 in 60s;
  Π(A₄) dim 20: deg 7 in 60s, deg 8 = 70s ⇒ slow); dim 42 T(A₆) reached
  deg ~32 (sparse structure helps — trust measured numbers only).
- bar/fast: degree ~5 ceiling on multi-loop local algebras — never a depth
  engine.
- decompose: fast to dim 16, deep to ~32 (GF(32003) safe: char > dim), slow
  at 48; τ-steps on Dynkin ≈ free.
- Budget: nothing in the deep bucket over ~60s/test; prefer session-scoped
  fixtures; anything beyond → `slow` marker (implies deep, opt-in) or the
  SUBMISSION.md step-4 cluster list.

## Design decisions

- **D1 — scale existing oracles; every pin cited or explicitly
  cross-engine/self-cert; honest-scope labels binding.** New tests carry the
  Plan-32 oracle-class markers correctly (literature pins →
  `oracle_literature`; engine-agreement-at-depth → `oracle_crossengine`;
  certificates → `oracle_selfcert`).
- **D2 — engine routing per the envelope above.** Bardzell depth tests wrap
  the computation in a temporarily raised recursion limit. Implement the bump
  as a **small engine-level guard** (context-managed raise/restore around the
  recursive Bardzell walk in `engine/resolutions_bardzell.py`) so users hit
  the same depth the tests do; pure/numba parity unaffected (Python-side
  only); a regression test pins deg-300 reachability on kZ₂₀/J¹¹.
- **D3 — `QuantumCI(q, a=2, b=2)`**: backward-compatible generalization of
  `families/quantum.py` (existing calls `QuantumCI(q)` byte-identical);
  relations `x^a, y^b, y*x − q·x*y` through the existing exact-coefficient
  tokens; catalog/discover entry updated; the raw-Quiver idiom stays valid.
- **D4 — Preprojective auto degree_bound**: a per-type bound table in
  `families/preprojective.py` from the verified builds (A₅:12, A₆:14, A₇:16,
  D₄:12, D₅:16 — use the prober's working values; keep the explicit kwarg
  override; E-types keep the loud AdmissibilityError with the bound hint —
  honestly documented as cluster-scale).
- **D5 — trace hygiene**: batteries/benchmarks run with `quiverlab.verbose =
  False` (check what tests/conftest.py already does — follow the existing
  pattern; also add `quiverlab_traces/` to .gitignore if absent).
- **D6 — papers**: JSC worked-examples section rebuilt around the research
  top-10 (flagships: Π(D₅) dim 60; QCI dim 16 cohomology-dies/homology-
  persists; Bardzell deg-300 at dim 220; Brauer star dim 55; canonical t=5;
  B₃ incidence dim 27; T(kD₄); 3-Kronecker wild; Taft Λ₅; Λ(k⁴) Koszul) with
  every number recomputed by replayable scripts in `paper-jsc/computations/`;
  the JSC INTERIOR gets the representation-theory-first pass (lead sections
  reordered/reframed; engines presented as serving representation theory).
  JOSS: fold 2–3 one-liners into Research impact WITHIN the 750–1750 word
  band (current 1669 — trim elsewhere as needed, keep gates green).
- **D7 — docs**: verification page gains the new oracle rows + honest-scope
  labels (xeng-only exterior/preprojective HH, the deferred list); backlog
  item + CLAUDE.md status + count resync at the end; SUBMISSION.md step-4
  cluster list replaced by the concrete feasibility list (preprojective HH at
  scale, Λ(kⁿ≥4) depth, ext_algebra/Koszul at scale, decompose ≳ 50,
  dim-30+ non-monomial past deg ~10, Π(D₆)/Π(E₆) builds).

## File ownership (three parallel agents)

- **C1 tests/**: extend `tests/engine/test_literature_p29.py`,
  `tests/resolutions_cs/test_battery_literature_p29.py`; new
  `tests/families/test_scale_batteries.py` (or split per family, agent's
  call); zoo-depth extensions; the Bardzell deg-300 regression test; qpa
  additions only where cheap. Markers per D1; buckets per envelope.
- **C2 src/**: `families/quantum.py` (D3), `families/preprojective.py` (D4),
  `engine/resolutions_bardzell.py` (D2 guard), `families/discover.py`
  catalog, `.gitignore` (D5 if needed). No other src.
- **C3 papers+docs**: `paper-jsc/` (worked examples + interior positioning +
  rebuild), `paper/paper.md` (one-liners, gates green),
  `docs/verification.md`, `docs/plans/DEEPER-ENGINES-BACKLOG.md` (add+tick
  the Plan-33 item), `SUBMISSION.md` step 4, CLAUDE.md status (counts left
  to the orchestrator).

## Acceptance

1. Batteries green in their stated buckets; nothing >60s in deep; `slow`
   used where needed; oracle-class markers correct (the Plan-32 audit gate
   recounts).
2. Zero unlabeled pins: every literature value cites its registry key in the
   docstring; xeng values marked as such.
3. `QuantumCI(q)` byte-identical; preprojective auto-bounds verified by
   building A₅/D₅ without kwargs; Bardzell deg-300 regression green.
4. Papers rebuilt (JSC PDF clean; JOSS gates green in-band); every paper
   number replayable.
5. Counts resynced everywhere; suites fast+deep+qpa green; strict docs green;
   merge only when Marco asks.
