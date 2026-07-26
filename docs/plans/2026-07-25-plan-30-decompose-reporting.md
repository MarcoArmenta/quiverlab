# Plan 30 — Krull–Schmidt decomposition + GUI/report completeness (Marco feedback)

**Date:** 2026-07-25. **Branch:** `plan-30-decompose-reporting` (worktree off
main; merges AFTER Plan 29 — the Tor-GUI wiring here is skip-guarded on
`quiverlab.modules.tor` so it activates in the merged tree). **Source:** Marco's
GUI/PDF feedback, 2026-07-25 (verbatim intent below), + coverage-program C1
(the decomposition is its flagship item, pulled forward).

**Marco's asks (binding):**
1. AR translation must be sure the module is indecomposable; given a module we
   can compute its decomposition into indecomposables.
2. Tor and Ext must accept ANY module for both arguments in the GUI (not only
   simples/injectives/projectives).
3. Projective/injective resolution tables: per term, specify the indecomposable
   projectives (e.g. `P₁² ⊕ P₃`); REMOVE the `#summands` column; REPLACE the
   dim-vector column by that decomposition.
4. The PDF includes descriptions of all projectives and injectives (not
   simples — obvious).
5. Worked steps are incomplete — show the RENDERED OBJECTS; make the `.tex`
   downloadable; the pdf/tex bundle is exhaustive: it shows at the algebraic
   level the computations the engine does.

## Part A — the decomposition engine (`src/quiverlab/modules/decompose.py`)

- `is_indecomposable(M) -> bool` and
  `decompose(M) -> list[tuple[Module, int]]` (indecomposable summands with
  multiplicities, `⊕ Mᵢ^{mᵢ} ≅ M` certified by an explicit iso — the existing
  `is_isomorphic` certificate machinery closes the loop).
- Algorithm (exact, deterministic-first, house-honest): compute `End(M)` from
  `hom_space(M, M)` (bases exist). Split via idempotents: for spanning
  elements φ (and bounded small combinations), compute the exact minimal
  polynomial over the Domain; a factorization into coprime factors yields a
  Fitting/idempotent splitting `M = ker e ⊕ im e` (`radtopsoc.submodule`
  internals build the summands); recurse. `M` indecomposable ⇔ no split found
  AND the local-ring certificate holds (`End/rad(End)` a division algebra —
  over CC: dim 1 after nilpotent quotient; over GF(p): every generator
  nilpotent-or-invertible with the unit group budget). **Loud beyond budget**
  (never a silent "indecomposable" that merely failed to split — mirror
  `is_isomorphic`'s refusal style). Domain-generic; float-free.
- Surface: `Module.decompose()`, `Module.is_indecomposable()`,
  `Algebra`-level nothing new. τ/τ⁻ additivity NOTE in docstrings:
  τ(⊕Mᵢ) = ⊕τMᵢ — τ stays total; the REPORTING layer certifies (Part B).
- Oracles: QPA `DecomposeModule`/`DecomposeModuleWithMultiplicities`/
  `IsIndecomposableModule` live crosscheck (`qpa/crosscheck.py` additions —
  this plan owns that file) across the zoo + constructed direct sums
  (P₁⊕P₁⊕S₂, τ-translates, radicals of projectives); theory: simples/
  indecomposable projectives/injectives are indecomposable; Krull–Schmidt
  uniqueness (decompose twice / permuted constructions agree up to iso);
  `⊕`-roundtrip (build an explicit direct-sum action matrix, decompose,
  recover the parts up to iso). Tests `tests/modules/test_decompose.py`
  (deep) + `tests/qpa/test_decompose_qpa.py` (qpa).

## Part B — GUI/webapp completeness

- **Second-argument module editor:** the GUI's Ext target (and the new
  `tor`/`tor_target`) gets the SAME no-code per-arrow matrix editor as the
  main module panel (schema already accepts explicit `{dims, maps, side}` for
  `ext_target` — this is GUI-side only). Webapp schema: add `tor` to
  MODULE_RANGE_KINDS + `tor_target` (left-side module spec; validated side =
  "left") gated to schema 2, canonicalizing with the Plan-25/26 rules
  (absent → dropped, byte-stable existing keys). Runner dispatch `tor` via
  `quiverlab.modules.tor` — LAZY import + typed 4xx "requires the Tor engine
  (Plan 29)" when absent; tests skip-guarded on the import so the branch is
  green standalone AND in the merged tree.
- **`decompose` compute kind** (schema 2): result block = summand list with
  multiplicities + per-summand dim vectors + indecomposability certificates.
- **τ/τ⁻ blocks certify:** each AR-translate result now includes
  `indecomposable: true` OR the decomposition `M ≅ ⊕ Mᵢ^{mᵢ}` with the note
  "τ computed summand-wise (τ is additive)" — bilingual EN/ES.
- **Resolution tables (Marco #3):** in webapp blocks, GUI rendering, AND the
  PDF: columns = term `n` | `⊕_v P_v^{b_{n,v}}` (LaTeX `P_1^{2} ⊕ P_3`;
  injective resolutions use `I_v`). `#summands` and dim-vector columns
  REMOVED. Data source: the resolution terms' summand vertices (already in
  the result payloads or one call away — keep result-JSON backward-compat by
  ADDING the `summands` field and dropping the removed columns only in the
  RENDERINGS; cache keys untouched).

## Part C — the exhaustive pdf/tex bundle

- **`.tex` downloadable:** the server writes `trace.tex` alongside
  `trace.pdf`/`trace_steps.html`; artifact whitelist (`pages.py::_mount_download`)
  + job-page link + GUI download button + `quiverlab-hpc render --format tex`.
- **Projectives & injectives section:** every report (trace renderers +
  `hpc/report.py`) gains "The projectives and injectives of A": for each
  vertex v, `P_v` and `I_v` with dimension vector, Loewy layers (radical
  series — iterate `radical()`), and socle/top; simples omitted (Marco).
- **Rendered objects in worked steps:** module-computation traces show the
  actual algebra: resolution differentials as matrices (`_latex_matrix`),
  Hom-space bases where computed, the per-step rank bookkeeping the engine
  does — with the EXISTING eliding rules extended (beyond a size threshold
  render shape + rank + a note, never silently omit a step). HH traces keep
  their current pipeline; the new depth applies to module kinds
  (resolutions, Ext/Tor, τ, decompose). The bundle must let an algebraist
  REPLAY the computation by hand on small examples — that is the acceptance
  bar (a golden-file test does exactly that on kA₂: every differential of
  the S₁ resolution appears verbatim in the .tex).

## Acceptance

Part A green incl. QPA; Parts B/C green across tests/webapp + tests/gui +
tests/hpc renderer goldens (updated for the new table format + tex artifact +
P/I section); byte-stability: existing non-module cache keys unchanged; result
blocks only ADD fields; bilingual completeness (ES no English leak);
verification page: decomposition oracles row + honest budget note; backlog
entry (this file's item) ticked on this branch; suites green. Merge after
Plan 29, only when Marco asks. Follow-ups spawned, not done here: GUI panel
for decompose-then-pick-a-summand workflows; Ext/Tor cocycle-level rendering.
