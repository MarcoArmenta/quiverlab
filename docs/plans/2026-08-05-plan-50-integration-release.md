# Plan 50 — Integration sweep + v0.2.0 release gate (W6, the last plan)

> Executes the metaplan §5 P50 card and §7 release-gate acceptance. Written
> LAST by design, against the real merged tree (dev after P36–P49, suite
> 3432). Everything here is verification, reconciliation and packaging — no
> new mathematics.

## Ground state (all inputs delivered)

All 15 subplans merged to dev with per-plan critic rounds and recounts:
P36 (M2 fifth oracle class), P37 (categorical glue), P38
(forms/recognition/Koszul exposure), P39 (complexes/hyper-Ext), P40 (C6
homological dimensions), P41 (AR completion), P42 (spectral sequences), P43
(derived category), P44 (C7 constructions), P45 (τ-tilting engine), P46
(gentle/strings), P47 (quasi-hereditary/recollements), P48 (marked
surfaces), P49 (C8 geometry). Suite 3432 (fast 1556 / deep 1676 / qpa 189 /
m2 11), classes lit 882 / xeng 528 / selfcert 1080 / union 2118, audit
gates green at every merge.

## Task A — Cross-tier integration sweep (worker A, worktree)

The per-plan batteries verify each kind in isolation; this sweep verifies
the RELEASE-level contracts across all v0.2.0 kinds at once:

1. **Kind reachability audit**: every new compute kind
   (`ext_algebra`, `recognizers`, `homological_profile`, `ss_hochschild`,
   `derived_fingerprint`, `strings`, `quasi_hereditary`, `tau_tilting`,
   `orbit_geometry`, `tilting_check`, `almost_split`, plus the P44/P48
   construction presets) is (a) dispatchable in `hpc/spec.py` AND its
   Pyodide twin with the same block shape, (b) listed in the GUI pick-list
   (both gui.js copies byte-identical), (c) rendered by both JS renderers
   and `trace/results_html.py`, (d) i18n-complete: every `inv.*`/block key
   the renderers reference exists in BOTH `en.json` and `es.json`.
   Deliverable: `tests/release/test_v020_integration.py` — a data-driven
   audit (introspect the dispatch tables / pick-lists / i18n dicts, assert
   set-equality, no hand-maintained list drift) + an end-to-end smoke: one
   small algebra runs EVERY v0.2.0 kind through `run_spec` and the Pyodide
   twin dict-agrees on math subkeys.
2. **ETA scalars audit**: every scalar kind in both runners' `ETA_MODEL`
   dicts agrees (same keys both sides — the P45 merge nearly broke this).
3. **README metagoal scorecard**: C1–C8 rows stating delivered-by-plan, with
   the honest re-scopes (C4 brick scope, C8 Dynkin-only canonical
   decomposition, P48 unpunctured v1) — one table in README.

## Task B — Docs reconciliation (worker B, worktree)

1. **GUI-deferral ledger on the verification page**: one subsection listing
   what v0.2.0 deliberately defers to successors, with the plan-doc
   pointers: CE/Grothendieck/radical SS GUI presets beyond `ss_hochschild`
   (P42), the derived compare panel (P43), the free-form draw-a-surface
   canvas (P48.1 with punctures/self-folded), DWZ potential
   right-equivalence (P48.1), σ_A/τ-Hochschild machinery (backlog Tier 2).
2. **Internals chapters**: short "under the hood" chapters for the four new
   subsystems (spectral sequences, derived category, τ-tilting fan,
   surfaces→gentle) wired into mkdocs nav (strict build must stay green).
3. **CHANGELOG.md**: the v0.2.0 entry — one paragraph per plan, user-facing
   wording, honest-scope notes inline.
4. **Metaplan status flip**: every P-card marked delivered with its merge
   commit; deferred ledger cross-checked against the verification page.

## Task C — Full gate matrix (orchestrator, sequential background runs)

On the final dev tree (after A+B merge):

1. fast (local macOS) — green, no skips beyond the known 3.
2. deep, numba path — green.
3. deep, pure path (`QUIVERLAB_NO_NUMBA=1`) — green (engine parity).
4. `-m qpa` live (local GAP+QPA) — green.
5. `-m m2` live (local Macaulay2; honest skip if absent locally — CI job
   enforces).
6. `mkdocs build --strict` (with `DISABLE_MKDOCS_2_WARNING=true`) — exit 0.
7. Final recount + Plan-32 audit; verification page/README badge updated.

## Task D — Version + release packaging (orchestrator)

1. Bump `0.1.0` → `0.2.0` (pyproject + `__version__` + release-metadata
   tests + the interface-freshness pins that encode the version).
2. Curated example refresh: `webapp/precomputed/manifest.yaml` gains at
   least one cached example per flagship new kind (τ-tilting kA3,
   ss_hochschild, derived_fingerprint, a surface preset, orbit_geometry) —
   regenerated + WAL-checkpointed.
3. `dev` → `main` merge prepared locally + annotated tag `v0.2.0`.
4. **STOP AND ASK MARCO** before any push: pushing main triggers docs
   deploy; the tag triggers PyPI trusted publishing; container/desktop
   assets rebuild on release. These are outward-facing and irreversible —
   the gate ends with a go/no-go summary, never an unasked push.

## Acceptance (metaplan §7)

- [ ] C1–C8 delivered or re-scoped with documented reason (Task A3/B4)
- [ ] Five oracle classes audited, table ≡ live collection (Task C7)
- [ ] Every new kind clickable end-to-end EN+ES (Task A1)
- [ ] mkdocs --strict green; fast+deep(both paths)+qpa green (Task C)
- [ ] v0.2.0 prepared: bump, tag, assets plan — awaiting Marco's push
      approval (Task D)
