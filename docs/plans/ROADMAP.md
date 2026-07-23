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
