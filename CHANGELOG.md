# Changelog

All notable changes to quiverlab are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org) (0.x during battle-testing; 1.0 at JOSS
acceptance).

## [Unreleased]

### Added

- **French and Chinese** join English and Spanish across the whole webapp/GUI:
  173-key catalogs `fr.json` / `zh.json`, every page mounted under `/fr` and
  `/zh`, a four-way header language menu (each language named in itself),
  localized big-job emails / verify links / job permalinks, and all
  family-catalog summaries and parameter help in all four languages. The i18n
  battery now gates key parity, placeholder-freedom, and `{url}`-slot survival
  for every catalog. (The worked-steps report and the HPC CLI remain
  English-only, as they were for Spanish.)
- **Search-first landing on `/draw`**: a search bar at the top of the GUI —
  type what you want to compute (in any of the four UI languages; the keyword
  index carries en/es/fr/zh synonyms), pick from the matching environments
  (related ones listed alongside), and a small curated example loads with the
  request pre-filled and computes on its own. 22 environments cover the whole
  compute surface (Hochschild/cyclic/products, invariants and recognizers,
  Ext-algebra/Koszulity, τ-tilting, gentle strings, module theory incl. AR
  translates, Krull–Schmidt, Ext/Tor, orbit geometry); every embedded example
  is live-validated against the real dispatch.

### Fixed

- `derived_fingerprint` killed the whole worked-steps report: the renderer
  built its HTML chunks but never returned them, so any request containing a
  `derived_fingerprint` block produced no report at all (the JSON record was
  fine). One-line fix plus a standing AST gate that every per-kind renderer
  ends in an explicit `return`.

## [0.2.0] — 2026-08-05

The *whole of representation theory* release, prepared at the P50 release gate (plans
P36–P50). It completes the C1–C8 coverage program and lands the five research axes
Marco named — a Macaulay2 oracle, spectral sequences, derived categories, the Coxeter
polynomial, and surfaces + category theory — while holding to the two metagoals: **no
code required** and **any computation in representation theory**. Verification grows a
fifth oracle class (live Macaulay2) alongside literature, cross-engine, self-cert, and
QPA.

### Added

- **Macaulay2 as an independent oracle (P36).** Single-vertex graded `kQ/I` (non-commutative
  Gröbner normal forms, dimensions, Hilbert series) and commutative examples (Ext/Tor, Betti
  tables, Poincaré series) are now cross-checked against Macaulay2. *Scope:* M2 sees no
  multi-vertex algebras and no Hochschild quantities — those stay with QPA and the theory
  oracles, and the bridge refuses out-of-scope inputs loudly.
- **First-class categorical glue (P37, C1).** Module morphisms (`ModuleHom`) with
  composition and kernel / image / cokernel as modules with connecting maps; short exact
  sequences with split tests, pushouts and pullbacks; direct sums with inclusions /
  projections; `End(M)` as an algebra; radical / socle *series* with Loewy display; and
  composition factors.
- **Forms, roots, recognition, Coxeter, Koszulity (P38, C2).** The Euler and Tits forms off
  the Cartan, with Dynkin / Euclidean / wild classification and type detection; a recognizer
  batch (`is_hereditary`, `is_gentle`, `is_string`, `is_special_biserial`, `is_nakayama`,
  `is_radical_square_zero`, and more); the exact Coxeter matrix and polynomial over ℚ (a loud
  refusal when the Cartan is singular); and a no-code Koszulity verdict.
- **Complexes, cones, and hyper-Ext (P39).** Bounded complexes of modules and of
  projectives, chain maps and homotopies, shift / truncation / homology, mapping cones and
  distinguished triangles, the derived-iso test (cone acyclicity), and `Hom_{D^b}(X, Y[n])`
  (hyper-Ext) via the Hom double complex — the substrate for the spectral-sequence and
  derived-category surfaces.
- **Homological dimensions (P40, C6).** Finitistic-dimension bounds, dominant dimension,
  Igusa–Todorov φ / ψ, Gorenstein dimension with `is_gorenstein`, and Ω / τ-periodicity
  certificates. *Scope:* the delooping level is deferred — its existential quantifier has no
  crisp bounded decision procedure.
- **Auslander–Reiten completion (P41, C3).** Almost-split sequences `0 → τM → E → M → 0`,
  irreducible maps and `rad / rad²` multiplicities, the Nakayama functor ν, stable Hom, and
  AR-quiver knitting. *Scope:* knitting is an honest semi-decision — complete iff
  representation-finite, loud budget cap otherwise; self-injective input is refused up front.
- **Spectral sequences (P42).** A double-complex / filtered-complex engine with exact pages,
  induced differentials, and a standing convergence certificate (`E_∞` equals total
  homology), plus four presets: the Hochschild `(b, B)` bicomplex (clickable as
  `ss_hochschild`), the radical / associated-graded filtration, and the Cartan–Eilenberg /
  Grothendieck change-of-rings sequence. *Scope:* only the `(b, B)` sequence has a no-code
  GUI this release (the other presets are library / HPC-config accessible); the general
  two-bimodule Grothendieck sequence is deferred (only the `U = B` change-of-rings case is
  implemented); Grothendieck acyclicity is a per-instance hypothesis check with a loud
  refusal.
- **Derived-category surface (P43).** `Hom_{D^b}(X, Y[n])` as honest chain maps, the derived
  AR translate `τ_{D^b} = ν[−1]` on perfect complexes (a loud refusal at infinite global
  dimension, per Happel), a tilting-complex *verifier* with `End(T)` as the Rickard
  derived-equivalent algebra, and a `derived_fingerprint` panel. *Scope:* the fingerprint
  speaks only in "distinguished / not distinguished by these invariants" — it is never a
  decider of derived equivalence (the 8-vertex cospectral trees are the pinned
  counterexample); the two-algebra compare panel is deferred.
- **Tilting and constructions (P44, C7).** `is_tilting` / cotilting with Bongartz
  completion, minimal add(M)-approximations, Gabriel-quiver recovery / basic-ization of
  structure-constant algebras (so `End(M)` / `End(T)` come back as quiver algebras),
  one-point (co)extensions, finite repetitive slices, and a Jacobian-algebra constructor from
  `(Q, W)` with a per-instance finiteness certificate. *Scope:* basic-ization is char-scoped
  (char 0 or char > dim) and split-only; tilting-complement *mutation* and the full
  (infinite) repetitive algebra are named successors.
- **τ-tilting engine (P45, C4).** Support τ-tilting pairs by mutation, g-vectors, the
  exchange graph and Hasse poset, the torsion-class lattice with brick / semibrick labels,
  King θ-stability with a **wall-and-chamber picture drawn live in the browser for n = 2, 3**,
  and maximal green sequences. *Scope:* the mutation search is an honest semi-decision —
  complete iff τ-tilting-finite, loud budget cap otherwise; QPA cannot compare.
- **Gentle / string subsystem (P46, C5).** Butler–Ringel string / band module
  classification, string-module syzygies and τ, the Avella-Alaminos–Geiss derived invariant,
  and a Brauer-graph-algebra constructor. *Scope:* string enumeration is complete only for
  the representation-finite (band-free) case; the AG invariant is a derived invariant,
  provably **not** complete.
- **Quasi-hereditary structure and recollements (P47).** Standard / costandard modules
  Δ(i) / ∇(i) for a chosen order, the Dlab–Ringel quasi-heredity test, good-filtration
  multiplicities with BGG reciprocity, the characteristic tilting module and Ringel dual, and
  recollements from an idempotent (`eAe`, `A/AeA`, the six functors). *Scope:* quasi-heredity
  is order-dependent (the no-code block reports the natural order); QPA has no comparable
  surface.
- **Marked surfaces (P48).** A marked-surface + triangulation datatype, the
  quiver-with-potential `(Q(T), W(T))` constructor, the gentle Jacobian algebra with
  finiteness certificates, and flip ↔ mutation certification — a new no-code *input* method
  (draw a surface, get the algebra), which then flows through every existing compute kind.
  *Scope:* v1 is unpunctured surfaces with non-empty boundary only — punctures, closed
  surfaces, and self-folded triangles refuse loudly (successor P48.1); flip is certified at
  the quiver level (full DWZ potential right-equivalence is P48.1); the free-form
  draw-a-surface canvas is deferred.
- **Representation geometry (P49, C8).** Orbit dimensions, Voigt rigidity (`is_rigid`), the
  Kac canonical decomposition, and the Zwara–Bongartz degeneration order with a Hasse
  diagram. *Scope:* the canonical decomposition is Dynkin-only (Euclidean / wild deferred);
  degeneration is representation-finite only; Hall numbers are out of scope.
- **Integration and the v0.2.0 gate (P50).** Every new compute kind is reachable end-to-end
  (GUI canvas → webapp → HPC → worked-steps report) in English and Spanish; the verification
  page is recounted with the `m2` class and all new oracle tables; four new "under the hood"
  chapters (spectral sequences, derived category, τ-tilting, surfaces) join the docs; a
  C1–C8 metagoal scorecard lands in the README; and the curated webapp / desktop examples are
  refreshed.

### Deferred (by choice, to a named successor)

- The σ_A Tamarkin–Tsygan-calculus automorphism (arXiv:2606.15595), the per-HH-degree
  "higher Coxeter polynomials", and τ-Hochschild (co)homology (arXiv:2607.10913) are held to
  a future release (`DEEPER-ENGINES-BACKLOG.md` Tier 2). The classical Coxeter matrix and
  polynomial themselves ship in P38. See the verification page's *v0.2.0 GUI-deferral
  ledger* for the full list.

## [0.1.0] — 2026-08-04

First public release (plans P01–P35, plus Marco's report-completeness passes). Published to
PyPI and GHCR.

### Added

- **Exact computation with finite-dimensional algebras `kQ/I`** over ℂ, GF(p), and GF(pⁿ):
  certified finiteness, Gröbner/Anick bases, the Cartan and Coxeter matrices, global
  dimension, and the centre — exact only, with floats failing loudly by design.
- **Hochschild (co)homology** through four resolutions — the normalized bar complex, the
  minimal `Aᵉ` resolution, the Bardzell resolution (monomial), and the Chouhy–Solotar
  resolution (any admissible presentation, any exact field) — reaching degrees the bar
  complex never can.
- **The Tamarkin–Tsygan product surface** — cup and cap products, the Gerstenhaber bracket,
  and the Connes differential `B` — with cup and cap computed natively past the bar window on
  the Chouhy–Solotar diagonal.
- **Modules and invariants** — simples / projectives / injectives, radical / top / socle,
  Hom / Ext / Tor, minimal projective and injective resolutions, the opposite algebra and
  duality, the Auslander–Reiten translates τ / τ⁻, Krull–Schmidt decomposition, the Yoneda
  Ext-algebra and Koszulity, symmetric / Frobenius / self-injective tests, and the trivial
  extension as a certified quiver presentation.
- **An in-browser no-code GUI (Pyodide)** on the docs landing page — draw a quiver, type
  relations, pick the field, and compute — with live wait estimates; a FastAPI **server
  tier** with a two-tier compute model and a result cache; a **containerized HPC batch tier**
  (`quiverlab-hpc`); and a fully **offline laptop app**.
- **Homework-grade worked-steps reports** (HTML + JSON) with definitions, indexed matrices,
  pivot / rank lines, and justifications for every computation.
- **An optional QPA / GAP backend** (`A.crosscheck(...)`) for independent recomputation, and
  the four-way oracle-class test taxonomy (literature / cross-engine / self-cert / QPA) as
  orthogonal pytest markers.
- Modernized packaging (PEP 639 SPDX license), a documentation site with executable tutorial
  notebooks, a JOSS paper, GitHub Actions CI (OS × Python matrix + engine-path legs), and
  community files.
