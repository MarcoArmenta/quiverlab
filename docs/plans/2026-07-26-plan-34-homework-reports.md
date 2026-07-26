# Plan 34 — homework-grade reports, honest artifacts, GUI/test gaps (devils mode)

- **Date:** 2026-07-26 · **Branch:** `plan-34-homework-reports` · **Release-gating**
  (Marco: "I don't think it is ready" — v0.1.0 held until this merges).
- **Method:** `/spawn /devils` — 4 workers, each shadowed by a paired
  devil's-advocate critic (read-and-execute access), adjudicated by the
  orchestrator, then a 3-agent fix round + a trace.json build. All critic
  findings below were REPRODUCED, not speculated.

## Marco's feedback (drivers)

1. GUI rad/top/soc: show the matrices at each arrow, not only the dim vector;
   drop the redundant total-dim column.
2. Test Ext and Tor with non-projective, non-injective, non-simple modules.
3. The print report does not give a rendered PDF, and it lacks detail — every
   object must be computed and justified at the level demanded of an
   undergraduate homework in a representation theory course.
4. (Follow-up decisions) The deliverables per computation are **PDF + HTML +
   JSON**, all cached: HTML = the complete human record (no elision,
   scrollable full-width matrices), JSON = the complete machine record (exact
   entries, deterministic, versioned envelope), PDF = the page-bounded
   homework document (scales 11–25-column matrices, states elision past 25
   with pointers to the HTML/JSON). Elision left the recorder: events carry
   full matrices (250 000-cell memory backstop only); renderers decide.

## What shipped (workers + fix round, condensed)

- **rad/top/soc as representations** (W1+fixes): one shared serializer
  `modules/qpa_module.py::module_blocks` → `{dims, maps}` in the Plan-26
  INPUT schema (results feed back in as inputs; left/right round-trip
  certified by `is_isomorphic`), consumed by both runners (structural
  lockstep), rendered with typeset matrices in GUI/webapp/CLI reports; dim
  column gone; GF(pⁿ) entries render via the domain's own notation with a
  `display_only` flag; legacy-shape cache hits show an honest recompute
  notice; request cache keys byte-identical.
- **Ext/Tor non-trivial battery** (W2+F3): both arguments range over
  in-test-certified interior modules (kA₅ intervals — the smallest Aₙ with
  ≥3 of them; kA₄ has exactly one, asserted as the load-bearing
  counterexample; square rad P₁ / P₁∕soc P₁; self-injective uniserials);
  pins: hereditary Euler form ⟨x,y⟩ = xᵀC⁻¹y arbitrated by an ASYMMETRIC
  pair + the independent arrow formula, Tor≡Ext∘D degreewise + resolve-side
  balance + live QPA anchors (square incl.), AR Ext¹(M,τM) ≥ 1 (literature)
  refined == 1 (engine); symmetric `ext_target.side` schema guard + refusal
  tests.
- **Rendered-PDF pipeline** (W3+F1+F2): root cause was amsmath's
  `MaxMatrixCols=10` aborting compiles → silent HTML fallback even with
  LaTeX; now the LaTeX preamble is computed from the events
  (`render_latex.py::matrix_preamble_lines`, shared with `hpc/report.py`),
  a shrink-only `\qlmat` resizebox handles 11–25 columns, >25 states elision
  pointing at HTML/JSON; empty matrices render as `0 (the zero map)`; the
  webapp worker records WHY a PDF is absent (`no_toolchain`/`compile`) and
  pages report the recorded reason (probe = legacy fallback only); no
  artifact is ever labeled PDF unless it is one; GUI print = typeset
  self-contained MathML report + window.print (popup-block message included);
  the MathML converter falls back to escaped source on unknown grammar
  (leak-guard test, proven by mutation).
- **Homework-depth content** (W4+F1): every traced object states the
  definition, shows the matrices (arrow actions, the assembled image matrix
  with pivot summary, the stacked socle kernel, resolution differentials,
  Hom-collapse with `dim = space − rank − rank` spelled out, τ through
  d₁*, Tr = coker, D — and τ⁻ symmetrically), side-aware narration for left
  modules, correct attributions (the rad = MJ justification via A/J
  semisimple, not "Nakayama"), decompose certificates incl. the honest
  char-p refusal; loud drift gates pin the narrated ranks to the engine's
  dims; acceptance samples exercise a genuinely NONZERO Ext differential and
  a dim-12 module (compiled, 0 overfull, read against the bar).
- **trace.json** (final build): schema-versioned, deterministic, exact,
  complete event stream; produced by the writer, promoted by the hpc spec,
  whitelisted/linked in the webapp, downloadable in the GUI,
  `quiverlab-hpc render --format json`.

## The adjudicated critique table (devils round)

| Worker | Critic verdict | Adjudication | Blocking finds (all reproduced, all fixed) |
|---|---|---|---|
| W1 rad/top/soc | NEEDS WORK | accepted | pmatrix content shipped into report.py without the width guard (dim-12 `ReportError`); GF(pⁿ) tuple entries; no display elision policy |
| W2 Ext/Tor | HOLDS UP | accepted + corrections | false "kA₄ has no interior indecomposables" claim; shared-resolution anchor overstated; missing ext_target side guard |
| W3 PDF pipeline | NEEDS WORK | accepted | converter garbled silently (dead fallback); unbounded width lift = off-page PDFs; `quiverlab-hpc render` left broken; split-deployment reason could lie |
| W4 homework depth | REJECTED | accepted with salvage | render_latex lacked the width counter (dim-11–20 dead zone); empty matrices as stray "("; degenerate self-judged samples; left modules mislabeled; Tor mislabeled in HTML |

Lesson recorded: three critics independently converged on centralizing the
matrix-width policy in the renderer preamble computed from events — width
handling was scattered across writer/renderers and each copy failed
differently. Second lesson: acceptance samples must be chosen to EXERCISE
the flagship feature (a nonzero differential), not merely to compile.

## The dispatch amendment (Marco, mid-plan)

Marco hit the stale bar guard ("deeper engines ... arrive in later phases") on
a dim-5 CC algebra. Two fixes, folded into this plan: the hint now tells the
truth (engine="cs" serves any admissible presentation over any exact field),
and `engine="auto"` FALLS BACK to Chouhy–Solotar at exactly the depth where
bar/fast raised DepthLimitError — presented algebras only, recorded in the
dispatch trace, in-window results byte-unchanged, explicit engines keep their
honest walls, presentation-less algebras still refuse. The Plan-04 Pillar-4
pin (`tests/resolutions_cs/test_dispatch.py`) is updated to the amended
contract; regression battery `tests/hochschild/test_auto_cs_fallback_p34.py`.

## Acceptance

1. Suites green (fast/deep/qpa + strict docs) after the trace.json build;
   oracle-class counts resynced (audit gate).
2. The three artifacts produced/served/cached on every surface with honest
   labels; the dim-12 and nonzero-Ext sample PDFs compile clean.
3. Marco's grading pass on the sample reports is the final gate — "homework
   standard" is his call.
