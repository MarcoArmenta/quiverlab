# quiverlab GUI on the docs landing page — Design (Plan 10)

Status: approved by Marco 2026-07-21 (brainstorming session).
Companion to: `2026-07-18-quiverlab-web-design.md` (Plan 09, the server tier —
unchanged, deferred). This document covers the **in-browser** GUI only.

## 1. Goal and audience

The docs landing page (`docs/index.md`, served at `/quiverlab/`) currently
explains how to use the library in code. Its primary audience is algebraists who
do **not** code. The landing page must let them compute examples immediately:
draw a quiver on a canvas, type relations, pick a field, tick invariants, press
Compute, read rendered mathematics, and print a worked-steps PDF — with zero
installation and zero server.

Everything runs in the visitor's browser: the GUI ships the repo's own pure-Python
wheel and executes it under Pyodide (WebAssembly CPython). This is also the local
test bed for the eventual Plan 09 server tier: both tiers speak the same request
schema and (eventually) the same runner.

## 2. Placement and page layout

`docs/index.md` becomes, top to bottom:

1. The one-line pitch (unchanged).
2. **The GUI block** — a raw-HTML container filled by `docs/gui/` assets.
3. "Prefer code?" — the current Python quickstart snippet and the navigation
   bullets (Tutorials / Under the hood / API Reference / Cite). The stale
   "**Web GUI** — … (Plan 09)" bullet is rewritten to point at the GUI block
   above ("or use the GUI at the top of this page").

The GUI shell (canvas, controls, buttons) is plain HTML/SVG/CSS and renders
instantly. The engine (Pyodide runtime + wheel) loads lazily in the background
after page idle; a status chip shows `engine loading…` → `engine ready`. Compute
is disabled until ready. All GUI styling follows the Material theme's light and
dark (slate) palettes via CSS variables.

## 3. Editor UX

- **Add vertex:** click empty canvas → numbered vertex (1, 2, 3, …).
- **Add arrow:** drag vertex → vertex; auto-named `a, b, c, …, z, a1, b1, …`.
  Drag onto the same vertex → loop. Parallel arrows curve apart. Arrowheads
  always visible.
- **Select / edit:** click selects a vertex or arrow; double-click renames an
  arrow (names must be valid library arrow names); `Delete`/`Backspace` removes
  the selection (removing a vertex removes its incident arrows). A Clear button
  empties the canvas.
- **Relations:** one text box, comma-separated, e.g. `a*b - c, d*d`. Validated
  **only** by the library's parser at compute time — never re-implemented in JS.
- **Field picker:** `CC` | `GF(p)` (p input) | `GF(p^n)` (p and n inputs).
  Primality etc. validated by the library (`FieldError` surfaces verbatim).
- **Presets:** a dropdown built from the library's `zoo()` catalog. Selecting a
  preset populates canvas, relations, and field, so a visitor computes an example
  in two clicks. The preset JSON (name → vertices/arrows/relations/field) is
  generated at docs build time by a script that builds each zoo algebra at small
  default parameters and extracts its quiver presentation from the library —
  never hand-maintained. If the library lacks a public accessor the script
  needs, the plan's freshness gate flags it before any GUI code is written.
- **Invariants checklist:** dimension (always computed and shown), `HH^n` and
  `HH_n` with a degree-range picker (UI cap: max degree 10), Cartan matrix,
  Coxeter polynomial, global dimension, center. Plus a "worked-steps report"
  toggle (§5).

## 4. Compute flow

- Pyodide runs in a **Web Worker**; the page never freezes.
- The UI serializes the Plan 09 request schema, `"schema": 1`, with
  `"algebra": {"kind": "quiver", "vertices": […], "arrows": […],
  "relations": […], "field": …}` — the v2 editor shape §5 of the web design
  already reserved. The runner accepts **only** `kind: "quiver"` and rejects
  anything else loudly (`kind: "family"` gets a message pointing at Plan 09);
  presets are expanded to quiver form client-side before submission, so every
  GUI request is a quiver request.
- The worker first **builds** the algebra (finite-dimensionality certified by
  the library; failures surface verbatim), reports the dimension, then computes
  the ticked invariants **one at a time**, posting each result to the page as it
  finishes — partial results render immediately.
- **Cancel** terminates the worker and restarts a fresh one (engine re-init from
  the browser cache is fast; the wheel stays cached).
- No cost estimator in v1: the compute happens on the visitor's own machine;
  Cancel plus the degree cap is the guardrail.

## 5. Results, PDF, artifacts

- Results render into a panel under the GUI: dims tables and polynomials as
  LaTeX typeset by the site's **existing MathJax** setup; each result block
  lists its citations from the library (`A.citations` / result references).
- **Print PDF:** renders the Plan 07 worked-steps **HTML** report (the no-JS
  fallback renderer) in a new tab and triggers the browser print dialog —
  "Save as PDF". No LaTeX toolchain in the browser.
- **Download TikZ:** `A.tikz()` output as a `.tex` file download.
- **Download JSON:** the result JSON (includes the request, the results, and
  `quiverlab.__version__`).
- **Copy Python snippet:** a copy-paste snippet reproducing the computation
  locally (the GUI-to-library bridge, mirroring web design §8).
- **Errors:** `QuiverlabError` subclasses (`RelationError`, `FieldError`,
  `NotFiniteDimensionalError`, `AdmissibilityError`, `DepthLimitError`, …)
  surface verbatim — type name + message — in an error panel. Unexpected errors
  show a generic engine-error message; details go to the browser console.

## 6. Engine payload

- A **mkdocs build hook** builds the pure-Python wheel from the checkout and
  copies it into the site under `gui/` (e.g. `gui/quiverlab-0.1.0.dev0-py3-none-any.whl`).
  The GUI installs it with `micropip` from that **relative URL** — the GUI always
  runs the exact code of the checkout, both under `mkdocs serve`/local builds and
  on GitHub Pages. No PyPI dependency.
- Pyodide itself plus its prebuilt `numpy`, `sympy`, `matplotlib` packages load
  from the **pinned** jsDelivr CDN (exact version chosen at implementation time;
  never a floating tag) — the same supply-chain posture as the site's pinned
  MathJax, with the same documented vendor-it-later option for fully-offline
  serving. First visit downloads ~60 MB once; the browser caches it.
- `numba` (`[fast]` extra) is absent under Pyodide; the engine's pure-Python
  path runs. This path is already CI-tested.

## 7. Where the logic lives (testability)

All request-execution logic lives in **one Python file**, `docs/gui/runner.py`,
shipped verbatim to the browser (mkdocs copies non-markdown docs files into the
site) and imported by the worker:

- `run_build(request_json: str) -> str` — validates the schema, builds the
  algebra into module state, returns basic info (dimension, counts) or a typed
  error.
- `compute_one(spec: str) -> str` — runs one Plan-09 compute string
  (`"hh_cohomology:0..4"`, `"cartan"`, …) against the built algebra; worked-steps
  events accumulate in module state.
- `trace_html() -> str`, `tikz() -> str`, `python_snippet() -> str`,
  `result_bundle() -> str` — the artifacts.

`runner.py` uses the **public** `quiverlab` surface (`import quiverlab`), the
stdlib, and exactly three sanctioned trace helpers for the worked-steps report
string (`quiverlab.trace.render_html.render_html`,
`quiverlab.trace.provenance.references_for`,
`quiverlab.trace.provenance.resolve_references`) — there is no public
HTML-report-as-string API; the freshness gate pins these three and
`quiverlab.engine.*` stays forbidden. It is unit-tested natively with pytest — no browser required —
covering: schema validation and loud rejection of unknown `kind`/fields, every
invariant in the checklist, error passthrough (type name + message), and golden
result JSONs for two known algebras expressed in the quiver schema (e.g. one
loop `x` with relation `x*x*x` over `GF(2)` — the truncated polynomial algebra —
and the Kronecker quiver over `CC`).

The JS stays thin — canvas editing, DOM, worker plumbing — with **no node
toolchain** in the repo. It is verified by a manual browser checklist (documented
in the plan) covering: add/rename/delete vertices and arrows, loops and parallel
arrows, preset load, compute with partial results, cancel, all four
artifact buttons, error display, dark mode.

When Plan 09's server is built, its worker adopts this same runner — one
execution semantics for both tiers.

## 8. Testing

- `tests/gui/test_runner.py` — the pytest suite of §7 (fast bucket).
- A docs-build test asserting the built site contains the GUI assets and exactly
  one wheel under `gui/`, and that `index.md` still passes the existing nav /
  strict-build gates.
- The float-literal AST gate (`tests/test_no_floats.py`) continues to scan
  `src/` only; `docs/gui/runner.py` computes no algebra values itself — it only
  calls the library — but lives outside the gate's scope by construction.
- Manual browser checklist (§7) run before merging; automated browser tests are
  out of scope (no node toolchain).

## 9. Non-goals (v1)

- No Spanish UI (the docs site is English; EN/ES belongs to the Plan 09 app).
- No server tier, job queue, permalinks, feedback form, accounts, or email —
  all remain Plan 09.
- No cost estimator; no mobile-first polish (touch is best-effort).
- No cup products / Gerstenhaber / module Ext / `sweep` in the checklist (the
  curated v1 list of §3 only; more can be added once the GUI exists).

## 10. Plan 09 alignment

- The GUI emits the exact versioned request schema Plan 09 §5 reserved for the
  v2 canvas editor (`kind: "quiver"`); when the server exists, a config switch
  can send big jobs to `POST /api/compute` instead of the local worker.
- Plan 09's plan gets a short amendment noting: the canvas editor exists (this
  GUI), and the server worker should reuse `runner.py`'s execution semantics.
