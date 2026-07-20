# 12 — Visualization and worked-steps traces

## The mathematics

A quiver `Q` is a directed graph: vertices, and arrows (directed edges, possibly
loops or parallel). An algebra `kQ/I` is that graph plus a list of relations
(parallel-path linear combinations). To *look* at it we draw the graph and print
the relations beneath. To *learn from a computation* — a Hochschild cohomology,
say — we record the resolution chosen, each term, each differential as a matrix
over the working field, each rank, and the resulting dimensions, then lay them out
as a short worked example. Nothing here is new mathematics; it is the same
computation the engine runs, written down so a human can read it.

## How it is represented

**Layout.** `viz.layout.layout(quiver, relations)` returns a `LayoutData`: a dict
`positions` from each vertex to an exact `(x, y)` (`x` an integer column depth, `y`
a `Fraction` row), a tuple of `EdgeRoute`s (straight or parallel, with a `Fraction`
bend), a tuple of `LoopRoute`s (integer base angle), and the relation strings.
Coordinates are `int`/`Fraction` on purpose: they never touch algebra, but keeping
them exact means the float-ban gate covers viz with no exemption, and the layout is
golden-testable to the last coordinate. For the commuting square
`kQ/(a*b - c*d)` the positions are

    1 -> (0, 0)   2 -> (1, 1/2)   3 -> (1, -1/2)   4 -> (2, 0)

**Trace events.** A computation records a flat list of typed events (all plain
dataclasses): a `Dispatch` (which resolution, and why), one `ResolutionTerm` per
degree (its generator count), one `RankStep` per degree (the differential matrix
over the field, or a one-line elision note above 400 cells, plus its rank), and —
for the Chouhy–Solotar engine and the cup/bracket operations — `AmbiguityEvent`,
`DifferentialEvent`, and `LiftStep`. The list is capped at 5000 events so a deep
computation cannot blow up memory.

## How the computation runs

1. **Layering** (`viz.layout.layer`). Compute strongly-connected components
   (iterative Tarjan), condense them, and take the longest path on the condensation
   — so a cycle or loop collapses to one column and every vertex gets an integer
   depth. Vertices in a column are centered on half-integer rows.
2. **Drawing** (`viz.draw.draw_quiver`). Circles at `positions`, `FancyArrowPatch`
   for arrows (parallel bundles fan out via `ConnectionStyle.Arc3(rad=<Fraction>)`),
   `Arc` for loops (integer angle), the relations as a text block below. `tikz`
   emits the *identical* layout as `\node`/`\draw`, so a paper and a screen agree.
3. **Recording** (`trace.recorder.Trace`). `A.hochschild_cohomology(top,
   verbose=…, trace=…)` resolves verbosity (per-call overrides the global
   `quiverlab.verbose`, default `True`), records the engine-choice `Dispatch`, and
   runs the engine with the recorder; the bar engine appends a `ResolutionTerm` and
   a `RankStep` at each degree.
4. **Rendering** (`trace.writer.write_trace`). If `pdflatex` or `tectonic` is on
   `PATH`, compile the LaTeX rendering to `./quiverlab_traces/HHc_<hash>.pdf` and
   print `Worked steps: quiverlab_traces/HHc_<hash>.pdf (N pp)`; otherwise write
   the self-contained, JavaScript-free HTML rendering (math as TeX source) and print
   a loud one-line explanation (no toolchain found, or — if one was found — that
   compilation failed). The resulting
   dimensions in every rendering are *derived from the recorded ranks*
   (`HH^n = dim C^n − rank_n − rank_{n-1}`), so a golden test can assert the
   document's claims equal the engine's own `.dims`.
5. **References.** The engine's `route` maps (via `trace.provenance`) to Plan 06
   REGISTRY keys (e.g. `bar`); those resolve through Plan 06's `bibliography()`
   (`.keys` tuple + entry-view iteration exposing `.key/.formatted/.bibtex_key`)
   into `(bibtex_key, formatted)` lines in the document's **References** section.
   The result itself carries the merged `table.references = self.citations()` (the
   family + engine key union), which Plan 07 does not modify.

## A worked micro-example — HH*(k[x]/(x²)) over ℂ

`A = truncated_polynomial(2, field=CC)` has dimension 2, so the normalized bar
complex has `C^n = 2` for every `n`. Running `A.hochschild_cohomology(2)` records,
per degree, the 2×2 coboundary matrix and its rank:

    d^0 = [[0,0],[0,0]]  rank 0
    d^1 = [[0,0],[2,0]]  rank 1
    d^2 = [[0,0],[0,0]]  rank 0

so `HH^0 = 2−0−0 = 2`, `HH^1 = 2−1−0 = 1`, `HH^2 = 2−0−1 = 1`: dims `[2, 1, 1]`.
Over `GF(2)` the `2` becomes `0`, `d^1` has rank 0, and the dims are `[2, 2, 2]` —
the classic characteristic pathology, visible directly in the traced matrix. These
numbers were produced by running the code (they are the golden trace).

## Where to look in the code

| concept | file | function / class |
|---|---|---|
| layered layout, exact coords | `src/quiverlab/viz/layout.py` | `layout`, `layer`, `LayoutData` |
| matplotlib rendering | `src/quiverlab/viz/draw.py` | `draw_quiver`, `Algebra.draw` |
| TikZ rendering | `src/quiverlab/viz/tikz.py` | `tikz_quiver`, `Algebra.tikz` |
| event taxonomy | `src/quiverlab/trace/events.py` | `Dispatch`, `ResolutionTerm`, `RankStep`, `DifferentialEvent`, `LiftStep`, `AmbiguityEvent`, `ReductionStep` |
| recorder + elision | `src/quiverlab/trace/recorder.py` | `Trace`, `rankstep`, `resolve_verbose` |
| engine emission | `src/quiverlab/hochschild/bar.py` | `hochschild_cohomology_dims` |
| renderers | `src/quiverlab/trace/render_text.py` / `render_latex.py` / `render_html.py` | `render_text`, `render_latex`, `render_html`, `derive_dims` |
| selection + output path | `src/quiverlab/trace/writer.py` | `write_trace`, `have_latex` |
| provenance + References | `src/quiverlab/trace/provenance.py` | `references_for`, `resolve_references` |
