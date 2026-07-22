# quiverlab GUI wait estimates — Design (Plan 11)

Status: approved by Marco 2026-07-22 (brainstorming session; coarse-buckets
option chosen over numeric countdowns and event-streamed progress).
Companion to: `2026-07-21-quiverlab-gui-landing-design.md` (Plan 10, the GUI
this extends).

## 1. Goal

When a visitor configures a computation in the landing-page GUI, show an
honest, coarse estimate of how long the results will take — live while
configuring, and again (with elapsed time) while computing. Estimates are
buckets, never numeric promises: compute cost spans orders of magnitude and
the visitor's machine adds another 2–5× under wasm.

## 2. What the visitor sees

- Once the engine is **ready**, a line under the invariant checklist
  (`#qlgui-eta`) updates as they draw, edit relations, or change field/boxes
  (debounced ~600 ms): `dim = 9 · estimated: under a minute`.
- Buckets (exact copy): `estimated: a few seconds` (< 15 s) ·
  `estimated: under a minute` (< 75 s) · `estimated: a few minutes` (< 6 min) ·
  `estimated: could be long — Cancel anytime` (≥ 6 min) ·
  `will hit the engine's cell cap near degree k` (see §4 — the computation
  would end in the library's own `DepthLimitError`).
- If the configured algebra itself is invalid, the same line shows the
  library error verbatim (type + message) — errors surface **while
  configuring**, before Compute is pressed.
- During a compute, the status chip shows `computing… · 12 s elapsed ·
  estimated: under a minute`; the estimate covers the *remaining* invariants
  and re-scales as each one finishes (§5).
- No confirmation gates for big requests — Cancel stays the only guardrail
  (Plan 10 §4 unchanged). Before the engine loads there is no estimate line:
  Compute is impossible then anyway.

## 3. The live probe is a real build

Each debounced config change sends the current request to the worker, which
runs the existing `run_build` (milliseconds at GUI scale) and then the
estimator. Bonus over a pure formula: the line shows the **true dimension**
and catches `RelationError` / `FieldError` / `NotFiniteDimensionalError`
early. Stated risk: a pathological relation set could make certification
slow during a probe and occupy the worker; Cancel (terminate + reload)
remains the escape hatch. Probe responses carry a sequence number so a stale
response never overwrites a newer line.

## 4. Cost model (in `runner.py`, pure arithmetic)

`estimate(factor) -> JSON` runs against the just-built algebra and the
request's compute list:

- **HH^n / HH_n:** replicate the engine's *exact* bar-basis size arithmetic
  (the same rows × cols numbers the library's `max_cells` guard prints —
  pinned by a golden: Kronecker/CC degree 6 is `8748 × 2916`). Per-degree
  cost `α_route · cells(n)^p`; total = sum over the requested range. Two
  routes, two fitted `α`: prime-field GF(p) rank path and the CC / GF(p^n)
  bar path. **Both α are fitted with numba disabled** — Pyodide has no
  numba, so the browser always runs the pure-Python kernels; a numba-timed
  fit would be wrong by an order of magnitude.
- **Cap prediction:** the same arithmetic decides, before computing, the
  first degree whose matrix exceeds `max_cells = 4_000_000`; if it lies
  inside the requested range the estimate becomes the cap bucket (with that
  degree), since the library will raise `DepthLimitError` there.
- **Cartan / Coxeter / center / gl.dim:** small fitted forms in the
  dimension (constants + a power term where the benchmarks demand one;
  gl.dim's bounded deepen-to-32 gets its own fitted constant).
- The model needs **no new imports** — it is arithmetic on values the runner
  already has (`A.dim`, vertex count, field kind, degree range), keeping the
  Plan-10 import policy intact.

## 5. Calibration and in-flight rescaling

- `calibrate() -> JSON` (runner): time a fixed small workload once at
  engine-ready (~1–2 s) with `time.monotonic()`; the workload's model units
  are computed by the same estimator, so `factor = measured_seconds /
  model_units` self-consistently absorbs machine + wasm speed in one number.
- The runner owns the whole chain `units × factor → seconds → bucket`:
  `estimate(factor)` returns the bucket id + display text, so pytest covers
  the mapping; JS passes the current factor in and displays what comes back.
- While computing, each finished invariant yields a measured
  `seconds / unit` sample; the factor updates (exponential moving average)
  and the *remaining* units are re-bucketed. The elapsed ticker is a 1 s JS
  interval.

## 6. Where the logic lives

All estimation logic — cell arithmetic, route choice, cap prediction, fitted
constants, bucket mapping, calibration workload — lives in `docs/gui/runner.py`
(JSON-string API, pytest-tested natively, shipped verbatim to the browser,
reusable later by the Plan-09 server estimator). `gui.js` adds only: the
`#qlgui-eta` element, the debounce, a worker `probe` command, the ticker,
and the factor bookkeeping. `worker.js` adds the `probe` and `calibrate`
commands and reports per-invariant wall times with each result message.

## 7. Where the constants come from

A committed, non-shipped script (`scripts/fit_eta_model.py`, run manually
like other repo scripts) benchmarks a grid — both HH routes across
dimensions ~2–12 and degrees 0–10 (stopping at the cap), plus the scalar
invariants — **with numba blocked**, and fits the constants by least squares
in log space. The fitted values are baked into `runner.py` as a constants
block with provenance (date, machine, commit). Refitting is rerunning one
script.

## 8. Testing

- Native pytest (`tests/gui/`): the cell-size arithmetic reproduces the
  pinned guard numbers (Kronecker/CC d^6 → 8748 × 2916); cap-degree
  prediction matches an observed `DepthLimitError`; monotonicity (higher
  degree / bigger dim / bar route is never cheaper); bucket mapping goldens
  with a fixed factor incl. the cap bucket; `calibrate()` returns positive
  sane seconds and its workload's declared units.
- Browser suite additions: the eta line appears only after engine-ready;
  shows `dim = 4` + a cheap bucket for the Kronecker preset; worsens (or
  goes to the cap bucket) when HH^0..10 is ticked; shows a verbatim
  `RelationError` inline without pressing Compute; elapsed ticker advances
  during a compute.

## 9. Non-goals

- No numeric countdowns or per-degree progress bars (the trace-event
  streaming option remains open for a later plan).
- No persistence of the calibration factor across visits.
- No request blocking or confirmation gates.
- No server-side estimation (Plan 09 may reuse `estimate` as-is later).
