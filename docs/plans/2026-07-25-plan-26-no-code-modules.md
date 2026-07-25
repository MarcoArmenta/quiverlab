# Plan 26 — no-code module input (GUI + webapp)

**Status:** DELIVERED 2026-07-25, branch `plan-26-no-code-modules` (stacked on
`plan-24-left-modules`, with `plan-25-webapp-result-cache` merged in for its
result cache + canonicalizer). Backlog: DEEPER-ENGINES-BACKLOG.md **Tier 1b item 3**
("No-code module input (GUI + webapp)") — the last unchecked Tier 1b item.

**Goal (Marco's vision):** *every representation theorist can use the tool without
writing code* — specify a module in the GUI (or the API) with zero code and read
off the classical module-level invariants. Plans 23/24 built the engine (`A^op`,
duality `D`, `Tr`, τ/τ⁻, injective resolutions, `side=`, `A.module(...)`); this
plan exposes that engine on both delivered surfaces — the Pyodide landing-page GUI
(`docs/gui/`) and the Plan-09 server tier (`webapp/`).

## Schema v2 — the module block

`webapp/server/schema.py` bumps its versioned request to accept `schema` ∈ {1, 2}.
A `module`/`ext_target` block requires schema 2 (the feature is versioned); every
existing schema-1 family/quiver request is untouched. A `ModuleSpec` is one of two
forms, sharing a top-level `side` (default `"right"`, always emitted):

- **explicit** — `{"dims": {v: n_v}, "maps": {arrow: [[…]]}, "side": …}`: a
  dimension vector + one exact-entry matrix per arrow.
- **builtin pick-list** — `{"builtin": {"kind": "simple"|"projective"|"injective",
  "vertex": v, "side": …}}`: the zero-typing S(v)/P(v)/I(v). A `side` nested inside
  `builtin` (the task's shape) is lifted to the canonical top-level `side`.

New module compute kinds (in `compute`, alongside the algebra kinds):
`dimension_vector`, `rad_top_soc`, `ext:0..n`, `tau`, `tau_minus`,
`projective_resolution:0..n`, `injective_resolution:0..n`,
`projective_dimension`, `injective_dimension`. `ext` also needs a second module
`ext_target` (specifiable the same two ways).

Matrix entries are **exact DATA** — JSON ints or exact strings (`"3"`, `"-2"`,
`"1/2"`); never floats, never evaluated. The schema vets each entry lexically at
request time; the exact parse into the chosen field (GF(p) residue / CC rational)
happens in the runner, where the domain is known, and any residual field
incompatibility surfaces as the library's loud `FieldError`.

### Three JSON examples (schema v2)

Explicit module — dim vector `{1: 2}`, arrow `x` acting nilpotently, τ and a
projective resolution:

```json
{"schema": 2,
 "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
             "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
 "compute": ["dimension_vector", "tau", "projective_resolution:0..4"],
 "module": {"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]}}}
```

Builtin pick-list — the injective `I(2)` as a left module:

```json
{"schema": 2,
 "algebra": {"kind": "quiver", "vertices": [1, 2], "arrows": {"a": [1, 2]},
             "relations": [], "field": {"kind": "GF", "p": 3, "n": 1}},
 "compute": ["rad_top_soc", "injective_dimension"],
 "module": {"builtin": {"kind": "injective", "vertex": 2, "side": "left"}}}
```

Ext with two modules — M explicit, N a builtin simple:

```json
{"schema": 2,
 "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
             "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
 "compute": ["ext:0..4"],
 "module": {"dims": {"1": 2}, "maps": {"x": [[0, 0], [1, 0]]}},
 "ext_target": {"builtin": {"kind": "simple", "vertex": 1}}}
```

## Two runners, one dispatch

The module dispatch lives twice — the two runners already duplicate the algebra
dispatch, and cannot import each other:

- `webapp/server/runner.py` (server tier): `run_spec` builds the module(s) from
  `req.module`/`req.ext_target` and dispatches each module kind. Blocks carry
  `references` (registry step-ids) + resolved `citations`, like the algebra kinds.
- `docs/gui/runner.py` (Pyodide client): `run_build` stores the blocks;
  `compute_one` builds the module lazily and dispatches. GUI blocks carry a MathJax
  `latex` display + `citations`.

**Block → full expansion.** The GUI/API present per-arrow BLOCK matrices sized
`dim[target] × dim[source]` (the natural rep-theory input). The runner expands each
block into the full vertex-ordered `n × n` action matrix `A.module` consumes, using
the **representation quiver's** own arrow directions — `A` for right, `A^op` for
left (a left A-module is a right A^op-module). This is correct for both sides by
construction: on the left the block orientation is transposed, and `A^op`'s
`source`/`target` place it correctly.

**Citations.** `dimension_vector`/`rad_top_soc`/`tau`/`tau_minus` cite ASS2006
(`assem_book`); `ext` cites GSZ2001 (`module_ext`); the resolutions/dimensions cite
GSZ2001 (`minimal_resolution`) (+ ASS2006 for the injective side's duality). The
top-level `references` aggregate these through the bibliography adapter.

## Exactness and loud failure

A module whose matrices violate the relations raises the library's loud
`QuiverlabError` (the base of `RelationError`) from `A.module` /
`from_arrow_action`. In the runner it becomes a typed `RunError`; the server's
`sanitize_error` recognises the QuiverlabError-subclass name as a **safe** type and
returns a clean **422** with the message verbatim — never a 500, never a silent
wrong answer. Shape errors (a block that is not `dim[t]×dim[s]`, an unknown
arrow/vertex, a non-exact entry) are the runner's own `SchemaError` tag — also a
safe 422. The GUI surfaces the same as a NAMED error block (never `InternalError`).
Message localization follows the existing mechanism: library errors surface with
their type + message (as every existing library error does); the surrounding
webapp UI is bilingual, and every NEW user-facing webapp string keeps EN/ES parity
(this plan adds none to the webapp — the no-code FORM lives in the English-only
docs GUI, whose existing labels are hardcoded English).

## Cache composability (Plan 25)

Module requests canonicalize deterministically through `canonical_key`:

- **Side default is explicit** — `ModuleSpec.side` defaults to `"right"` and is
  always emitted, so an omitted `side` and an explicit `"right"` hit the SAME key;
  a builtin's nested `side` is lifted, so both module forms dump one canonical
  shape. The same module typed twice collides; dict key order is irrelevant.
- **Non-module requests are byte-unchanged** — `ComputeRequest.model_dump` drops an
  absent `module`/`ext_target`, so every existing family/quiver request keeps its
  exact pre-Plan-26 cache key (no cache invalidation, no inherited test weakened).
- All three tiers serve module requests; a big-job module request whose key is
  already cached is served with NO email/token/pending row (the Plan-25 crux),
  asserted with the fake mailer in `test_big_module_cache_hit_sends_no_email`.

## Validation limits / tier sizing

`estimator.sizing_dim(algebra_dim, req)` sizes a module request on the larger of
the algebra dimension and the DECLARED module dimension (summed from the spec,
without building the module), so a big module over a small algebra routes off the
instant tier to queued/big — and an absurd module hits the existing beyond-big-cap
reject — exactly like an oversized family request. A builtin pick-list is bounded
by the algebra, so it adds nothing. Wired into `/api/compute` and `/api/jobs/big`.
The instant wall-net and the library's own `max_term_dim` depth guard back the
heuristic.

## GUI (docs/gui/)

A "Module (no code)" fieldset on the existing quiver canvas: an enable toggle, a
build-mode selector (explicit / S(v) / P(v) / I(v)), a right/left side toggle
(right default), a dynamic body that renders per-vertex dimension inputs and
per-arrow exact-entry matrix grids whose dimensions follow the arrow's
source/target dims (transposed for the left side; rebuilt when a dimension is
committed so multi-digit typing never loses focus), and a compute-kind row with
the module kinds (Ext carries a degree select + a builtin target-N pick-list).
`buildRequest` emits the schema `module`/`ext_target` blocks; results render like
existing GUI blocks (MathJax `latex` for scalars, small tables for rad/top/soc,
Ext, and resolutions) with citations. Everything stays client-side (the same
Pyodide engine), no new external assets; the engine loads on first intent
(enabling the panel counts). CSS reuses the theme-aware control styling plus a
compact matrix grid.

## Tests

- `tests/webapp/test_schema.py` — schema v2 parsing, side lifting/collision,
  schema-2 guard, ext-needs-target, float/bool/non-rectangular refusals,
  explicit-vs-builtin exclusivity.
- `tests/webapp/test_runner_modules.py` — every module kind (explicit, builtin,
  multi-vertex block placement, left side), the relation-violation typed refusal,
  shape/unknown-arrow refusals, references aggregation, reproduce-snippet exec.
- `tests/webapp/test_estimator.py` — `sizing_dim` uses the module, ignores
  builtins, is unchanged without a module, routes a big module off instant.
- `tests/webapp/test_cache.py` — module canonicalization units (typed-twice
  collide, side omitted == right, builtin nested side, dict-order, distinctness,
  non-module key unchanged).
- `tests/webapp/test_acceptance.py` — instant module compute, Ext with two
  modules, queued module replayed from cache across users, the big-module cache
  hit with NO email (fake mailer), relation violation as a clean 422.
- `tests/gui/test_runner_modules.py` — the client runner's block goldens, the
  gui.js request shape (string entries) contract, loud error paths, snippet exec.
- `tests/gui/test_interface_freshness.py` — pins the module surface (sided
  builders, `A.module` signature, `Module` methods, `ext_dims`, the
  `assem_book`/`minimal_resolution`/`module_ext` citation step-ids).

The full fast suite, `-m "deep or fast" tests/modules` (Plans 23/24 unregressed),
and the strict docs build (`mkdocs build --strict`, which packs the engine wheel
for Pyodide, exercising the GUI changes) all stay green.

## Verification-page note (merge-time TODO)

`docs/verification.md` does not exist on this branch (it lives on the concurrent
`plan-22-verification-transparency` branch). When the two merge, add Plan 26's
guarantees to that page: module matrices are exact DATA (never evaluated, floats
refused); a relation-violating module is a certified loud refusal surfaced as a
clean 4xx (never a 500, never a silent wrong answer); module invariants carry the
same GSZ2001/ASS2006 provenance as the engine; and module requests canonicalize
through the Plan-25 cache key with the side default explicit.

## Out of scope (unchanged by this plan)

The webapp's own server-rendered pages gain no module FORM (the no-code form is the
docs GUI, per the backlog item); the pre-existing Tier-3 webapp-polish items
(error-envelope unification, `/es` dynamic-label localization, `test_instant_compute`
de-flake) are untouched. Ext's target N in the GUI is a zero-typing S/P/I
pick-list (the API accepts an explicit N too).
