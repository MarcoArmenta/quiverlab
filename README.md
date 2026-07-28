# quiverlab

[![CI](https://github.com/MarcoArmenta/quiverlab/actions/workflows/ci.yml/badge.svg)](https://github.com/MarcoArmenta/quiverlab/actions/workflows/ci.yml)
[![Docs](https://github.com/MarcoArmenta/quiverlab/actions/workflows/docs.yml/badge.svg)](https://marcoarmenta.github.io/quiverlab/)
[![Tests](https://img.shields.io/badge/tests-2295_oracle--pinned-brightgreen)](https://marcoarmenta.github.io/quiverlab/verification/)
[![PyPI](https://img.shields.io/pypi/v/quiverlab.svg)](https://pypi.org/project/quiverlab/)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue)](https://github.com/MarcoArmenta/quiverlab/blob/main/pyproject.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
<!-- [![DOI](https://joss.theoj.org/papers/<id>/status.svg)](https://doi.org/<doi>) -- added at JOSS acceptance -->

**Exact representation theory of quivers with relations, for algebraists** —
modules and Auslander–Reiten theory, resolutions, Ext-algebras and Koszulity,
Hochschild (co)homology with its Gerstenhaber calculus, cyclic homology, and
Cartan/Coxeter/spectral invariants, all exactly.

## The two metagoals

quiverlab is built toward two long-term goals, and every release is measured
against them:

1. **No code required.** Every computation the library can do should be
   reachable without writing a single line of code: draw the quiver and the
   relations in the browser GUI, specify modules entry-by-entry in the no-code
   panel, export a config file for a cluster, and read the results as rendered
   mathematics (or a PDF report). Python is a power-user option, never a
   prerequisite.
2. **Any computation done in representation theory.** The aim is that whatever
   a representation theorist of finite-dimensional algebras computes in a paper
   — homological invariants, module-theoretic constructions, Auslander–Reiten
   data, Ext algebras, spectral/Coxeter data, and beyond — can be computed
   here, exactly and with certified, oracle-tested results. The gap between
   this goal and the current surface is tracked openly as the coverage program
   in [`docs/plans/ROADMAP.md`](docs/plans/ROADMAP.md); if your computation is
   missing, it belongs on that list.

quiverlab computes with finite-dimensional algebras `kQ/I` over the complex numbers
(exactly — no floating point, ever) and over all finite fields: certified
finite-dimensionality, Hochschild (co)homology with cup products and Gerstenhaber
brackets, the first full Chouhy–Solotar resolution, module Ext, and Cartan/Coxeter
invariants. Floats fail loudly by design.

## Run the GUI

**Zero install — it runs in your browser:** open
<https://marcoarmenta.github.io/quiverlab/> and use the canvas at the top of the
page. Draw vertices and arrows, add relations, pick a field, and compute Hochschild
(co)homology exactly; the full quiverlab engine runs client-side via Pyodide (the
first computation loads it — give it a few seconds).

To launch the same GUI locally from a clone (it always runs the exact code you
checked out):

```bash
git clone https://github.com/MarcoArmenta/quiverlab.git && cd quiverlab
python -m venv .venv && source .venv/bin/activate
pip install -e ".[docs]"
mkdocs serve
```

then open <http://127.0.0.1:8000> — the GUI is at the top of the landing page.
The first start packs the engine wheel and executes the tutorial notebooks; give
it a minute or two.

## Install

```bash
pip install quiverlab                 # pure-Python core, no external systems
pip install "quiverlab[fast]"         # + numba GF(p) acceleration (optional)
pip install "quiverlab[qpa]"          # + GAP/QPA cross-check backend (macOS/Linux)
```

## Three lines to a Hochschild table

```python
from quiverlab import Quiver, CC

Q = Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3), "c": (1, 3)})
print(Q.algebra(relations=["a*b"], field=CC).hochschild_cohomology(3))
```

## Learn more

- **Documentation:** <https://marcoarmenta.github.io/quiverlab/>
- **Tutorials:** [executable notebooks](docs/tutorials/) — start here.
- **Under the hood:** [internals chapters](docs/internals/) — how each number is produced.
- **Web GUI:** the [docs landing page](https://marcoarmenta.github.io/quiverlab/) computes in your browser today; a self-hostable server tier (`webapp/`) adds queued and email-verified big jobs — see [Web interface](#web-interface).
- **Cite:** see the JOSS paper (`paper/paper.md`) and [`CITATION.cff`](CITATION.cff).

## The classic characteristic pathology, in one loop

```python
from quiverlab import truncated_polynomial, CC, GF

for field in (CC, GF(2), GF(3)):
    print(field, truncated_polynomial(2, field=field).hochschild_cohomology(4).dims)
# CC     [2, 1, 1, 1, 1]
# GF(2)  [2, 2, 2, 2, 2]
# GF(3)  [2, 1, 1, 1, 1]
```

## General quivers with relations (kQ/I)

```python
from quiverlab import Quiver, CC

Q = Quiver(vertices=[1, 2, 3, 4],
           arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
A = Q.algebra(relations=["a*b - c*d"], field=CC)   # commutative square, exact
print(A.dim)                                        # 9
print(A.hochschild_cohomology(1))                   # HH^0 = 1  HH^1 = 0
```

Non-monomial relations are completed with an exact noncommutative Gröbner
(Buchberger–Mora overlap) engine and certified finite-dimensional; a
non-admissible or infinite presentation fails loudly with `AdmissibilityError`
or `NotFiniteDimensionalError`, never a hang.

## Modules and invariants

```python
from quiverlab import Quiver, CC

A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)}
           ).algebra(relations=["a*b - c*d"], field=CC)   # commutative square

S1, S4 = A.simple(1), A.simple(4)
A.projective(1).dimension_vector()      # {1: 1, 2: 1, 3: 1, 4: 1}
A.ext(S1, S4, 2)                        # 1     (Ext^2 of simples)
int(A.global_dimension())              # 2
A.loewy_length()                       # 3
A.simple(1).projective_resolution(4)   # P_1 <- P_2(+)P_3 <- P_4 <- 0
```

Every module is a right A-module over the stated exact field; Ext, Hom, and the
projective resolution are exact. Exact `spectral_radius`/`mahler_measure`, `center()`,
`complexity()` (a lower-bound estimate — can under-report, exact only on local /
single-vertex inputs), and `sweep()` (invariant × field) round out the invariant surface.

## Families and citations

```python
from quiverlab import NakayamaAlgebra, QuantumCI, families, bibliography

A = NakayamaAlgebra([3, 2, 2])          # cyclic Nakayama, dim 7
print(A.hochschild_cohomology(0))       # HH^0 = 1
print(A.citations())                    # ('nakayama', 'assem_book', 'bar')

print(families())                       # the whole v1 catalog with signatures
print(bibliography(A.citations()))      # grouped, annotated references
```

## How quiverlab is verified

Every shipped feature is unit tested (the suite is 2295 tests over the
`[dev,fast,docs,web,qpa,hpc]` extras), and the mathematics is pinned by **two classes
of oracle** — surfaced since Plan 32 as four orthogonal, runnable marker classes
(`oracle_literature` / `oracle_crossengine` / `oracle_selfcert` / `qpa`), audited
against live collection:

- **Theory and literature, on constructed examples.** We build many algebras the
  literature (or a theorem we know) has already resolved and assert quiverlab
  reproduces the published value exactly — Happel's hereditary vanishing, the
  Buchweitz–Green–Madsen–Solberg / Bergh–Erdmann quantum complete intersection,
  the classical `k[x]/(x^n)` and Künneth commutative-CI values, and more. Where no
  single published vector is at hand we cross-check an *independent* path in the
  library and say so inline. The read-only hanlab bank supplies byte-level
  closed-form oracles.
- **Cross-engine and external agreement.** The bar complex, the minimal `A^e`,
  Bardzell, and Chouhy–Solotar resolutions are independent engines; where two
  overlap they must agree degreewise over the primes `{32003, 2, 3, 5}`. And
  wherever the GAP package **QPA** implements a feature we recompute with it and
  demand equality (`A.crosscheck(...)`). QPA does not implement everything
  quiverlab does; the docs page names exactly where it is used and which theory
  oracle stands in where it cannot.

Exactness is enforced structurally: an AST gate bans every float from `src/`, and
the entire deep suite runs twice in CI — once on the numba kernels, once on the
pure-Python path (`QUIVERLAB_NO_NUMBA=1`) — with the two required to agree exactly.
The full methodology, a subsystem → oracles → test-file table, the CI matrix, and
an honest-scope section live in **[How quiverlab is verified](https://marcoarmenta.github.io/quiverlab/verification/)**
(`docs/verification.md`).

## Status

Engine, module, and families phase (Plans 01–06 delivered, together with the
Plan-04 Chouhy–Solotar resolution). On top of the foundations — monomial presentations,
exact fields, bar-complex Hochschild (co)homology — the hanlab deep engine is now
ported and wired in:

- **A fast GF(p) engine** behind the field interface: `hochschild_cohomology`
  and `hochschild_homology` take `engine="auto" | "bar" | "fast"`. `auto` picks
  the numpy mod-p rank engine over prime fields and the exact bar path
  everywhere else; both agree exactly where both can run. The fast engine still
  builds the exponential bar basis, so it guards its depth loudly (raise
  `max_cells` deliberately) — the depth *unlock* lives in the resolutions below.
- **Deep monomial resolutions.** The minimal (Bardzell) and periodic bimodule
  resolutions reach degrees the bar complex never could — k[x]/(x^a) and cyclic
  Nakayama to depth 40 instantly — and certify structural facts (a finite global
  dimension shows up as vanishing generators), cross-checked exactly against the
  bar oracle over primes {32003, 2, 3, 5} on the overlap range.
- **The Chouhy–Solotar resolution** (`resolutions_cs`, `engine="cs"`). The
  domain-generic CS projective bimodule resolution for admissible kQ/I — its
  HH•/HH^• dimensions and representative (co)cycles reach Hochschild degrees the
  bar oracle cannot, with CS↔bar comparison maps; it specializes to Bardzell's
  minimal resolution on monomial algebras (operation transport is certified
  inside the bar-buildable window).
- **Tamarkin–Tsygan calculus** at the engine level: cup product, cap product,
  and the Gerstenhaber bracket; plus **cyclic homology** (Connes' mixed complex).
- **Invariants:** the integer **Cartan** matrix, the **Coxeter** matrix and its
  characteristic polynomial (all fields, exact via sympy); and, over GF(p), the
  **Nakayama** automorphism with the **Frobenius** and **symmetric** tests
  (loud `FieldError` off a prime field).
- **Modules, scalar invariants, and the exact spectral layer.** Right A-modules
  with exact **Ext**, **Hom**, and minimal **projective resolutions**; the scalar
  invariants **Loewy length**, **center**, and **complexity** (GF(p); the last a
  lower-bound estimate that can under-report, exact only on local / single-vertex
  inputs); and the
  exact **spectral radius** / **Mahler measure** of the Coxeter polynomial as
  sympy algebraic numbers — no floats, ever.
- **Algebra families and citations.** A curated catalog of named families
  (`NakayamaAlgebra`, `QuantumCI`, `ExteriorAlgebra`, `IncidenceAlgebra`,
  `PreprojectiveAlgebra`, `TrivialExtension`, `TensorProduct`, …) with `families()`
  discovery and the `zoo` iterator, each stamped with the literature it comes from;
  `A.citations()` and `bibliography(...)` resolve those keys to grouped, annotated
  references, plus a batch scan surface for family sweeps.
- **In-browser GUI** — the [docs landing page](https://marcoarmenta.github.io/quiverlab/)
  computes examples with nothing installed (Pyodide running the same exact engine).

Everything is exact — no floating point, ever — and the full test suite runs
green on both the numba kernel path and the pure-Python path
(`QUIVERLAB_NO_NUMBA=1`).

Honest scope note: the calculus lives at the *engine* level today. A classy
`A.cup(u, v)` on named cohomology classes awaits the cohomology-classes
machinery of a later phase (see `docs/plans/ROADMAP.md`).

Coming next (see `docs/plans/ROADMAP.md`): full operation transport, drawing and
TikZ export, worked-steps PDFs, and an optional QPA backend.

## Draw it, and read the worked steps

```python
from quiverlab import Quiver, CC

Q = Quiver(vertices=[1, 2, 3, 4],
           arrows={"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
A = Q.algebra(relations=["a*b - c*d"], field=CC)

A.draw(file="square.svg")     # matplotlib PNG/SVG: loops, parallels, relations below
print(A.tikz())               # same layout, paste-into-paper TikZ

A.hochschild_cohomology(2)    # writes quiverlab_traces/HHc_<hash>.pdf (or .html) and
                              # prints: Worked steps: quiverlab_traces/HHc_3f2a.pdf (N pp)
```

Worked-steps documents are on by default (`quiverlab.verbose = True`); every claim
in them is a golden-file-tested equality with the value the engine computed. Turn
them off per call (`A.hochschild_cohomology(2, verbose=False)`) or globally
(`quiverlab.verbose = False`). Reports are delivered as a self-contained,
JavaScript-free HTML document (math shown as TeX source) plus an exact JSON event
stream; the browser's Print-to-PDF turns the HTML into a page-ready document when
one is needed.

## Web interface

A no-code web GUI (`webapp/`) exposes the library for algebraists who prefer not
to write Python: pick a family, a field, and invariants; read exact results with
rendered mathematics; download the worked-steps PDF. Small computations run
instantly; deep ones become queued jobs with a permalink; very large ones run as
email-verified **big jobs** (a single-use magic link; requires an outbound SMTP
relay, disabled otherwise). Every result carries a References block (the
literature the computation stands on, from the library's citations subsystem),
and `/literature` shows the full curated bibliography. The UI is bilingual
(English at `/`, Spanish at `/es/`) with a public feedback form at `/feedback`
(including a "suggest literature" category).

Results are cached: because every computation is exact and deterministic, a
previously computed example is never recomputed — an identical request is served
instantly from the cache, across users. Email verification gates only the *cost* of
computing a new big example, not access to the mathematics, so a big example that
someone already computed is served immediately, with no email needed.

Each finished computation exposes downloadable artifacts under
`/download/<job-id>/…`: `result.json` (exact dimensions, references, and a
copy-paste reproduction snippet), the worked-steps `trace.pdf` (or a
self-contained `trace_steps.html` when no LaTeX toolchain is present), and
`tikz.tex` when a drawing was requested. Every number is exact — the server never
approximates, and an out-of-scope request fails loudly rather than silently
truncating.

Run it locally:

```bash
pip install -e ".[web,fast]"
uvicorn webapp.server.app:create_app --factory --reload      # terminal 1
python -m webapp.worker.run_loop                             # terminal 2
# open http://127.0.0.1:8000
```

The web tier is two processes sharing one SQLite database: the FastAPI app
(instant computations under a hard wall-time net; everything larger is enqueued)
and one or more worker loops (each job runs in a resource-capped subprocess). A
full-stack local smoke driving the real processes over HTTP lives at
[`scripts/webapp_smoke.py`](scripts/webapp_smoke.py); the equivalent flow runs
in-process (no ports) as `tests/webapp/test_acceptance.py`.

Deploy (DRAC Arbutus, Docker Compose + Caddy TLS): see
[`webapp/deploy/PROVISIONING.md`](webapp/deploy/PROVISIONING.md).

## HPC and offline use (container)

The same library ships as **one container** (`ghcr.io/marcoarmenta/quiverlab`) with
a `quiverlab-hpc` CLI, serving two stories from the one image.

**Run a big example on a SLURM cluster in 5 steps** (only `ssh`/`scp`/`sbatch`
needed; Apptainer is rootless):

```bash
apptainer pull quiverlab.sif docker://ghcr.io/marcoarmenta/quiverlab:latest   # 1. pull
apptainer run quiverlab.sif sample-config > my-config.yaml                    # 2. config (or export from the GUI)
sbatch slurm/quiverlab-drac.sbatch my-config.yaml result.json                 # 3. submit
scp you@cluster:result.json .                                                 # 4. fetch
apptainer run --bind "$PWD" quiverlab.sif render result.json -o report.html   # 5. render locally (HTML/JSON)
```

Very large examples become reachable via **atomic per-degree checkpoints**: a job
that runs out of wall time exits 75, requeues, and resumes from `$SCRATCH` on the
next submit — just `sbatch` again. **quiverlab is CPU-only** — request cores
(`--cpus-per-task`) and RAM (`--mem`), never a GPU; the arithmetic is exact
(integers mod p / rationals) and a GPU would sit idle. `quiverlab-hpc estimate
my-config.yaml` suggests the resources.

**Offline laptop app.** Pull the image once with internet, then run
`apptainer run quiverlab.sif gui` (or `docker run -p 8000:8000 … gui`) and open
`http://localhost:8000` — the zero-code GUI computes locally with no network, showing
your machine's detected cores/RAM, memory/time estimates, and the limits you are
computing under, and ships precomputed examples.

Full instructions: [Run on your HPC cluster](docs/hpc.md) and
[Offline laptop app](docs/offline-app.md).

MIT © 2026 Marco Armenta
