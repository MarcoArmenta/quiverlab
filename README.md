# QuiverLab

[![CI](https://github.com/MarcoArmenta/quiverlab/actions/workflows/ci.yml/badge.svg)](https://github.com/MarcoArmenta/quiverlab/actions/workflows/ci.yml)
[![Docs](https://github.com/MarcoArmenta/quiverlab/actions/workflows/docs.yml/badge.svg)](https://marcoarmenta.github.io/quiverlab/)
[![Tests](https://img.shields.io/badge/tests-2981_oracle--pinned-brightgreen)](https://marcoarmenta.github.io/quiverlab/verification/)
[![PyPI](https://img.shields.io/pypi/v/quiverlab.svg)](https://pypi.org/project/quiverlab/)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue)](https://github.com/MarcoArmenta/quiverlab/blob/main/pyproject.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
<!-- [![DOI](https://joss.theoj.org/papers/<id>/status.svg)](https://doi.org/<doi>) -->

# **Exact representation theory of quivers with relations, for algebraists** 
Modules and Auslander–Reiten theory, resolutions, Ext-algebras and Koszulity,
Hochschild (co)homology with its calculus (cup, Gerstenhaber, cap, Connes), cyclic homology, and
Cartan/Coxeter/spectral invariants, all exactly.

## ⬇️ DOWNLOAD APPLICATION HERE

> ### **[⬇ Download the QuiverLab app](https://github.com/MarcoArmenta/quiverlab/releases/tag/app-latest)** — one file, no install, no code.
>
> Double-click it and the GUI opens in your browser: draw a quiver, pick a
> field, read exact results. Fully offline.
>
> | OS | Download |
> |---|---|
> | **macOS** (Apple Silicon: M1–M4) | [QuiverLab-macos-arm64.zip](https://github.com/MarcoArmenta/quiverlab/releases/download/app-latest/QuiverLab-macos-arm64.zip) |
> | **Windows** | [QuiverLab-windows.exe](https://github.com/MarcoArmenta/quiverlab/releases/download/app-latest/QuiverLab-windows.exe) |
> | **Linux** (x86-64) | [QuiverLab-linux-x86_64.tar.gz](https://github.com/MarcoArmenta/quiverlab/releases/download/app-latest/QuiverLab-linux-x86_64.tar.gz) |
>
> First-open notes (the app is not code-signed yet, so each OS warns **once**):
> 
> **macOS** — unzip and double-click; when the *"Apple could not verify…"*
> dialog appears click **Done** (not "Move to Trash"), then System Settings →
> Privacy & Security → scroll to Security → **Open Anyway** → Open. (Terminal
> alternative: `xattr -d com.apple.quarantine ./QuiverLab`.)
>
> **Windows** — if
> SmartScreen appears, choose *More info* → *Run anyway*.
>
> **Linux** — `tar xzf`,
> then run `./QuiverLab`.
>
> **Intel Mac** — no one-file build (GitHub retired its
> Intel-mac builders); use
> `docker run -p 8000:8000 ghcr.io/marcoarmenta/quiverlab:latest gui`
> or the [pip path](https://marcoarmenta.github.io/quiverlab/offline-app/).


## The two metagoals

QuiverLab is built toward two long-term goals, and every release is measured
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

QuiverLab computes with finite-dimensional algebras `kQ/I` over the complex numbers
(exactly — no floating point, ever) and over all finite fields: certified
finite-dimensionality, Hochschild (co)homology with cup products and Gerstenhaber
brackets, the first full Chouhy–Solotar resolution, module Ext, and Cartan/Coxeter
invariants. Floats fail loudly by design.

## Get QuiverLab

Most users want one of these, in this order:

**1. Download the desktop app** — one file, double-click it, and the zero-code
GUI opens in your browser on localhost, fully offline, using your machine's
real cores and RAM. Grab the binary for your OS from the
[**download box at the top of this page**](#️-download-application-here)
(macOS / Windows / Linux), which also carries the one-time first-open steps for
the unsigned binaries.

**2. Download the containerized application** — one image, the full exact
engine, no Python setup (registry paths are lowercase-only):

```bash
docker pull ghcr.io/marcoarmenta/quiverlab:latest      # or: apptainer pull quiverlab.sif docker://ghcr.io/marcoarmenta/quiverlab:latest
docker run --rm -p 8000:8000 ghcr.io/marcoarmenta/quiverlab:latest gui
# open http://localhost:8000 — the zero-code GUI, fully offline, using your
# machine's cores and RAM. The same image runs batch configs; see
# "Writing and running config files" below.
```

**3. Clone the repo and build the container yourself:**

```bash
git clone https://github.com/MarcoArmenta/quiverlab.git && cd quiverlab
docker build -f container/Dockerfile -t quiverlab:local .
docker run --rm -p 8000:8000 quiverlab:local gui
```

**4. Use the web interface** — the self-hostable server tier (`webapp/`):
instant answers for small examples, queued jobs with permalinks for deep ones,
and a shared exact-result cache — see [Web interface](#web-interface).

**5. Prefer code?** - Python-library installs and SLURM clusters are covered
[at the bottom](#install-the-python-library).

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
- **No-code interfaces:** the containerized app ships an offline GUI (`quiverlab-hpc gui`), and the self-hostable server tier (`webapp/`) adds queued and email-verified big jobs — see [Web interface](#web-interface).
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

Every shipped feature is unit tested (the suite is 2981 tests over the
`[dev,fast,docs,web,qpa,hpc]` extras), and the mathematics is pinned by **two classes
of oracle** — surfaced since Plan 32 as five orthogonal, runnable marker classes
(`oracle_literature` / `oracle_crossengine` / `oracle_selfcert` / `qpa` / `m2`), audited
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
  oracle stands in where it cannot. And a live **Macaulay2** bridge recomputes nc
  graded dimensions and commutative Ext data (single-vertex scope; `-m m2`) — a
  genuinely different computer-algebra system as a second external oracle.

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
- **Tamarkin–Tsygan calculus**, as a public product surface: **cup/cap products,
  the Gerstenhaber bracket, and the induced Connes differentials** on `HH^•`/`HH_•`
  (`A.cup_products`, `A.cap_products`, `A.gerstenhaber_brackets`,
  `A.connes_differentials`) — exact structure-constant tables on the recorded HH
  basis, with worked-steps reports; plus **cyclic homology** (Connes' mixed complex).
- **Invariants:** the integer **Cartan** matrix, the **Coxeter** matrix and its
  characteristic polynomial (all fields, exact via sympy); **Euler / Tits forms**
  with exact finite/tame/wild definiteness, orientation-blind **Dynkin/Euclidean
  type detection**, **positive-root** enumeration for Dynkin type, and the
  **structural recognizers** (`is_semisimple` … `is_gentle`, with a live QPA
  crosscheck); **Koszulity** and the Yoneda Ext-algebra clickable in the no-code
  GUI; and, over GF(p), the **Nakayama** automorphism with the **Frobenius** and
  **symmetric** tests (loud `FieldError` off a prime field).
- **Modules, scalar invariants, and the exact spectral layer.** Right A-modules
  with exact **Ext**, **Hom**, and minimal **projective resolutions**; the scalar
  invariants **Loewy length**, **center**, and **complexity** (GF(p); the last a
  lower-bound estimate that can under-report, exact only on local / single-vertex
  inputs); and the
  exact **spectral radius** / **Mahler measure** of the Coxeter polynomial as
  sympy algebraic numbers — no floats, ever.
- **Homological dimensions (C6).** Public **syzygy/cosyzygy** operators,
  **finitistic / dominant / Gorenstein dimensions**, the **Igusa–Todorov φ/ψ**
  functions, and **Ω/τ-periodicity certificates** — the C6 homological-dimensions
  family, each result carrying the `GlobalDimension`-style certified-value-or-honest-bound
  honesty (never a bare number when unresolved, `is_gorenstein` three-valued
  True/None), and clickable end-to-end via the no-code `homological_profile`.
- **Algebra families and citations.** A curated catalog of named families
  (`NakayamaAlgebra`, `QuantumCI`, `ExteriorAlgebra`, `IncidenceAlgebra`,
  `PreprojectiveAlgebra`, `TrivialExtension`, `TensorProduct`, …) with `families()`
  discovery and the `zoo` iterator, each stamped with the literature it comes from;
  `A.citations()` and `bibliography(...)` resolve those keys to grouped, annotated
  references, plus a batch scan surface for family sweeps.
- **Zero-code GUI** — the containerized app serves the full-engine GUI offline
  on localhost (`quiverlab-hpc gui`), with your machine's real cores and RAM.

Everything is exact — no floating point, ever — and the full test suite runs
green on both the numba kernel path and the pure-Python path
(`QUIVERLAB_NO_NUMBA=1`).

Honest scope note: the calculus is now public as **structure-constant tables** over
the whole HH basis (`A.cup_products(top)` and friends). A classy `A.cup(u, v)` on
two *named* cohomology-class representatives still awaits the cohomology-classes
machinery of a later phase (see `docs/plans/ROADMAP.md`); and the Gerstenhaber
bracket is GF(p)-only and window-bounded.

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


## Writing and running config files (the containerized app)

Everything the container computes is driven by **one YAML document** — the same
schema the webapp and the browser GUI speak, so a config exported from the GUI
runs unchanged on a cluster. Run it with any of the three installs:

```bash
# Docker (make the output dir writable for the in-image uid first)
mkdir -p out && chmod 777 out
docker run --rm -v "$PWD:/cfg:ro" -v "$PWD/out:/out" quiverlab:local \
    run /cfg/my-config.yaml -o /out/result.json
docker run --rm -v "$PWD/out:/out" quiverlab:local \
    render /out/result.json -o /out/report.html --format html

# Apptainer (clusters; rootless)
apptainer run --bind "$PWD" quiverlab.sif run my-config.yaml -o result.json
apptainer run --bind "$PWD" quiverlab.sif render result.json -o report.html

# Plain pip install (no container)
pip install "quiverlab[fast,hpc]"
quiverlab-hpc run my-config.yaml -o result.json
quiverlab-hpc render result.json -o report.html
```

`quiverlab-hpc sample-config` prints an annotated template and
`quiverlab-hpc estimate my-config.yaml` suggests `--time/--cpus-per-task/--mem`
before you submit. The rendered report shows the quiver presentation (labeled
arrows), every requested invariant with rendered matrices, and a resources
footer (wall time, peak RSS, cores).

### Anatomy of a config

```yaml
schema: 2                  # 1 = algebra-only; 2 required for module blocks
algebra:                   # EITHER a named family ...
  kind: family
  family: QuantumCI        # discover names: python -c "import quiverlab; print(quiverlab.families())"
  params: {q: 2, a: 2, b: 2}
  field: {kind: GF, p: 32003, n: 1}    # GF(p^n), or {kind: CC} for exact char 0
compute:                   # any subset; ranged kinds take "kind:lo..hi"
  - "hh_cohomology:0..8"
  - cartan
artifacts: {tikz: true}    # optional; tikz.tex written beside result.json
hpc:                       # optional; CLI-only budgets
  time_limit_s: 3600
  max_mem_bytes: 4294967296
```

**Compute kinds.** Algebra-level: `hh_cohomology:lo..hi`, `hh_homology:lo..hi`,
`cartan`, `coxeter_polynomial`, `global_dimension`, `center`, `dimension`.
Module-level (need a `module` block, schema 2): `dimension_vector`,
`rad_top_soc`, `decompose`, `tau`, `tau_minus`, `projective_resolution:0..n`,
`injective_resolution:0..n`, `projective_dimension`, `injective_dimension`,
`ext:0..n` (needs `ext_target`), `tor:0..n` (needs `tor_target`, a **left**
module).

**Module blocks.** A module is either a builtin pick
(`module: {builtin: {kind: simple|projective|injective, vertex: 3, side: right}}`)
or an explicit representation: `dims` maps **string** vertex labels to
dimensions (missing vertices are 0), `maps` gives one `dim_target x dim_source`
matrix per arrow (arrows touching a 0-dimensional vertex may be omitted).
Entries are exact data — integers or fraction strings like `"1/2"`; floats are
refused loudly. `side: left` means a representation of the opposite quiver.

### Worked configs

A hereditary path algebra over **exact characteristic 0** — no proxy prime:

```yaml
schema: 1
algebra:
  kind: family
  family: PathAlgebra
  params: {type_or_quiver: "A5"}
  field: {kind: CC}
compute: [cartan, coxeter_polynomial, global_dimension, dimension]
# dim 15, gl.dim = 1 (exact), the A5 Coxeter polynomial
```

The exterior algebra in char 0 — Hochschild cohomology grows linearly:

```yaml
schema: 1
algebra:
  kind: family
  family: ExteriorAlgebra
  params: {n: 2}
  field: {kind: CC}
compute: ["hh_cohomology:0..4", center, dimension]
# HH^0..4 = [2, 4, 6, 8, 10]
```

A truncated path algebra over the **non-prime field GF(9)**:

```yaml
schema: 1
algebra:
  kind: family
  family: TruncatedPathAlgebra
  params: {type_or_quiver: "A6", r: 3}
  field: {kind: GF, p: 3, n: 2}
compute: [cartan, global_dimension, "hh_cohomology:0..4"]
# gl.dim = 3 (exact)
```

An **explicit quiver** (the Kronecker quiver, no relations) with a no-code
module given by matrices — the regular representation `R_2` (`a` acts by 1,
`b` by 2):

```yaml
schema: 2
algebra:
  kind: quiver
  vertices: [1, 2]
  arrows: {a: [1, 2], b: [1, 2]}
  relations: []
  field: {kind: GF, p: 5, n: 1}
compute: [dimension, cartan, global_dimension, dimension_vector,
          rad_top_soc, decompose, tau, "projective_resolution:0..3"]
module:
  side: right
  dims: {"1": 1, "2": 1}
  maps:
    a: [[1]]
    b: [[2]]
```

An explicit quiver with a **non-monomial relation** — the commutative square,
over CC:

```yaml
schema: 1
algebra:
  kind: quiver
  vertices: [1, 2, 3, 4]
  arrows: {a: [1, 2], b: [1, 3], c: [2, 4], d: [3, 4]}
  relations: ["a*c - b*d"]
  field: {kind: CC}
compute: [dimension, global_dimension, center, "hh_cohomology:0..3"]
# dim 9, gl.dim = 2 (exact)
```

Larger ready-to-run configs live in
[`container/examples/`](container/examples/): the quantum complete intersection
with the full invariant surface ([`qci-q2.yaml`](container/examples/qci-q2.yaml)),
a cyclic Nakayama algebra with a decomposable module
([`nakayama-kz4.yaml`](container/examples/nakayama-kz4.yaml)), the 3x3
commutative grid with interior modules paired by the Auslander-Reiten
translate — `Ext^1(M, tau M) = 1` ([`grid3x3.yaml`](container/examples/grid3x3.yaml)),
and a dim-220 deep-degree run ([`nakayama-kz20-deep.yaml`](container/examples/nakayama-kz20-deep.yaml)).
Every one computes byte-identically in the container and from the wheel.

## Install the Python library

```bash
pip install quiverlab                 # pure-Python core, no external systems
pip install "quiverlab[fast]"         # + numba GF(p) acceleration (optional)
pip install "quiverlab[qpa]"          # + GAP/QPA cross-check backend (macOS/Linux)
pip install "quiverlab[fast,hpc]"     # + the quiverlab-hpc CLI (configs, reports)
```

MIT © 2026 Marco Armenta
