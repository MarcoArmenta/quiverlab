# Plan 28 — quiverlab as a container: HPC batch tier + offline laptop app

**Date:** 2026-07-25. **Branch:** `plan-28-hpc-container`. **Backlog:** new item
(Marco, 2026-07-25) — added to `DEEPER-ENGINES-BACKLOG.md` on this branch.

**Goal (Marco's vision).** Two user stories, one image:

1. **HPC:** download a quiverlab container (Apptainer preferred), send it to a
   SLURM cluster with EXTREMELY simple instructions, edit a sample config YAML
   naming the example (or design it in the GUI, which prints/exports the config),
   `sbatch` it, a result file lands on disk, download it, render a PDF report
   with the same container locally. Very large examples become reachable.
2. **Offline laptop app:** download the image, open it with NO internet (e.g.
   traveling), zero code — the GUI appears in the browser, computing works
   locally, with **memory limitations and time estimates shown**. Ships with
   precomputed examples (**which ones = Marco's open decision**; this plan builds
   the mechanism and a placeholder manifest).

## Architecture decisions (from the 2026-07-25 research)

- **The CLI lives in the WHEEL, the container is a convenience wrapper.** New
  console script **`quiverlab-hpc`** (`[project.scripts]`) + `python -m
  quiverlab.hpc`. New extra **`[hpc] = ["pyyaml>=6"]`**. The same CLI from
  `pip install quiverlab[fast,hpc]` is the no-container fallback — and the ONLY
  path inside the drac-local emulator (verified: no apptainer in its compute
  images).
- **Spec core promoted into the wheel:** `src/quiverlab/hpc/` — `spec.py`
  (parse dict → validate with STDLIB only (no pydantic in the base wheel) →
  dispatch over the PUBLIC `import quiverlab` surface → result dict),
  `report.py` (result → LaTeX/HTML/text via the existing `quiverlab.trace`
  pipeline), `cli.py`, `__main__.py`. `webapp/server/runner.py` DELEGATES its
  dispatch to the spec core (pydantic stays at the HTTP boundary;
  **byte-identical results and Plan-25 cache keys pinned by tests**).
  `docs/gui/runner.py` delegation is a stretch goal (the wheel is present under
  Pyodide), NOT required this plan; its interface-freshness pins must stay green.
- **Config = the request schema as YAML.** Same shape as `ComputeRequest`
  (schema 1/2: `algebra` family|quiver, `compute`, `artifacts`, `module`,
  `ext_target`) plus an `hpc:` block (`checkpoint_dir`, `time_limit_s`,
  `max_mem_bytes`, `allow_large`, `prime`, engine knobs `max_cells`/`engine`
  threaded through to `hochschild_*`). JSON accepted wherever YAML is.
  A GUI-exported config drops in unchanged (parity test).
- **Long-job story = `engine/deepen.py`** (currently orphaned): the `run` verb
  routes big HH jobs through a thin public wrapper over `deepen(A, ckpt_dir,
  time_limit_s, max_transient_bytes, ...)` — atomic per-degree checkpoints,
  resume on rerun, `finalize_only` recovery. Exit codes (sysexits): 0 done;
  **75 = clean checkpoint stop, resumable** (sbatch requeues); 65 bad
  config/relation error; 64 usage; 73 cannot write; 70 internal.
- **Render = reuse the library's own LaTeX.** `quiverlab.trace`
  (`render_latex`/`render_html`/`render_text` + `writer.py` engine detection)
  extended to consume the full `result.json` (all invariant blocks, references).
  Ladder: tectonic/pdflatex → PDF; else self-contained HTML; `--format txt`
  always works. NO KaTeX/weasyprint/typst/chromium. The image pre-warms
  tectonic's bundle at build time (fully offline at runtime).
- **Result envelope:** existing runner dict + **`result_schema: 1`** (int, new).
  `render` refuses newer `result_schema` loudly (65) and warns on version skew.
  Dims/metadata by default; heavy payloads (bases) only on request to a sidecar
  file with a size guard (`--allow-large`).

## The verbs

```
quiverlab-hpc run config.yaml -o result.json [--checkpoint-dir D] [--time-limit S]
quiverlab-hpc render result.json -o report.pdf [--format pdf|html|txt|auto]
quiverlab-hpc sample-config            # annotated YAML to stdout
quiverlab-hpc estimate config.yaml     # tier + suggested --time/--mem (estimator)
quiverlab-hpc gui [--port 8000]        # offline local webapp (below)
quiverlab-hpc version | selftest
```
Progress goes to stderr (`deepen` log callback); stdout stays clean. The CLI
sets `NUMBA_NUM_THREADS`/`OMP_NUM_THREADS` from `$SLURM_CPUS_PER_TASK` (fallback
2) unless already set.

## Offline GUI mode (Marco, 2026-07-25 addition)

`quiverlab-hpc gui` serves the Plan-09 webapp on localhost with an embedded
worker — fully offline by construction (server-side compute, vendored KaTeX, no
CDN): `create_app` + a worker loop in-process (or child), SQLite under a user
data dir (`~/.quiverlab/` or `$QUIVERLAB_DATA`), big-job/email tier disabled
(everything is local; no SMTP), and a startup banner "open
http://localhost:8000".

- **Memory + time visibility:** the estimator's `{cells, minutes}` estimate is
  already surfaced; add a **memory estimate** to `estimator.classify` output and
  the result pages, plus the container/host memory context (worker caps —
  `RLIMIT_AS` config — shown in the GUI footer/status so a laptop user sees the
  limits they're computing under). Exact wording bilingual (EN/ES) like existing
  pages.
- **Precomputed examples:** build-time seeding of the Plan-25 `result_cache` —
  `webapp/precomputed/manifest.yaml` (list of request specs), a seeding script
  that runs them at IMAGE BUILD and ships the seeded cache DB in the image
  (copied to the data dir on first run, version-keyed so a library bump
  invalidates per Plan-25 semantics). The GUI's "previously computed" note
  already exists. **The manifest ships with 3–5 placeholder examples clearly
  marked `# placeholder — Marco to curate`;** in the [web]-less wheel nothing
  changes.
- GUI export: a **"Config (YAML)"** button in `docs/gui/` (next to
  TikZ/JSON/Copy-Python; serialize `buildRequest()` + `# quiverlab-hpc run
  this-file -o result.json` header comment) and the same on the webapp result
  page ("Reproduce on a cluster" block printing the YAML next to the existing
  "Reproduce locally" Python snippet).

## Image build & distribution

- **Primary: OCI → GHCR.** `docker build` in CI on tag push; tag = library
  version (fail on mismatch with the git tag), `:latest` floats; base
  `python:3.12-slim` pinned by digest; labels + `QUIVERLAB_IMAGE_VERSION`; user
  installs with `apptainer pull docker://ghcr.io/<owner>/quiverlab:<ver>`
  (rootless, no SIF hosting; same image runs under docker).
- **Fallback: SIF release asset** (CI `apptainer build` from the OCI, assert
  < 2 GiB, attach to the GitHub Release) for air-gapped clusters.
- Contents: `pip install .[fast,hpc,web]` (web included FOR THE OFFLINE GUI;
  still no uvicorn needed for `run`/`render` paths — import-boundary test keeps
  `quiverlab.hpc` free of webapp imports except inside the `gui` verb), tectonic
  + pre-warmed bundle, `ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 MPLBACKEND=Agg
  PYTHONUNBUFFERED=1`, `ENTRYPOINT ["quiverlab-hpc"]`. **No `[qpa]`** (test-only
  oracle, huge). Optional `container/quiverlab.def` mirrors it for native SIF
  builds (`%runscript exec quiverlab-hpc "$@"`).
- Assets live in `container/` (Dockerfile, def, seed manifest hook, README) —
  the webapp deploy assets in `webapp/deploy/` stay untouched.

## SLURM story

`slurm/quiverlab.sbatch` (generic) + `slurm/quiverlab-drac.sbatch` (account +
`module load apptainer` header) + array variant in comments:

- Auto-detect launcher: `QUIVERLAB_LAUNCH` env override, else `apptainer` →
  `singularity` → `python -m quiverlab.hpc` (the emulator/venv path).
- `--bind $PWD,$SCRATCH`; checkpoint dir on `$SCRATCH`; `#SBATCH
  --signal=B:USR1@180` + handler → clean checkpoint stop (exit 75) → `scontrol
  requeue`; plain re-`sbatch` also resumes (deepen reads `latest.txt`).
- `APPTAINER_CACHEDIR=$SCRATCH/.apptainer` in the instructions (quota safety).
- **Docs page `docs/hpc.md`: the "5 steps" instructions** (pull/scp; sample-config
  or GUI export; sbatch; scp result back; render locally — with the three local
  render options: apptainer / docker / pip CLI). Genuinely minimal, assumes only
  ssh+scp+sbatch.
- **Local smoke test (acceptance):** drac-local (`slurm-local prepare` with
  `-e .[fast,hpc]` in `requirements-cluster.txt`, `activate narval`, sbatch the
  template with a tiny config via the venv fallback, assert `result.json`
  appears and renders). Documented in the plan + a `container/SMOKE.md` runbook;
  not a pytest (needs Docker).

## Tests (host pytest, no Docker needed)

- `tests/hpc/` (fast bucket by default): CLI≡API parity on a tiny config (run
  result equals direct public-API computation); `sample-config` output parses +
  validates; config-parity fixture = a GUI-shaped request dict validates
  identically in spec core and webapp pydantic (field names/semantics);
  renderer goldens (tex/html/txt token asserts; PDF only if tectonic on PATH —
  skip cleanly); `render` refuses future `result_schema` (exit 65); import
  boundary (`import quiverlab.hpc` pulls no fastapi/webapp/sqlite-server, and
  base `import quiverlab` unchanged — extend
  `tests/webapp/test_dependency_isolation.py` pattern); exit codes.
- Deep: checkpoint-resume end-to-end (run with tiny `time_limit_s` → 75 +
  checkpoint → rerun → completion, final result identical to uninterrupted).
- `tests/webapp/`: runner-delegation byte-stability (result dicts AND
  `canonical_key` unchanged on the existing fixture corpus); estimator memory
  field; gui-verb app factory smoke (TestClient, offline flags, seeded-cache
  lookup path); seeding script unit test (manifest → cache rows, math-only).
- `tests/release/test_workflows.py` conventions for the new
  `.github/workflows/container.yml` (build → smoke run+render inside docker →
  push GHCR → optional SIF release asset).
- `tests/gui/`: export-button request serialization (js-side logic kept thin;
  test the runner-side helper if any); freshness pins stay green.

## Acceptance

Wheel: `[hpc]` extra, console script, `quiverlab/hpc/` spec core + report +
CLI; runner delegation with byte-stability pins; deepen wired with exit-75
contract; estimator memory estimate; offline `gui` verb + seed mechanism +
placeholder manifest; GUI/webapp config-export affordances; `container/` +
`slurm/` assets; `container.yml` workflow; docs (`docs/hpc.md` + offline-app
section, mkdocs nav, README section); `docs/verification.md` additions (CLI≡API
parity oracle, renderer goldens, honest scope: real-cluster Apptainer manual,
emulator covers the venv-fallback orchestration only, `--mem` OOM validated by
host deepen tests); backlog entry added + ticked; CLAUDE.md/ROADMAP updates.
All suites green (fast + relevant deep + webapp) in this worktree's venv.
Merge/push only when Marco asks. **Open decision left for Marco: the curated
precomputed-example list** (placeholder manifest ships).

## Addendum (Marco, 2026-07-25): host resources & accelerators

- **`quiverlab/hpc/resources.py::detect_resources()`** — stdlib-only, exact-int
  detection of cores (cpu_count / sched_getaffinity / cgroup v1+v2 quotas /
  `$SLURM_CPUS_PER_TASK`), RAM (`/proc/meminfo` | macOS sysctl | cgroup limits |
  SLURM mem vars), and GPUs (`nvidia-smi -L`, report-only). Thread caps and the
  deepen memory budget default from it (env/SLURM always win); `estimate`
  compares the job against the host honestly; the offline GUI shows detected
  resources + configured caps.
- **GPUs: detected but UNUSED, by design.** quiverlab's engines are exact CPU
  arithmetic (int64 mod p, exact rationals); GPUs would idle. Exact mod-p linear
  algebra on GPUs (FFLAS-style delayed-reduction or integer kernels) is a real
  research direction but a large engineering+certification lift for uncertain
  gains on memory-bound ranks — recorded as a Tier-2 performance exploration,
  NOT in this plan. All resource guidance (sbatch templates, docs, DRAC-cloud
  hosting note) says: request cores+RAM, never GPUs.

## Out of scope

Real-cluster execution (manual release checklist); QPA in the image; Pyodide
canvas inside the container (the webapp surface serves the offline GUI; porting
the canvas to the webapp is a follow-up); in-flight job dedup; Windows-native
container guidance beyond docker; result-bases streaming formats.
