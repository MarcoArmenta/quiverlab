# quiverlab-web — Design (Plan 09)

**Status:** Approved by Marco 2026-07-18 (brainstorm dialogue, sections A/B/C).
**Implements after:** Plans 03–08 (library complete). The implementation plan
(`docs/plans/…-plan-09-web.md`) is written when Plan 09 starts, per roadmap convention.
**Location:** `webapp/` in this repo. Web dependencies live in the `[web]` optional
extra so `pip install quiverlab` stays dependency-light.

## 1. Goal and audience

A public web GUI so algebraists can use quiverlab without touching code: build a
finite-dimensional algebra, pick invariants, compute on the server, read the results
with rendered mathematics, and download artifacts (JSON, worked-steps PDF, TikZ).
Public means no accounts: anyone can compute within limits; the page always shows the
`pip install quiverlab` escape hatch for going further locally.

V1 interaction is **catalog-driven**: choose a family from the library's `zoo()`
catalog, set its parameters, choose the field (CC, GF(p), GF(p^n)), tick invariants.
The **quiver canvas editor is v2**, but the v1 request schema is already
editor-shaped (§5): the editor will emit the same JSON the catalog form emits.

## 2. Architecture

Three components in `webapp/`, all pure Python plus static assets; no node toolchain.

| Component | Role |
|---|---|
| `webapp/server` | FastAPI app: server-rendered pages (Jinja2, KaTeX, vanilla JS) + JSON API under `/api/` |
| `webapp/worker` | Worker process pool consuming the jobs table; runs quiverlab in-process with thread caps |
| `webapp/deploy` | `docker-compose.yml` (app + worker + Caddy with automatic TLS), provisioning doc |

State is **one SQLite database** (WAL mode) plus an **artifacts directory** —
both on a mounted volume (`/data`). No Redis, no message broker: the jobs table is
the queue (SQLite WAL supports the single-writer/many-reader pattern this needs).

## 3. Pages

- `/` — build form: family dropdown (from `zoo()`), parameter fields, field picker,
  invariant checklist (HH^n/HH_n ranges, Coxeter polynomial, Cartan matrix, gl.dim,
  center, …: exactly the library's public invariant surface), "worked-steps PDF"
  toggle, Compute button. Instant-tier results render on this page directly.
- `/job/{id}` — permalink: queued/running (live progress, §7) → results (dims tables,
  polynomials in LaTeX via KaTeX, download buttons) or failure (§9).
- `/about` — what this is, citation (JOSS), library version, local-install pointer.
- `/api/…` — the same functionality as documented JSON endpoints:
  `POST /api/compute` (sync tier; 422 with `{queued_hint: true}` if too big),
  `POST /api/jobs`, `GET /api/jobs/{id}`, `GET /api/catalog` (families + parameter
  schemas + field constraints, generated from the library, never hand-maintained).

## 4. Two-tier compute

Every request is **cost-estimated server-side** before running, using the library's
own data (algebra dimension from the finiteness certificate, requested degree, field):

- **Instant tier:** estimated ≤ ~2 s CPU → run synchronously in the request
  (Cartan, Coxeter polynomial, bar HH on small algebras to modest degree).
- **Queued tier:** everything else → job row + `/job/{id}` permalink.

The estimator is a heuristic; a hard safety net backs it: instant-tier execution is
killed and converted to a queued job if it exceeds ~5 s wall.

## 5. Request schema (editor-ready)

One JSON schema for both tiers, versioned (`"schema": 1`):

```json
{
  "schema": 1,
  "algebra": {
    "kind": "family",           // v2 adds: "quiver"
    "family": "truncated_polynomial",
    "params": {"n": 3},
    // v2 "quiver" kind: {"vertices": [...], "arrows": [...], "relations": ["a*b - c*d", ...]}
    "field": {"kind": "GF", "p": 5, "n": 1}   // or {"kind": "CC"}
  },
  "compute": ["hh_cohomology:0..6", "coxeter_polynomial", "cartan"],
  "artifacts": {"pdf": true, "tikz": false}
}
```

Relation strings (v2) go only through the library's parser. There is no code
execution path anywhere in the app (§11).

## 6. Job model and limits

Jobs table columns: `id` (ULID), `spec` (JSON), `status`
(pending/running/done/failed), `progress` (JSON, §7), `created/started/finished`,
`quiverlab_version`, `artifact_dir`, `error`.

- **Worker pool:** default `vCPUs − 2` single-job worker processes; each sets the
  engine thread caps (`NUMBA_NUM_THREADS`/`OMP_NUM_THREADS`) so a job uses ~1 core.
- **Per-job caps** (configurable): wall-time 15 min, memory 4 GB via
  `resource.setrlimit`, result-size cap.
- **Public protection** (no accounts): 1 running job per IP + small daily per-IP
  budget; global queue cap that refuses politely with the local-install hint.

## 7. Progress

The engine's checkpointed deepen already yields per-degree checkpoints; the worker
writes them to `progress` and `/job/{id}` polls — real progress ("degree 7 of 12
certified"), not a spinner.

## 8. Artifacts, retention, reproducibility

`artifacts/{job_id}/`: always `result.json`; `trace.pdf` when the PDF toggle is on
(Plan 07 worked-steps renderer); `tikz.tex` on request. Retention 90 days
(configurable; the 10 TB volume quota makes this comfortable), stated on the page.
Every result page shows the exact quiverlab version and a **copy-paste Python
snippet reproducing the computation locally** — the GUI-to-library bridge.

## 9. Failure behavior

The library fails loudly by design; the app never softens that. `ExactnessError`,
`RelationError`, `DepthLimitError`, … surface verbatim (message + type) in a
clearly-framed failure box; queued failures mark the job `failed` with the same.
No stack traces to the public; full tracebacks go to server logs.

## 10. Deployment (DRAC cloud, RAS)

Target: an Arbutus **persistent instance** (RAS quota: 25 vCPU / 50 GB / 10
instances) — provision ~16 vCPU / 50 GB for app + workers, ports 80/443/22, one
floating IP, an attached volume at `/data`. Marco points a domain/subdomain at the
floating IP; Caddy handles TLS automatically. `webapp/deploy/PROVISIONING.md` walks
the OpenStack steps (flavor, security groups, floating IP, volume, compose up).
Optional config-flagged backups: nightly SQLite `.backup` + artifact sync to object
storage. The design is host-generic: any Linux box with Docker runs the same compose
file. **v2+ burst tier** (recorded, not built): dispatch oversized jobs to
short-lived RAS compute instances (80 vCPU / 300 GB, 1-month wall-time) via the
OpenStack API.

## 11. Security

- **No user code execution.** Inputs: family id + typed parameters + field spec
  (+ v2 relation strings through the library parser only). No eval, no pickle.
- Request-size caps; JSON-only API; strict CSP; no cookies, no accounts.
- Rate limiting at both Caddy and app layers; IP-hashed logs.
- ULID job ids: unlisted-but-public permalinks (results contain no personal data).
- SQLite and artifacts are never directly exposed; downloads stream through the app
  with `Content-Disposition`.

## 12. Testing

- API contract, tier decision, and limit enforcement via FastAPI `TestClient`.
- End-to-end worker test: real worker process, temp SQLite, tiny jobs through
  pending→running→done and →failed paths.
- Golden tests for the `result.json` schema and `/api/catalog` output.
- Page assertions (results render, KaTeX present, downloads linked).
- `docker compose config` validated in CI; full-stack smoke at deploy time.

## 13. Non-goals (v1)

Accounts/auth, the canvas editor (v2, schema-ready), burst compute (v2+), GPU
anything, result sharing beyond permalinks, non-catalog algebras (v2), HPC/SLURM
integration.
