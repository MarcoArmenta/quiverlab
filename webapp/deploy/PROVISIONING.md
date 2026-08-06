# Provisioning quiverlab-web on DRAC Arbutus (RAS)

Target: a **persistent instance** (RAS quota: 25 vCPU / 50 GB RAM / 10 instances).

1. **Instance.** In the Arbutus Horizon dashboard, launch an instance with a
   flavor around `p8-16gb` … `c16-60gb` (≈16 vCPU / 50 GB), the **current Ubuntu
   LTS image (Ubuntu 24.04 LTS)**, and your SSH key. Boot from a new volume (40 GB+).
2. **Volume.** Create and attach a data volume (e.g. 200 GB) for `/data`
   (artifacts + SQLite); format and mount it at `/data`.
3. **Security group (VM security baseline — DRAC "Security considerations").**
   Create a dedicated security group (do **not** edit the `default` group) with
   exactly three ingress rules and nothing else:
   - `80/tcp` from `0.0.0.0/0` — the sanctioned public web port.
   - `443/tcp` from `0.0.0.0/0` — the sanctioned public web port.
   - `22/tcp` from **`<YOUR.IP.CIDR/32>`** — replace with Marco's own IP/CIDR.
     **Never open SSH to `0.0.0.0/0`, and never put SSH in the `default` group.**
   Add **no** other ports, no port ranges, no RDP/VNC/database ports.
4. **Floating IP.** Allocate one floating IP (RAS quota: 2) and associate it.
   Point a DNS A record (e.g. `quiverlab.<domain>`) at it; put that name in the
   Caddyfile.
5. **Host hardening.**
   - SSH: key auth only (the cloud image default); never enable password auth.
   - `sudo apt-get update && sudo apt-get install -y fail2ban unattended-upgrades`,
     then `sudo systemctl enable --now fail2ban` and enable weekly security
     updates: `sudo dpkg-reconfigure -plow unattended-upgrades` (or set
     `APT::Periodic::Unattended-Upgrade "7";`). Document a manual cadence too:
     roughly monthly `sudo apt-get update && sudo apt-get dist-upgrade && sudo reboot`.
6. **Docker.** `sudo apt-get install -y docker.io docker-compose-plugin`.
7. **Secrets (REQUIRED — do this before deploying).** The stack refuses to start
   without two secrets, wired into `docker-compose.yml` with the required-var
   syntax `${VAR:?...}`. Clone the repo, then create `webapp/deploy/.env` from
   the template and fill both in:
   ```bash
   git clone https://github.com/MarcoArmenta/quiverlab.git
   cd quiverlab/webapp/deploy
   cp .env.example .env
   printf 'QLWEB_IP_HASH_SALT=%s\nQLWEB_TOKEN_SECRET=%s\n' \
       "$(openssl rand -hex 32)" "$(openssl rand -hex 32)" > .env
   ```
   - `QLWEB_IP_HASH_SALT` — client IPs are stored only as `sha256(salt+":"+ip)`;
     rotating the salt anonymizes past rate-limit keys.
   - `QLWEB_TOKEN_SECRET` — signs the single-use big-job magic-link tokens (§17).
   Keep `.env` off version control (it holds live secrets). `docker compose`
   auto-reads it from this directory. Optional secrets (SMTP, `QLWEB_ADMIN_TOKEN`,
   `QLWEB_EMAIL_HASH_SALT` — salt for the per-email rate-limit hash, defaults to
   `QLWEB_TOKEN_SECRET` when unset) are set the same way — see "Secrets and
   feedback" below.
8. **Deploy.**
   ```bash
   # from quiverlab/webapp/deploy (repo already cloned in step 7)
   # edit Caddyfile: set your domain
   sudo docker compose up -d --build
   ```
   `compose up` fails fast if either required secret from step 7 is missing.
   Caddy obtains TLS automatically. Verify: `https://quiverlab.<domain>/`.
   KaTeX is vendored in the repo (`webapp/static/katex/`) and copied into the
   image, so the build needs no CDN and the app serves math with a strict CSP.
9. **Web listeners.** HTTPS only. Caddy binds 80 solely to redirect to 443
   (the `Caddyfile` site block on `:443`; the automatic HTTP→HTTPS redirect is
   Caddy's default) and 443 for the app. Nothing else listens on the host — no
   mail server (all email is the outbound SMTP relay, §17), no database port
   (SQLite is a file), no BitTorrent, no extra daemons. Confirm with
   `sudo ss -tlnp` — only sshd (bound to your CIDR via the security group),
   plus Docker-published 80/443.
10. **Backups (optional).** Cron a nightly `sqlite3 /data/quiverlab_web.sqlite3
   ".backup /data/backup.sqlite3"` and sync `/data/artifacts` to object storage
   (RAS object quota: 10 TB).

Retention: finished jobs and artifacts are swept after `QLWEB_RETENTION_DAYS`
(default 90). Tune limits via `QLWEB_*` env in `docker-compose.yml`.

## Cloud capacity tuning

The deployed website runs on real hardware and is meant to do HPC-grade work the
downloadable laptop app cannot. `docker-compose.yml` therefore ships a **cloud
profile** tuned for the RAS persistent instance (~**16 vCPU / 50 GB RAM**, RAS
quota 25 vCPU / 50 GB, plus a **200 GB `/data` volume**). Every value is a
`${VAR:-tuned_default}` in the compose `environment:` — applied out of the box,
overridable from `.env` (see `.env.example`) without editing compose. The
**instant tier is deliberately left at the library code defaults** — it is the
anonymous, unauthenticated DoS surface, so widening it is off the table.

### The RAM arithmetic (the one hard rule)

Each job runs in a `spawn` child that pins every math runtime to **one thread**,
and the queue is a **single shared FIFO** (`JobStore.claim_next`, oldest-first,
no per-tier priority). So the worst case is **every worker running a big job at
once**, and the binding constraint is

```
QLWEB_WORKER_PROCESSES  ×  QLWEB_BIG_JOB_MEM_BYTES   ≤   RAM − headroom
```

where `QLWEB_BIG_JOB_MEM_BYTES` is the per-job `RLIMIT_AS` the worker enforces
(hard on Linux). Budgeting **~10 GB of the 50 GB** for the OS, the FastAPI app
process, Caddy, the SQLite/WAL page cache, and the app-tier instant children
leaves **~40 GB** for the worker fleet:

```
2 workers × 20 GiB (big)  = 40 GiB   ≤ 40 GiB usable      ✔ (worst case, all big)
2 workers ×  8 GiB (queued)= 16 GiB                        ✔
1 big (20 GiB) + 1 queued (8 GiB) = 28 GiB                 ✔
```

`RLIMIT_AS` caps **virtual** address space, and numba/LLVM reserve large virtual
regions up front, so real resident memory sits well under these ceilings — the
arithmetic above is deliberately conservative (it treats the ceiling as if fully
resident). That conservatism is the margin.

### Chosen values (knob → old → new → why)

| Knob | Code default | Cloud profile | Rationale |
|---|---|---|---|
| `QLWEB_WORKER_PROCESSES` | `cpu−2` (≈14 here) | **2** | RAM-bound, not core-bound. 2 × 20 GiB big = 40 GiB fits; keeping a second loop free means a big job never fully starves the queued tier. |
| `QLWEB_JOB_WALL_SECONDS` (queued) | 900 (15 min) | **3600 (1 h)** | The cloud accepts computations past a laptop's patience. |
| `QLWEB_JOB_MEM_BYTES` (queued) | 4 GiB | **8 GiB** | Deeper resolutions the raised thresholds admit need room. Also the **instant** child ceiling (see below); harmless, instant jobs are tiny. |
| `QLWEB_BIG_JOB_WALL_SECONDS` | 14400 (4 h) | **86400 (24 h)** | Email-verified HPC jobs get a full day (48 h available via `.env` if wanted). |
| `QLWEB_BIG_JOB_MEM_BYTES` | 16 GiB | **20 GiB** | Largest the RAM rule allows at fleet size 2 (see below on 40 GiB). |
| `QLWEB_QUEUED_OPS_THRESHOLD` | 5e8 | **2e10** | Order-of-magnitude aligned to a 1 h wall via the estimator's 500M ops/min hint (~30e9 ops in 1 h). |
| `QLWEB_QUEUED_MAX_DEGREE` | 20 | **30** | Deeper resolutions. |
| `QLWEB_BIG_OPS_THRESHOLD` | 5e10 | **5e11** | Aligned to a 24 h wall (~700e9 ops). |
| `QLWEB_BIG_MAX_DEGREE` | 40 | **60** | Deeper campaigns. |
| `QLWEB_CACHE_MAX_ENTRIES` | 1000 | **10000** | 200 GB volume; the cache is an entry-count LRU — see the byte caveat below. |
| `QLWEB_RETENTION_DAYS` | 90 | **365** | A year of artifacts for reproducibility on the 200 GB volume. |
| `QLWEB_GLOBAL_QUEUE_MAX` | 200 | **1000** | Deeper backlog; concurrency is still the fleet size (jobs drain 2 at a time). |
| `QLWEB_BIG_QUEUE_MAX` | 20 | **50** | Deeper big backlog. |
| `QLWEB_PER_EMAIL_WEEKLY_MAX` | 5 | **10** | More generous on the bigger box. |
| `QLWEB_PER_IP_DAILY_MAX` | 100 | **300** | More generous daily budget. |
| `stop_grace_period` (compose) | 930s | **3630s** | = queued wall + 30s slack (see below). |
| INSTANT knobs, `*_RUNNING_MAX` | — | **unchanged** | Instant is the DoS surface; per-IP / per-email *running* caps stay 1 (the concurrency + RAM bound). |

The per-IP and per-email **running** caps stay at **1** on purpose: they bound
how many jobs one identity can have in flight, which (with the fleet size) is
what keeps the RAM arithmetic honest. The *backlog* and *weekly/daily* budgets
are what we widened.

### Why not a single 40 GiB big job?

A 40 GiB big cap was evaluated. On the shared FIFO queue it forces
`QLWEB_WORKER_PROCESSES = 1` (1 × 40 GiB = 40 GiB), which means a 24 h big job
would **freeze the anonymous queued tier for up to a day** — the second loop that
normally keeps serving queued work while a big job runs would not exist. That
breaks the tier philosophy (email gates the *cost* of computing, never a queued
user's access). So the profile keeps 20 GiB at fleet size 2. If you knowingly
want to dedicate the whole box to one huge run, the `.env` "single-huge-job
campaign profile" (`QLWEB_WORKER_PROCESSES=1`, `QLWEB_BIG_JOB_MEM_BYTES=40 GiB`)
is documented — but the right home for campaign-scale work is the burst instance.

### Burst instance for campaign-scale work

RAS also offers short-lived **compute** instances (**80 vCPU / 300 GB RAM**,
1-month wall) for campaigns. These are NOT where the web service lives; they are
for batch computation, driven by the in-wheel **`quiverlab-hpc`** CLI (Plan 28):
`quiverlab-hpc run` executes the same spec the webapp runner does (byte-stable
results, Plan-25 cache keys pinned by frozen goldens), with checkpoint/resume
(exit 75 = clean checkpoint stop) and SLURM templates. Spin one up, run the
campaign with `quiverlab-hpc`, harvest the artifacts to object storage (RAS
object quota 10 TB), tear it down. A finished result can then be seeded into the
persistent site's Plan-25 result cache so it replays instantly for everyone. The
persistent web instance handles interactive + moderate jobs; the burst instance
handles the 40+ GiB / many-parallel campaigns.

### What still will not fit (honest scope)

- **The Plan-35 `dim ≥ 220` product examples** (`nakayama-kz20-deep`,
  `nakayama-kz24-deep`): the bar/TT `to_engine` step alone is ~290 s and the
  degree-2 cochain basis is ~10.5M cells (over `max_cells`), so a `cup`/`cap`
  probe does not finish even in a 24 h / 20 GiB big slot (a direct 25-min kZ20
  `cup:0..2` timed out — evidence in `webapp/precomputed/manifest.yaml`). These
  carry no products by design; they are a burst-instance / offline-`quiverlab-hpc`
  job, not a website job.
- **Deep non-monomial HH past ~degree 10 at `dim ≳ 30`**, `Π(D₅)` HH-at-scale,
  `Λ(kⁿ≥4)` depth, and `decompose` past ~`dim 50` are the standing deferred
  cluster (see the verification page honest-scope section) — the box makes them
  *more* reachable but does not change their asymptotics.
- The **result cache is an entry-count LRU, not a byte budget**: 10000 entries is
  safe only while cached results stay small (KB–MB, which curated examples are).
  Worst case per entry ≈ `QLWEB_RESULT_MAX_BYTES` (32 MB) plus the trace HTML, so
  watch `df -h /data` and lower `QLWEB_CACHE_MAX_ENTRIES` / `QLWEB_RETENTION_DAYS`
  if the volume fills.

### Graceful stop with the longer walls

`stop_grace_period` is **3630s** = the queued wall (3600) + 30s slack, kept just
above `run_loop.main()`'s own join deadline (queued wall + 20s) so `main` reaps
its loops cleanly before Docker's SIGKILL. A **big job (24 h) is not covered** —
a deploy must not block for a day. Any job still running at the cut (a big job, or
a queued job that outlives the grace) is SIGKILLed and **requeued** from
`running` → `pending` at the next startup (`JobStore.requeue_stale_running`), so
nothing is lost and no running slot leaks. For an immediate restart use
`docker compose kill` and let the startup requeue recover the in-flight job. If
you raise `QLWEB_JOB_WALL_SECONDS` past 1 h, raise `stop_grace_period` to match.

## Result cache (Plan 25)

Finished results are cached and replayed for any later identical request, across
users and tiers — a previously computed example is never recomputed, and a *cached*
big example is served with no email/token at all (email verification gates the cost
of NEW computations, not access to the mathematics). The cache lives in a
`result_cache` table inside the same SQLite database, keyed by the canonicalized
request plus the library version (a version bump invalidates every entry naturally),
and carries mathematics only — no email, IP, or token. A cache entry "pins" its
finished job against the retention sweep so the artifacts survive to back replays;
when the entry is evicted (LRU size cap or a version bump) the pin lifts and ordinary
retention reclaims the job once it is also past `QLWEB_RETENTION_DAYS`.

- `QLWEB_CACHE_ENABLED` (default `1`) — set to `0` to disable the cache entirely
  (every request recomputes; the `result_cache` table is simply left empty).
- `QLWEB_CACHE_MAX_ENTRIES` (default `1000`) — LRU size cap. The least-recently-hit
  entries beyond this many are evicted by the hourly sweep (which runs alongside the
  retention sweep on worker loop 0). Higher = more replays kept warm (and more
  artifacts pinned on disk); lower = tighter disk footprint.

## Worked-steps PDF (optional TeX engine)

The default image ships **without** a TeX engine. This is deliberate: the compute
runner degrades honestly. When no `pdflatex`/`tectonic` is on `PATH`, the
worked-steps artifact is written as a self-contained, no-JS **HTML** transcript
instead of `trace.pdf` (see `webapp/server/runner.py::_PDF_HTML_FALLBACK`); the
mathematics is identical, only the container format differs. Users still get
correct, complete worked steps.

Rationale for the default: a TeX engine adds real image cost (a full TeX Live is
hundreds of MB; even Tectonic is ~50 MB plus a build-time download that can fail
offline) for an artifact most requests do not ask for (`artifacts.pdf` defaults
to `false`). Keeping it out yields a lean, network-light, reproducible build.

To render real `trace.pdf`, enable **Tectonic**: uncomment the `TECTONIC_*`
`ARG`/`RUN` block in `webapp/deploy/Dockerfile` (it downloads a pinned static
release into `/usr/local/bin` and self-tests with `tectonic --version`), then
`docker compose up -d --build`. No other change is needed — the runner picks up
`tectonic` automatically once it is on `PATH`. Pin `TECTONIC_VERSION`/
`TECTONIC_ARCH` to the release matching your build architecture.

## Build context and image notes

`docker compose build` uses the **repo root** as the build context (the compose
`build.context` is `../..`). To keep the context small and avoid shipping the
local virtualenv, git history, or agent scratch into the daemon, a root
`.dockerignore` excludes `.venv/`, `.git/`, `.claude/`, `docs/`, caches, and test
trees. Edit it if you add large top-level directories that the image does not
need. The image installs `.[web,fast]`, so numba is present in production.

## Operational limits

- **No per-service CPU/memory limits by default.** `docker-compose.yml` sets no
  `deploy.resources` / `mem_limit` caps; the host's RAS flavor bounds the whole
  stack. Per-**job** compute is bounded instead by the worker: each job runs in a
  resource-capped `spawn` child with `RLIMIT_CPU` (Linux + macOS) and `RLIMIT_AS`
  (Linux, hard in the Docker image) plus a parent wall-time kill, so a runaway
  job cannot exhaust the VM. Add `deploy.resources.limits` if you want a hard
  ceiling per container as well.
- **Graceful worker stop.** `docker compose stop` SIGTERMs the worker; the poll
  loops finish their in-flight job then exit, bounded by the worker's
  `stop_grace_period` (930s ≈ one job wall + slack). Any job still running when
  Docker SIGKILLs — or lost to a reboot/crash — is requeued from `running` back
  to `pending` at the next worker startup (`JobStore.requeue_stale_running`), so
  no job is stranded and no per-IP running slot leaks.
- **`/data` volume ownership.** The image runs as non-root uid `10001`. A FRESH
  `qldata` named volume inherits writable ownership from the image. A volume that
  already exists from an older root-based image is root-owned; fix it once on the
  host with `sudo chown -R 10001:10001 /var/lib/docker/volumes/deploy_qldata/_data`
  (adjust the volume name to `docker volume ls`).
- **Caddyfile is mounted read-only** (`:ro` in compose); the app container never
  writes it. Caddy's ACME state lives in the separate `caddy_data` volume.

## Secrets and feedback

`QLWEB_IP_HASH_SALT` and `QLWEB_TOKEN_SECRET` are **required** (step 7): compose
wires them with `${VAR:?...}` and refuses to start if either is missing, so no
dev-default secret can ship. Put them in `webapp/deploy/.env` (see `.env.example`).
The **optional** secrets below go in the same `.env`, or in the `app`/`worker`
`environment:` blocks:

- `QLWEB_IP_HASH_SALT` — a long random string. Client IPs are stored only as
  `sha256(salt + ":" + ip)`; rotating the salt anonymizes past rate-limit keys.
- `QLWEB_ADMIN_TOKEN` — enables the feedback admin view. **When unset, the
  `/admin/feedback` route does not exist at all.** With it set, pass the token in
  the `X-Admin-Token` **header** (never a query string — uvicorn's access log
  records query strings, so a `?token=...` would leak the admin secret into the
  logs):

      curl -H "X-Admin-Token: $QLWEB_ADMIN_TOKEN" https://quiverlab.<domain>/admin/feedback

  The comparison is constant-time. Alternatively, query the SQLite directly on
  the VM (no token needed, you are on the box):

      sqlite3 /data/quiverlab_web.sqlite3 \
        "SELECT created_at, category, message, contact, job_ref FROM feedback ORDER BY created_at DESC LIMIT 50;"

- `QLWEB_FEEDBACK_DAILY_MAX` (default 5) caps submissions per hashed IP per day.

Feedback is public and no-account; the page also links to the GitHub Issues
tracker for users who prefer that path.

### Big jobs and email (spec §17)

The big-job tier needs an **outbound SMTP relay** (the VM runs no mail server).
Set, in the `app` and `worker` `environment:` blocks:

- `QLWEB_SMTP_HOST`, `QLWEB_SMTP_PORT` (default 587), `QLWEB_SMTP_USER`,
  `QLWEB_SMTP_PASS`, `QLWEB_SMTP_FROM` — the relay credentials and From address.
  **If `QLWEB_SMTP_HOST`/`QLWEB_SMTP_FROM` are unset, the big-job tier is
  disabled** and the app tells users to run locally.
- `QLWEB_TOKEN_SECRET` — **required** (already set in step 7): the random string
  signing the single-use magic-link tokens. Rotating it invalidates outstanding
  links.
- `QLWEB_PUBLIC_BASE_URL` — e.g. `https://quiverlab.<domain>`; used to build the
  verification and completion links in emails.
- `QLWEB_BIG_JOB_WALL_SECONDS`, `QLWEB_BIG_JOB_MEM_BYTES`,
  `QLWEB_PER_EMAIL_WEEKLY_MAX`, `QLWEB_BIG_QUEUE_MAX` — the library **code**
  defaults are 4 h / 16 GB / 5 / 20, but the shipped cloud profile in
  `docker-compose.yml` raises them to **24 h / 20 GiB / 10 / 50** (see
  "Cloud capacity tuning" above for the arithmetic and how to override).
- `QLWEB_DOCS_URL` — optional; when set, a "Docs" link appears in the nav.

Emails are used only for verification + a completion notice, hashed for
rate-limiting, deleted right after the completion email, and never shown in the
admin feedback view.

**Campaign-scale / burst tier.** Jobs that exceed even the 24 h / 20 GiB big
slot belong on a short-lived RAS *compute* instance (**80 vCPU / 300 GB**,
1-month wall-time), driven **today** by the in-wheel `quiverlab-hpc` batch CLI
(Plan 28) — see "Cloud capacity tuning → Burst instance" above. Fully *automatic*
dispatch of oversized web jobs to such an instance via the OpenStack API is
recorded as a v2+ item (not built); for v1 the burst path is the manual
`quiverlab-hpc` run, and its finished results can be seeded into the site's
Plan-25 result cache to replay instantly.
