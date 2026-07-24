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
   auto-reads it from this directory. Optional secrets (SMTP, `QLWEB_ADMIN_TOKEN`)
   are set the same way — see "Secrets and feedback" below.
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
  `/admin/feedback` route does not exist at all.** With it set, read submissions at:

      https://quiverlab.<domain>/admin/feedback?token=<QLWEB_ADMIN_TOKEN>

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
- `QLWEB_BIG_JOB_WALL_SECONDS` (default 14400 = 4 h), `QLWEB_BIG_JOB_MEM_BYTES`
  (default 16 GB), `QLWEB_PER_EMAIL_WEEKLY_MAX` (5), `QLWEB_BIG_QUEUE_MAX` (20).
- `QLWEB_DOCS_URL` — optional; when set, a "Docs" link appears in the nav.

Emails are used only for verification + a completion notice, hashed for
rate-limiting, deleted right after the completion email, and never shown in the
admin feedback view.

**v2+ burst tier (recorded, not built):** oversized jobs can be dispatched to
short-lived RAS *compute* instances (80 vCPU / 300 GB, 1-month wall-time) via
the OpenStack API; out of scope for v1.
