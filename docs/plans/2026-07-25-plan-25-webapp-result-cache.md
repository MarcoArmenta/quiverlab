# Plan 25 — the webapp result cache

**Status:** DELIVERED 2026-07-25, branch `plan-25-webapp-result-cache` (from
`main`; the webapp tier is merged to `main`). Backlog: DEEPER-ENGINES-BACKLOG.md
Tier 3, "Result cache — never recompute a known example" (Marco, 2026-07-25).

**Goal (Marco's directive):** cache the examples computed once the webpage is set
up, so if another user asks for the same example we don't recompute. If it is a big
example, a user verifies email, we do the computation and cache it; if another user
asks for the same example, since it is cached we don't need to verify their email.

Exact computations are **deterministic**: the same canonical request always yields
the same result. A finished result is therefore a correctness-safe replay for any
later identical request — across users, across tiers — with zero recompute. Email
verification gates the **cost** of computing, not access to the mathematics.

## Canonical key

`webapp/server/cache.py::canonical_key(spec, library_version)` =
`sha256` of a canonical JSON encoding of `{"lib": <library version>, "req": <the
versioned compute-request spec>}`. The spec is the pydantic
`ComputeRequest.model_dump(by_alias=True)` — schema version + request kind
(family/quiver) + family/quiver spec + field + `compute` invariants/params +
`artifacts` flags. Canonicalisation rules (all pinned by unit tests,
`tests/webapp/test_cache.py`):

- **Sorted keys** — dict key order (`params`, `arrows`, `field`, the top-level
  dict) is irrelevant; the same mathematics always hashes identically.
- **Through JSON** — `model_dump()` yields arrow *tuples*, a store round-trip
  yields *lists*; `json.dumps` erases the distinction, so the API (in-memory tuple
  form) and the worker (reloaded `job.spec` list form) compute ONE key. This is the
  load-bearing invariant that lets the check side and the record side agree.
- **Library version is in the key** — a version bump changes every key, so old
  entries can never be replayed for a new build (natural invalidation). The schema
  version rides along inside the spec, so a schema bump invalidates too.
- **Artifact flags are in the key** — a `pdf:true` request produces a `trace.pdf` a
  `pdf:false` run does not, so they are distinct entries (a hit never 404s on a
  missing artifact).
- **Compute-list order is significant** (deliberate, conservative): the order
  changes the produced artifacts (`result.json` key order, the `reproduce` snippet,
  and — with two HH kinds — which one's worked-steps PDF is rendered). Treating a
  permuted list as a distinct entry guarantees a hit replays EXACTLY what this
  request would have produced, never a subtly-different artifact. Dict key order is
  the only thing normalised away.

Nothing user-identifying enters the key: it is a hash of mathematics + version
only. Big-job requests key on the spec with `email`/`lang` **stripped** (the same
stripped spec `pending_big`/`jobs` persists), so an anonymous **queued** run and an
email-verified **big** run of the same mathematics share one cache entry — a big
request is served from a prior anonymous computation and vice versa.

## Storage — new table, not a jobs column

A dedicated `result_cache` table (`webapp/server/store.py`), inside the single
SQLite/WAL state store, WAL/`_conn` idioms unchanged:

```
result_cache(key PK, job_id, quiverlab_version, created_at, last_hit_at, hits)
```

Chosen over a `cache_key` column on `jobs` because:

- **Separation of concerns / privacy.** The cache is an *index over finished jobs*
  plus LRU metadata. A separate table keeps the `jobs` schema stable and makes the
  "mathematics only, no PII" contract self-evident and testable (the table simply
  has no `email`/`ip`/`token` column — asserted in
  `test_cache_row_is_math_only_no_pii`).
- **Independent lifetimes.** The cache entry's lifetime (LRU/size cap) is decoupled
  from a job's own retention window; a column would fuse them.
- **No jobs-schema change** — which also settles the in-flight-dedup question below
  ("within the existing schema").

`cache_put` is idempotent (`ON CONFLICT(key) DO NOTHING`, first writer wins).
`cache_get` bumps `hits`/`last_hit_at` and self-heals a dangling row (defensive —
should never occur while pinned). `cache_sweep` version-purges then LRU-evicts.

## Artifact-lifetime rule (the subtle part) — the cache *pins* its job

A `result_cache` row references a `done` job whose artifact directory backs the
replay. The rule:

> A finished job's row and artifacts live until it is **both** older than the
> retention cutoff **and** unreferenced by a live cache entry.

Mechanism: `purge_older_than` (the ordinary retention sweep's DB step) now excludes
any job referenced by `result_cache` (`... AND id NOT IN (SELECT job_id FROM
result_cache)`). So a cached job's row + artifacts survive retention, and the
existing `/job/{id}` page and `/download/{id}/...` routes keep working unchanged
(they resolve through the still-present job row — no second artifact copy, no
download-route change). When the cache entry is later evicted (`cache_sweep`:
version purge or LRU size cap), the pin lifts and a subsequent retention pass
reclaims the now-unreferenced job if it is also past the cutoff. Eviction itself
never deletes a job or artifacts directly — a young-but-cold evicted entry's job
simply rejoins normal retention. Sweep ordering (run_loop, hourly on loop 0):
`sweep_cache_once` **first** (lift stale pins), then `sweep_once` (reclaim
now-unpinned old jobs) — same-tick convergence. Validated end-to-end by
`test_cache_lru_sweep_evicts_and_retention_reclaims`.

## Serve-from-cache on all tiers

All three request entry points check the cache **first** — before building the
algebra, classifying, or rate-limiting — so a hit costs nothing (not even a
rate-limit charge):

- `POST /api/compute` (`app.py`): `cache.lookup(...)` at the top → `_cached_response`
  (`200 {"tier":"cached","job_id","cached_at"}`). Covers the instant tier too: an
  instant-classified request that was previously forced to queue by a wall-net
  overflow is now cached, so the doomed instant attempt is skipped.
- `POST /api/jobs` (`app.py`): same check → replays instead of enqueuing a duplicate.
- `POST /api/jobs/big` (`bigjobs.py`): the crux. The cache check is the **first**
  statement in `submit_big`, **before** `big_jobs_enabled`, so a cached big example
  is served even with SMTP off (a cached result needs no relay). On a hit it returns
  `200 {"status":"cached","job_id","cached_at"}` with **no token minted, no email
  sent, no `pending_big` row, no `email_hash` stored** — the key strips email/lang,
  so it collides with the anonymous/big cache regardless of who computed it first.
  Only a genuinely new big example falls through to the magic-link flow (and enters
  the cache for everyone after). Asserted via the fake mailer in
  `test_big_job_cache_hit_sends_no_email_and_mints_no_token`.

Misses behave exactly as before. On completion the **worker**
(`worker.py::run_one_job`, success path only) calls `cache.record` →
`cache_put(canonical_key(job.spec, version), job.id, version, now)`. Failed jobs
are never cached (a failure may be transient — e.g. a wall-time kill under load).
A `cache_put` failure is swallowed with a type-only warning: the result is already
saved and must never be lost. The whole feature is gated by `cfg.cache_enabled`.

Frontend (`app.js`, CSP-safe, no new inline JS): the compute form navigates to
`/job/{id}` on `tier:"cached"`; the big-email button navigates on
`status:"cached"` (so a cached big response never renders as an error).

## In-flight dedup (races) — benign, documented, skipped

Two identical requests racing to compute the same math is **benign for
correctness** (exact/deterministic), only wasteful. Deduplicating an *in-flight*
(pending/running) job would require the canonical key on the `jobs` row (a new
indexed column + populating it at three call sites + version-consistency edge cases
across a redeploy-mid-job) — i.e. **not cheap "within the existing schema"**, the
requirement's own bar for skipping. So we skip it and rely on:

1. `cache_put`'s `ON CONFLICT DO NOTHING` idempotency — concurrent identical
   completions never crash; exactly one job stays pinned, the other ages out
   normally (`test_put_is_idempotent_first_writer_wins`).
2. The self-limiting window — the first completion caches the result, so every
   later identical request hits.

## Privacy

The cache row is mathematics only (key, job ref, version, timestamps, hit
counter) — no email hash, no IP, no token. A big job's cached result served to a
second user leaks nothing about the first: `result.json` and the `/job/{id}` page
render the algebra spec, dims, references, and reproduce snippet only (never
`email`/`ip`, which live on the `Job` row but are never templated). Verified by
`test_cache_row_is_math_only_no_pii` and the privacy asserts in the big-job
acceptance test.

## Config (env, PROVISIONING.md)

- `QLWEB_CACHE_ENABLED` (default 1) — the whole feature; 0 disables lookup +
  record.
- `QLWEB_CACHE_MAX_ENTRIES` (default 1000) — LRU size cap; bounds the pinned
  artifacts.

## UX / honesty (bilingual)

The `/job/{id}` page shows, for `done` jobs, the finish timestamp and a bilingual
note (`job.finished`, `job.cached_note` in en/es): *this exact result is cached —
identical requests are served instantly without recomputing; every number is exact,
this is the answer, not an approximation.* EN/ES parity gated by
`test_en_es_key_parity`.

## Tests

- `tests/webapp/test_cache.py` (new, unit): canonicalizer (determinism, dict-order
  invariance, tuple/list collision, version-in-key, mathematics-distinctness,
  compute-order significance) + store (`cache_get/put`, LRU + version sweep,
  retention pin, dangling self-heal, no-PII columns).
- `tests/webapp/test_acceptance.py` (extended): queued replay across users, instant
  tier checks cache first, big-job cache hit sends no email / mints no token,
  big cache hit served with SMTP off, version-bump invalidation, LRU sweep +
  retention reclaim, bilingual cached note.

**Verification-page note (merge-time TODO):** `docs/verification.md` does not exist
on this branch (concurrent branch `plan-22-verification-transparency`). When the two
merge, add Plan 25's guarantees to that page: the result cache is a
correctness-safe replay (exact ⇒ deterministic ⇒ THE answer), keyed by
canonical-request + library version (version bump invalidates), carrying no
user-identifying data, with the pin/LRU artifact-lifetime rule above.

## Out of scope (unchanged by this plan)

Other Tier-3 webapp-polish items (error-envelope unification, `app.js` /es dynamic
label localization, verify-page transient-vs-terminal, tier-ordered `claim_next`,
real-SMTP/TLS/concurrency acceptance, `test_instant_compute` de-flake) are NOT
touched here.
