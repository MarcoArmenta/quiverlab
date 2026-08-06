"""SQLite in WAL mode as the single source of truth and the job queue.
claim_next() uses an atomic UPDATE ... WHERE status='pending' guarded by a
short IMMEDIATE transaction so multiple worker processes never double-claim.

The `ip` column holds an OPAQUE salted hash, never a raw client address — the
app layer (Task 9) hashes before calling create_job/count_*_for_ip. The store
is agnostic; it just stores and counts whatever string it is given."""
from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path

from ulid import ULID

# The error a cancelled job carries. `mark_failed` records it, so the job ends in
# the terminal `failed` state the /draw poller (webapp/static/gui/worker.js) already
# recognises -- no new status value the read-only client would hang on.
CANCEL_ERROR = "Cancelled: you cancelled this job."


@dataclass
class Job:
    id: str
    spec: dict
    status: str
    progress: dict
    created_at: str
    started_at: str | None
    finished_at: str | None
    quiverlab_version: str | None
    artifact_dir: str | None
    error: str | None
    ip: str
    tier: str
    email: str | None
    email_hash: str | None
    wall_seconds: int | None
    mem_bytes: int | None
    lang: str


_SCHEMA = """
CREATE TABLE IF NOT EXISTS jobs (
  id TEXT PRIMARY KEY,
  spec TEXT NOT NULL,
  status TEXT NOT NULL DEFAULT 'pending',
  progress TEXT NOT NULL DEFAULT '{}',
  created_at TEXT NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%SZ','now')),
  started_at TEXT,
  finished_at TEXT,
  quiverlab_version TEXT,
  artifact_dir TEXT,
  error TEXT,
  ip TEXT NOT NULL,
  tier TEXT NOT NULL DEFAULT 'anonymous',
  email TEXT,
  email_hash TEXT,
  wall_seconds INTEGER,
  mem_bytes INTEGER,
  lang TEXT NOT NULL DEFAULT 'en'
);
CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status);
CREATE INDEX IF NOT EXISTS idx_jobs_ip ON jobs(ip);
CREATE INDEX IF NOT EXISTS idx_jobs_email_hash ON jobs(email_hash);

-- Big-job (spec §17) pending verifications: a signed magic-link consumes one
-- of these exactly once, turning it into a `jobs` row with big-job caps.
CREATE TABLE IF NOT EXISTS pending_big (
  id TEXT PRIMARY KEY,
  spec TEXT NOT NULL,
  email TEXT NOT NULL,
  email_hash TEXT NOT NULL,
  lang TEXT NOT NULL DEFAULT 'en',
  created_at TEXT NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%SZ','now'))
);

CREATE TABLE IF NOT EXISTS feedback (
  id TEXT PRIMARY KEY,
  category TEXT NOT NULL,
  message TEXT NOT NULL,
  contact TEXT,
  ip TEXT NOT NULL,
  created_at TEXT NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%SZ','now')),
  job_ref TEXT,
  extra TEXT
);
CREATE INDEX IF NOT EXISTS idx_feedback_ip ON feedback(ip);

-- Result cache (Plan 25): a finished computation is replayable for any later
-- identical request, keyed by the CANONICALISED request + library version (see
-- webapp/server/cache.py). Mathematics ONLY -- no email/ip/token ever lives here.
-- `job_id` references a `done` job whose artifacts back the replay; a live cache
-- row "pins" that job against the ordinary retention purge (see purge_older_than),
-- so the artifacts survive to be replayed. `hits`/`last_hit_at` drive the LRU
-- size-cap sweep (cache_sweep); `quiverlab_version` lets a version bump purge stale
-- rows even before the key mismatch would have hidden them.
CREATE TABLE IF NOT EXISTS result_cache (
  key TEXT PRIMARY KEY,
  job_id TEXT NOT NULL,
  quiverlab_version TEXT NOT NULL,
  created_at TEXT NOT NULL,
  last_hit_at TEXT NOT NULL,
  hits INTEGER NOT NULL DEFAULT 0
);
CREATE INDEX IF NOT EXISTS idx_result_cache_job ON result_cache(job_id);
CREATE INDEX IF NOT EXISTS idx_result_cache_recency ON result_cache(last_hit_at);

-- Cooperative cancellation (OFFLINE GUI only). A row here asks the worker parent
-- loop (webapp/worker/worker.py::run_one_job) to kill a RUNNING job's child
-- promptly. The ONLY writer is `request_cancel_for_ip`, called solely by the
-- offline-only cancel endpoint (webapp/server/offline.py); on the deployed server
-- nothing writes it, so the table stays empty and the per-tick `cancel_requested`
-- check is inert. A PENDING job needs no signal -- it is finalised in place. The
-- signal is torn down when the worker acts on it (`clear_cancel_request`).
CREATE TABLE IF NOT EXISTS cancel_requests (
  job_id TEXT PRIMARY KEY,
  requested_at TEXT NOT NULL
);
"""


class JobStore:
    """The whole app's queue + state, backed by one SQLite database in WAL mode.

    A fresh short-lived connection is opened per operation (see `_conn`): WAL
    allows many concurrent readers alongside a single writer, and per-call
    connections keep the store safe to share across worker processes without any
    Python-side locking. The web tier and the worker(s) each open their own
    connection; correctness comes from SQLite's transactions, not from a shared
    handle."""

    def __init__(self, db_path: Path):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)

    def _conn(self) -> sqlite3.Connection:
        """Open a WAL connection with autocommit + a long busy timeout.

        `isolation_level=None` puts the connection in autocommit mode so that
        `claim_next`/`consume_pending_big` can drive their own explicit
        `BEGIN IMMEDIATE` transactions; every other method is a single
        auto-committed statement. `busy_timeout` lets a writer wait instead of
        failing immediately when another process holds the write lock."""
        conn = sqlite3.connect(self.db_path, timeout=30, isolation_level=None)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA busy_timeout=30000")
        return conn

    def init_schema(self) -> None:
        """Create the `jobs`, `pending_big`, `feedback`, `result_cache`, and
        `cancel_requests` tables (idempotent -- CREATE ... IF NOT EXISTS, so it
        upgrades an older DB in place by adding only the missing tables on the next
        start)."""
        with self._conn() as conn:
            conn.executescript(_SCHEMA)

    def _row_to_job(self, row: sqlite3.Row) -> Job:
        """Hydrate a `jobs` row into a `Job`, decoding the JSON spec/progress."""
        return Job(
            id=row["id"], spec=json.loads(row["spec"]), status=row["status"],
            progress=json.loads(row["progress"]), created_at=row["created_at"],
            started_at=row["started_at"], finished_at=row["finished_at"],
            quiverlab_version=row["quiverlab_version"],
            artifact_dir=row["artifact_dir"], error=row["error"], ip=row["ip"],
            tier=row["tier"], email=row["email"], email_hash=row["email_hash"],
            wall_seconds=row["wall_seconds"], mem_bytes=row["mem_bytes"],
            lang=row["lang"])

    def create_job(self, spec: dict, ip: str, tier: str = "anonymous",
                   email: str | None = None, email_hash: str | None = None,
                   wall_seconds: int | None = None, mem_bytes: int | None = None,
                   lang: str = "en") -> str:
        """Insert a new `pending` job and return its freshly minted ULID id."""
        jid = str(ULID())
        with self._conn() as conn:
            conn.execute(
                "INSERT INTO jobs (id, spec, ip, tier, email, email_hash, "
                "wall_seconds, mem_bytes, lang) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (jid, json.dumps(spec), ip, tier, email, email_hash,
                 wall_seconds, mem_bytes, lang))
        return jid

    def get_job(self, job_id: str) -> Job | None:
        """Return the `Job` with this id, or None if it does not exist."""
        with self._conn() as conn:
            row = conn.execute("SELECT * FROM jobs WHERE id=?", (job_id,)).fetchone()
        return self._row_to_job(row) if row else None

    def claim_next(self) -> Job | None:
        """Atomically claim the oldest pending job, flipping it to `running`.

        The read-then-update runs inside a single `BEGIN IMMEDIATE` transaction,
        which takes SQLite's write lock up front, so two concurrent workers can
        never select and claim the same row. Returns the claimed `Job`, or None
        when no job is pending."""
        conn = self._conn()
        try:
            conn.execute("BEGIN IMMEDIATE")
            row = conn.execute(
                "SELECT id FROM jobs WHERE status='pending' "
                "ORDER BY created_at LIMIT 1").fetchone()
            if row is None:
                conn.execute("COMMIT")
                return None
            conn.execute(
                "UPDATE jobs SET status='running', "
                "started_at=strftime('%Y-%m-%dT%H:%M:%SZ','now') WHERE id=?",
                (row["id"],))
            conn.execute("COMMIT")
            return self.get_job(row["id"])
        finally:
            conn.close()

    def requeue_stale_running(self) -> list[str]:
        """Requeue every `running` job back to `pending` (clearing started_at);
        return the requeued ids.

        Call this ONLY at worker-fleet startup -- the single-writer moment,
        before any poll loop begins claiming. `claim_next` only ever claims
        `pending`, so a worker that dies mid-job (`docker compose stop`, reboot,
        OOM-kill) strands its row in `running` FOREVER, permanently consuming the
        per-IP running slot (`count_running_for_ip` counts `running`). Running
        this once before the loops start adopts those orphans back into the
        queue. It is NOT safe to call while loops are live: it would yank a
        genuinely in-flight job's row out from under its worker, letting a second
        worker double-run it."""
        with self._conn() as conn:
            ids = [r["id"] for r in conn.execute(
                "SELECT id FROM jobs WHERE status='running'").fetchall()]
            conn.execute(
                "UPDATE jobs SET status='pending', started_at=NULL "
                "WHERE status='running'")
        return ids

    def update_progress(self, job_id: str, progress: dict) -> None:
        """Overwrite the job's progress blob (a small JSON-serialisable dict)."""
        with self._conn() as conn:
            conn.execute("UPDATE jobs SET progress=? WHERE id=?",
                         (json.dumps(progress), job_id))

    def mark_running(self, job_id: str) -> None:
        """Force a job to `running` and stamp its start time."""
        with self._conn() as conn:
            conn.execute(
                "UPDATE jobs SET status='running', "
                "started_at=strftime('%Y-%m-%dT%H:%M:%SZ','now') WHERE id=?",
                (job_id,))

    def mark_done(self, job_id: str, artifact_dir: str) -> None:
        """Mark a job `done`, recording its artifact dir, finish time, and the
        quiverlab version that produced it (imported lazily so the store module
        stays importable without quiverlab present)."""
        import quiverlab
        with self._conn() as conn:
            conn.execute(
                "UPDATE jobs SET status='done', artifact_dir=?, "
                "quiverlab_version=?, "
                "finished_at=strftime('%Y-%m-%dT%H:%M:%SZ','now') WHERE id=?",
                (artifact_dir, getattr(quiverlab, "__version__", "unknown"), job_id))

    def mark_failed(self, job_id: str, error: str) -> None:
        """Mark a job `failed`, storing the error string and finish time."""
        with self._conn() as conn:
            conn.execute(
                "UPDATE jobs SET status='failed', error=?, "
                "finished_at=strftime('%Y-%m-%dT%H:%M:%SZ','now') WHERE id=?",
                (error, job_id))

    # --- cooperative cancellation (offline GUI only) ------------------------
    def request_cancel_for_ip(self, ip: str, now_iso: str) -> str | None:
        """Cancel THIS caller's active (pending or running) job -- the offline
        GUI's "cancel the running job" action. Returns the cancelled job id, or
        None when the caller has nothing pending/running (an idempotent no-op).

        A PENDING job is finalised directly (no worker/child is involved) so the
        per-IP running slot frees immediately. A RUNNING job -- or a pending one
        the worker claims in the very same instant -- is instead SIGNALLED via
        `cancel_requests`; the worker parent loop notices within its ~1s poll tick,
        kills the child, and marks the job failed itself (so the child's late
        result can never overwrite the cancellation).

        Runs inside `BEGIN IMMEDIATE`, so it serialises against `claim_next`: the
        pending-vs-running decision can never tear, and a job claimed between the
        SELECT and the guarded UPDATE falls through to the signal path."""
        conn = self._conn()
        try:
            conn.execute("BEGIN IMMEDIATE")
            row = conn.execute(
                "SELECT id, status FROM jobs WHERE ip=? AND status IN "
                "('pending','running') ORDER BY created_at LIMIT 1", (ip,)).fetchone()
            if row is None:
                conn.execute("COMMIT")
                return None
            jid = row["id"]
            if row["status"] == "pending":
                cur = conn.execute(
                    "UPDATE jobs SET status='failed', error=?, finished_at=? "
                    "WHERE id=? AND status='pending'", (CANCEL_ERROR, now_iso, jid))
                if cur.rowcount:              # finalised before any worker claimed it
                    conn.execute("COMMIT")
                    return jid
            # running (or just-claimed): the worker loop must kill the child.
            conn.execute("INSERT INTO cancel_requests (job_id, requested_at) "
                         "VALUES (?, ?) ON CONFLICT(job_id) DO NOTHING", (jid, now_iso))
            conn.execute("COMMIT")
            return jid
        finally:
            conn.close()

    def cancel_requested(self, job_id: str) -> bool:
        """Whether a cooperative-cancel signal is pending for this job. The worker
        parent loop polls this each tick for a running job (offline GUI only; the
        table is always empty on the deployed server)."""
        with self._conn() as conn:
            return conn.execute("SELECT 1 FROM cancel_requests WHERE job_id=?",
                                (job_id,)).fetchone() is not None

    def clear_cancel_request(self, job_id: str) -> None:
        """Drop a job's cancel signal once the worker has acted on it."""
        with self._conn() as conn:
            conn.execute("DELETE FROM cancel_requests WHERE job_id=?", (job_id,))

    def count_running_for_ip(self, ip: str) -> int:
        """Count this ip's jobs that are still pending or running (its live load)."""
        with self._conn() as conn:
            return conn.execute(
                "SELECT COUNT(*) c FROM jobs WHERE ip=? AND status IN "
                "('pending','running')", (ip,)).fetchone()["c"]

    def count_created_today_for_ip(self, ip: str, now_iso: str) -> int:
        """Count this ip's jobs created on the same UTC day as `now_iso`."""
        day = now_iso[:10]
        with self._conn() as conn:
            return conn.execute(
                "SELECT COUNT(*) c FROM jobs WHERE ip=? AND substr(created_at,1,10)=?",
                (ip, day)).fetchone()["c"]

    def count_pending(self) -> int:
        """Count jobs waiting to be claimed (queue depth)."""
        with self._conn() as conn:
            return conn.execute(
                "SELECT COUNT(*) c FROM jobs WHERE status='pending'").fetchone()["c"]

    def purge_older_than(self, cutoff_iso: str) -> list[str]:
        """Delete finished jobs older than `cutoff_iso`; return their ids so the
        caller can clean up the matching artifact directories.

        A job referenced by a live `result_cache` row is EXCLUDED (the cache "pins"
        it): its row and artifacts must survive to back a cached replay. When the
        cache entry is later evicted (`cache_sweep`), the pin lifts and a subsequent
        retention pass reclaims the now-unreferenced job. This is the Plan-25
        artifact-lifetime rule: a finished job lives until it is BOTH older than the
        retention cutoff AND unreferenced by the cache."""
        with self._conn() as conn:
            rows = conn.execute(
                "SELECT id FROM jobs WHERE finished_at IS NOT NULL AND "
                "finished_at < ? AND id NOT IN (SELECT job_id FROM result_cache)",
                (cutoff_iso,)).fetchall()
            ids = [r["id"] for r in rows]
            conn.execute(
                "DELETE FROM jobs WHERE finished_at IS NOT NULL AND "
                "finished_at < ? AND id NOT IN (SELECT job_id FROM result_cache)",
                (cutoff_iso,))
        return ids

    # --- result cache (Plan 25) ---------------------------------------------
    def cache_get(self, key: str, now_iso: str) -> str | None:
        """Return the cached done-job id for `key` (a cache HIT), or None on a miss.

        A hit bumps the entry's `hits` and `last_hit_at=now_iso` (LRU recency). If
        the referenced job has vanished or is no longer `done` (should not happen
        while pinned -- defensive), the stale row is dropped and the call reports a
        miss. Non-transactional: a concurrent double-hit may lose one increment of
        the popularity counter, which is immaterial; the pin makes correctness
        independent of the counter."""
        with self._conn() as conn:
            row = conn.execute("SELECT job_id FROM result_cache WHERE key=?",
                               (key,)).fetchone()
            if row is None:
                return None
            job = conn.execute("SELECT status FROM jobs WHERE id=?",
                               (row["job_id"],)).fetchone()
            if job is None or job["status"] != "done":
                conn.execute("DELETE FROM result_cache WHERE key=?", (key,))
                return None
            conn.execute("UPDATE result_cache SET hits=hits+1, last_hit_at=? "
                         "WHERE key=?", (now_iso, key))
            return row["job_id"]

    def cache_put(self, key: str, job_id: str, version: str, now_iso: str) -> None:
        """Record a finished `job_id` under `key`. First writer wins: an identical
        request racing to compute the same math (benign -- results are exact and
        deterministic) records the same key on completion; the second `cache_put` is
        a no-op via ON CONFLICT DO NOTHING, so exactly one job stays pinned and the
        other ages out under normal retention. Never crashes on a duplicate key."""
        with self._conn() as conn:
            conn.execute(
                "INSERT INTO result_cache (key, job_id, quiverlab_version, "
                "created_at, last_hit_at, hits) VALUES (?, ?, ?, ?, ?, 0) "
                "ON CONFLICT(key) DO NOTHING",
                (key, job_id, version, now_iso, now_iso))

    def cache_row(self, key: str) -> dict | None:
        """Return the raw cache row (math-only: key, job_id, version, timestamps,
        hits) as a dict, or None. Used by tests and by callers that want the
        first-computed timestamp for the 'previously computed' UX."""
        with self._conn() as conn:
            row = conn.execute("SELECT * FROM result_cache WHERE key=?",
                               (key,)).fetchone()
        return dict(row) if row else None

    def cache_sweep(self, max_entries: int, current_version: str) -> int:
        """Evict cache rows and return how many were removed. Two passes:

        1. **Version purge** -- drop every row not from `current_version`
           (a library bump invalidates them; the key would already hide them, this
           reclaims the slots).
        2. **LRU size cap** -- keep the `max_entries` most-recently-hit rows, evict
           the rest (ties broken by created_at then key for determinism).

        Eviction only removes the fast-lookup row + the retention pin; it never
        deletes a job or its artifacts directly -- the next retention pass reclaims
        a now-unpinned job once it is also older than the cutoff."""
        removed = 0
        with self._conn() as conn:
            cur = conn.execute("DELETE FROM result_cache WHERE quiverlab_version != ?",
                               (current_version,))
            removed += cur.rowcount
            cur = conn.execute(
                "DELETE FROM result_cache WHERE key IN ("
                "  SELECT key FROM result_cache "
                "  ORDER BY last_hit_at DESC, created_at DESC, key DESC "
                "  LIMIT -1 OFFSET ?)", (max_entries,))
            removed += cur.rowcount
        return removed

    # --- feedback -----------------------------------------------------------
    def create_feedback(self, category: str, message: str, contact: str | None,
                        ip: str, job_ref: str | None, extra: str | None = None) -> str:
        """Insert a feedback row and return its ULID id.

        `extra` is a JSON string for category-specific fields (e.g. the
        `literature` category's reference + why-relevant). A single nullable
        column keeps the schema stable as new structured categories appear —
        cleaner than sparse per-category columns."""
        fid = str(ULID())
        with self._conn() as conn:
            conn.execute(
                "INSERT INTO feedback (id, category, message, contact, ip, job_ref, extra) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                (fid, category, message, contact, ip, job_ref, extra))
        return fid

    def count_feedback_today_for_ip(self, ip: str, now_iso: str) -> int:
        """Count feedback rows from this ip on the same UTC day as `now_iso`."""
        day = now_iso[:10]
        with self._conn() as conn:
            return conn.execute(
                "SELECT COUNT(*) c FROM feedback WHERE ip=? AND substr(created_at,1,10)=?",
                (ip, day)).fetchone()["c"]

    def list_feedback(self, limit: int = 200) -> list[dict]:
        """Return up to `limit` feedback rows, newest first, each including `extra`."""
        with self._conn() as conn:
            rows = conn.execute(
                "SELECT id, category, message, contact, created_at, job_ref, extra "
                "FROM feedback ORDER BY created_at DESC LIMIT ?", (limit,)).fetchall()
        return [dict(r) for r in rows]

    # --- big jobs (spec §17) ------------------------------------------------
    def create_pending_big(self, spec: dict, email: str, email_hash: str,
                           lang: str = "en") -> str:
        """Stash a big-job request awaiting email verification; return its id."""
        pid = str(ULID())
        with self._conn() as conn:
            conn.execute(
                "INSERT INTO pending_big (id, spec, email, email_hash, lang) "
                "VALUES (?, ?, ?, ?, ?)", (pid, json.dumps(spec), email, email_hash, lang))
        return pid

    def consume_pending_big(self, pid: str) -> dict | None:
        """Atomically fetch-and-delete a pending row: the magic link is
        single-use because the second consume finds nothing. Returns
        {spec, email, email_hash, lang} or None if already used/never existed."""
        conn = self._conn()
        try:
            conn.execute("BEGIN IMMEDIATE")
            row = conn.execute("SELECT spec, email, email_hash, lang FROM pending_big "
                               "WHERE id=?", (pid,)).fetchone()
            if row is None:
                conn.execute("COMMIT")
                return None
            conn.execute("DELETE FROM pending_big WHERE id=?", (pid,))
            conn.execute("COMMIT")
            return {"spec": json.loads(row["spec"]), "email": row["email"],
                    "email_hash": row["email_hash"], "lang": row["lang"]}
        finally:
            conn.close()

    def count_big_running_for_email_hash(self, email_hash: str) -> int:
        """Count this email-hash's big jobs still pending or running (concurrency cap)."""
        with self._conn() as conn:
            return conn.execute(
                "SELECT COUNT(*) c FROM jobs WHERE tier='big' AND email_hash=? "
                "AND status IN ('pending','running')", (email_hash,)).fetchone()["c"]

    def count_big_since_for_email_hash(self, email_hash: str, since_iso: str) -> int:
        """Count this email-hash's big jobs created at/after `since_iso` (rolling-window cap)."""
        with self._conn() as conn:
            return conn.execute(
                "SELECT COUNT(*) c FROM jobs WHERE tier='big' AND email_hash=? "
                "AND created_at >= ?", (email_hash, since_iso)).fetchone()["c"]

    def count_big_pending(self) -> int:
        """Count big-tier jobs waiting to be claimed."""
        with self._conn() as conn:
            return conn.execute(
                "SELECT COUNT(*) c FROM jobs WHERE tier='big' AND status='pending'"
                ).fetchone()["c"]

    def clear_email(self, job_id: str) -> None:
        """Delete the plaintext address; keep email_hash for weekly rate-limiting."""
        with self._conn() as conn:
            conn.execute("UPDATE jobs SET email=NULL WHERE id=?", (job_id,))

    def purge_pending_big(self, cutoff_iso: str) -> int:
        """Delete unverified big-job requests older than `cutoff_iso`; return the count."""
        with self._conn() as conn:
            cur = conn.execute("DELETE FROM pending_big WHERE created_at < ?",
                               (cutoff_iso,))
            return cur.rowcount
