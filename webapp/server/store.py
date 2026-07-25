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
        """Create the `jobs`, `pending_big`, and `feedback` tables (idempotent)."""
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
        caller can clean up the matching artifact directories."""
        with self._conn() as conn:
            rows = conn.execute(
                "SELECT id FROM jobs WHERE finished_at IS NOT NULL AND "
                "finished_at < ?", (cutoff_iso,)).fetchall()
            ids = [r["id"] for r in rows]
            conn.execute("DELETE FROM jobs WHERE finished_at IS NOT NULL AND "
                         "finished_at < ?", (cutoff_iso,))
        return ids

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
