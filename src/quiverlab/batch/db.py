"""SQLite results database for batch scans (quiverlab.batch).

Ported verbatim (schema, migration, queries) from hanlab's ``labdb.ResultsDB``
-- Marco Armenta's ``HansConjecture`` codebase (MIT), the "B4" results-database
component.  Only the module docstring and the surrounding package layout are
quiverlab's; the SQLite schema, the ``_migrate_open_columns`` back-compat path,
``_nullable_bool`` and every query method are the bank's, kept faithful so a
quiverlab scan DB round-trips byte-for-byte with a hanlab one.

``ResultsDB`` persists ``scan.analyze`` records to SQLite, indexed by complexity
and structural data, and exposes the counterexample-search queries the Han
programme needs -- in particular ``bounded_homology_candidates()`` ("bounded HH
but apparently infinite global dimension", the shape of a Han counterexample) and
``asymmetric()`` ("HH_* and HH^* have different complexity").  Records are stored
as exact JSON (ints/strings; no floats): ``record_json`` is the full record and
the promoted columns are exact copies of its scalar fields.
"""

import json
import sqlite3

_SCHEMA = """
CREATE TABLE IF NOT EXISTS results (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    name          TEXT,
    builder       TEXT,
    args          TEXT,
    dim           INTEGER,
    N             INTEGER,
    associative   INTEGER,
    cx_homology   TEXT,
    cx_cohomology TEXT,
    asymmetric    INTEGER,
    homology_bounded      INTEGER,
    homology_nonzero_tail INTEGER,
    kind          TEXT,                 -- 'open' for reduction_system rows, else NULL
    han_verdict   TEXT,                 -- open zone: nonzero-tail / COUNTEREXAMPLE? / ...
    growth        TEXT,                 -- open zone: bounded / growing / ?
    depth_reached INTEGER,              -- open zone: last HH_* degree computed
    truncated_at  INTEGER,              -- open zone: degree the budget stopped at, or NULL
    record_json   TEXT,
    UNIQUE(builder, args, N)
);
CREATE INDEX IF NOT EXISTS idx_cx ON results(cx_homology, cx_cohomology);
CREATE INDEX IF NOT EXISTS idx_asym ON results(asymmetric);
CREATE INDEX IF NOT EXISTS idx_verdict ON results(han_verdict);
"""


def _nullable_bool(v):
    """None -> None (SQL NULL); else 0/1.  Open-zone records have asymmetric=None,
    which must not be lossily collapsed to 0 ('symmetric')."""
    return None if v is None else int(bool(v))


class ResultsDB:
    """SQLite store for scan records, with the counterexample-search queries."""

    def __init__(self, path=":memory:"):
        self.conn = sqlite3.connect(path)
        self.conn.row_factory = sqlite3.Row
        self.conn.executescript(_SCHEMA)
        self._migrate_open_columns()        # backward-compat for old on-disk DBs
        self.conn.commit()

    def _migrate_open_columns(self):
        """Add the open-zone columns to a pre-existing table that lacks them.

        No-op on a DB created with the current _SCHEMA -- the supported path always
        builds the DB fresh, so this is belt-and-suspenders for ad-hoc reuse of a
        stale on-disk file (its pre-existing rows backfill as NULL)."""
        have = {r["name"] for r in self.conn.execute("PRAGMA table_info(results)")}
        for col, decl in (("kind", "TEXT"), ("han_verdict", "TEXT"),
                          ("growth", "TEXT"), ("depth_reached", "INTEGER"),
                          ("truncated_at", "INTEGER")):
            if col not in have:
                self.conn.execute(f"ALTER TABLE results ADD COLUMN {col} {decl}")

    def close(self):
        self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    def insert(self, rec):
        """Insert (or replace) one analyze() record (structural or open-zone).

        The five open-zone columns are read with `.get(...)`, so a structural record
        (which lacks them) writes NULL -- backward compatible."""
        self.conn.execute(
            """INSERT OR REPLACE INTO results
               (name, builder, args, dim, N, associative, cx_homology, cx_cohomology,
                asymmetric, homology_bounded, homology_nonzero_tail,
                kind, han_verdict, growth, depth_reached, truncated_at, record_json)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (
                rec.get("name"),
                rec.get("builder"),
                json.dumps(rec.get("args", [])),
                rec.get("dim"),
                rec.get("N"),
                int(rec.get("associative", False)),
                str(rec.get("cx_homology")),
                str(rec.get("cx_cohomology")),
                _nullable_bool(rec.get("asymmetric")),    # open: None -> SQL NULL
                int(rec.get("homology_bounded", False)),
                int(rec.get("homology_nonzero_tail", False)),
                rec.get("kind"),                          # 'open' or None
                rec.get("han_verdict"),                   # str or None
                rec.get("growth"),                        # str or None
                rec.get("depth_reached"),                 # int or None
                rec.get("truncated_at"),                  # int or None
                json.dumps(rec),
            ),
        )
        self.conn.commit()

    def insert_many(self, recs):
        for r in recs:
            self.insert(r)

    def all(self):
        return [json.loads(r["record_json"])
                for r in self.conn.execute("SELECT record_json FROM results")]

    def query(self, where, params=()):
        rows = self.conn.execute(
            f"SELECT record_json FROM results WHERE {where}", params)
        return [json.loads(r["record_json"]) for r in rows]

    def asymmetric(self):
        """Algebras whose HH_* and HH^* have different apparent complexity."""
        return self.query("asymmetric = 1")

    def bounded_homology_candidates(self):
        """Counterexample-shaped rows (associative): (i) bounded HH_* -- not obviously
        gl.dim finite (structural + open), OR (ii) open-zone rows whose HH_* vanished in
        high degree (han_verdict='COUNTEREXAMPLE?').  `han_verdict` is NULL on structural
        rows, so the second disjunct never changes their behavior."""
        return self.query(
            "associative = 1 AND "
            "(homology_bounded = 1 OR han_verdict = 'COUNTEREXAMPLE?')"
        )

    def open_counterexample_candidates(self):
        """Open-zone rows whose HH_* appears to vanish in high degree (the literal
        Han-counterexample shape); indexed by idx_verdict."""
        return self.query("han_verdict = 'COUNTEREXAMPLE?'")

    def by_complexity(self, cx_homology):
        return self.query("cx_homology = ?", (str(cx_homology),))

    def __len__(self):
        return self.conn.execute("SELECT COUNT(*) FROM results").fetchone()[0]
