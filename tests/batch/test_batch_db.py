"""quiverlab.batch ResultsDB: SQLite persistence for scans (labdb lift)."""
import quiverlab.batch as batch
from quiverlab.batch.db import ResultsDB


def test_schema_and_upsert_idempotent():
    rec = {"name": "kA3", "builder": "linear_path_algebra", "args": [3], "N": 2,
           "dim": 6, "associative": True, "kind": "structural"}
    with ResultsDB(":memory:") as db:
        db.insert(rec)
        db.insert(dict(rec))                 # same (builder, args, N)
        assert len(db) == 1                  # INSERT OR REPLACE dedup


def test_query_roundtrip_and_json():
    import json
    rec = {"name": "qCI", "builder": "quantum_ci", "args": [3], "N": 4, "dim": 4,
           "associative": True, "asymmetric": True, "cx_homology": 1, "cx_cohomology": 0,
           "kind": "structural"}
    with ResultsDB(":memory:") as db:
        db.insert(rec)
        (got,) = db.query("builder = ?", ("quantum_ci",))
        assert got["name"] == "qCI"
        json.dumps(db.all())
        assert db.asymmetric()[0]["name"] == "qCI"
