"""Plan 28 -- build-time cache seeding tests (fast).

Validates the shipped ``webapp/precomputed/manifest.yaml`` against the request
schema, and drives ``container/seed_cache.py`` end-to-end on a tiny manifest:
the cache keys it writes must equal ``canonical_key`` of the requests, it must be
idempotent on a rerun, and its ``result_cache`` rows must carry mathematics only
(no email/ip/token -- Plan-25 invariant).
"""
import importlib.util
import pathlib
import sqlite3

from webapp.server.cache import canonical_key, library_version
from webapp.server.schema import ComputeRequest
from webapp.server.store import JobStore

_ROOT = pathlib.Path(__file__).resolve().parents[2]
_SEED_CACHE_PY = _ROOT / "container" / "seed_cache.py"
_MANIFEST = _ROOT / "webapp" / "precomputed" / "manifest.yaml"

LOOP_CARTAN = {
    "schema": 1,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
    "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False},
}
MODULE_BUILTIN = {
    "schema": 2,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}},
    "module": {"builtin": {"kind": "simple", "vertex": 1}, "side": "right"},
    "compute": ["dimension_vector"], "artifacts": {"pdf": False, "tikz": False},
}


def _load():
    spec = importlib.util.spec_from_file_location("seed_cache", _SEED_CACHE_PY)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _key(spec):
    dump = ComputeRequest.model_validate(spec).model_dump(by_alias=True)
    return canonical_key(dump, library_version())


# --------------------------------------------------------------------------- #
# The shipped manifest
# --------------------------------------------------------------------------- #

def test_shipped_manifest_is_schema_valid():
    seed = _load()
    manifest = seed.load_manifest(_MANIFEST)
    assert 3 <= len(manifest) <= 5, "plan asks for 3-5 placeholder examples"
    for i, spec in enumerate(manifest):
        ComputeRequest.model_validate(spec)      # raises on any invalid example
    kinds = {spec.get("kind") or spec["algebra"]["kind"] for spec in manifest}
    assert "quiver" in kinds                     # spans quiver + family + module


def test_manifest_marked_placeholder():
    text = _MANIFEST.read_text(encoding="utf-8")
    assert "placeholder -- Marco to curate" in text


# --------------------------------------------------------------------------- #
# Seeding end to end
# --------------------------------------------------------------------------- #

def test_seed_writes_rows_matching_canonical_keys(tmp_path):
    seed = _load()
    out = tmp_path / "seed" / "seed-cache.db"
    manifest = [LOOP_CARTAN, MODULE_BUILTIN]
    ok, skipped, total = seed.seed(manifest, out)
    assert (ok, skipped, total) == (2, 0, 2)

    store = JobStore(out)
    for spec in manifest:
        row = store.cache_row(_key(spec))
        assert row is not None, "seed did not record the request under its canonical key"
        job = store.get_job(row["job_id"])
        assert job is not None and job.status == "done"
        # artifacts live in the sibling artifacts/ tree (the offline copy bundle)
        assert (out.parent / "artifacts" / row["job_id"] / "result.json").exists()
        assert job.ip == ""                       # no PII on the backing job row


def test_seed_is_idempotent_on_rerun(tmp_path):
    seed = _load()
    out = tmp_path / "seed" / "seed-cache.db"
    ok1, skipped1, _ = seed.seed([LOOP_CARTAN], out)
    assert (ok1, skipped1) == (1, 0)
    ok2, skipped2, _ = seed.seed([LOOP_CARTAN], out)      # already seeded
    assert (ok2, skipped2) == (0, 1)
    conn = sqlite3.connect(out)
    try:
        n = conn.execute("SELECT COUNT(*) FROM result_cache").fetchone()[0]
    finally:
        conn.close()
    assert n == 1                                          # exactly one row, not two


def test_seed_cache_rows_are_math_only(tmp_path):
    seed = _load()
    out = tmp_path / "seed" / "seed-cache.db"
    seed.seed([LOOP_CARTAN], out)
    conn = sqlite3.connect(out)
    try:
        cols = {r[1] for r in conn.execute(
            "PRAGMA table_info(result_cache)").fetchall()}
    finally:
        conn.close()
    assert cols == {"key", "job_id", "quiverlab_version", "created_at",
                    "last_hit_at", "hits"}
    assert cols & {"email", "email_hash", "ip", "token", "contact", "lang"} == set()


def test_seed_skips_failures_and_reports(tmp_path):
    seed = _load()
    out = tmp_path / "seed" / "seed-cache.db"
    bad = {"schema": 1,
           "algebra": {"kind": "family", "family": "NoSuchFamily", "params": {},
                       "field": {"kind": "GF", "p": 2, "n": 1}},
           "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False}}
    ok, skipped, total = seed.seed([bad, LOOP_CARTAN], out)
    assert (ok, total) == (1, 2)                  # bad skipped, good seeded


def test_main_exits_nonzero_when_all_fail(tmp_path):
    seed = _load()
    bad_manifest = tmp_path / "bad.yaml"
    bad_manifest.write_text(
        "examples:\n"
        "  - schema: 1\n"
        "    algebra: {kind: family, family: NoSuchFamily, params: {},\n"
        "              field: {kind: GF, p: 2, n: 1}}\n"
        "    compute: [cartan]\n"
        "    artifacts: {pdf: false, tikz: false}\n", encoding="utf-8")
    rc = seed.main(["--manifest", str(bad_manifest),
                    "--out", str(tmp_path / "seed.db")])
    assert rc == 1


def test_main_succeeds_on_a_good_manifest(tmp_path):
    seed = _load()
    good = tmp_path / "good.yaml"
    good.write_text(
        "examples:\n"
        "  - schema: 1\n"
        "    algebra: {kind: quiver, vertices: [1], arrows: {x: [1, 1]},\n"
        "              relations: ['x*x*x'], field: {kind: GF, p: 2, n: 1}}\n"
        "    compute: [cartan]\n"
        "    artifacts: {pdf: false, tikz: false}\n", encoding="utf-8")
    rc = seed.main(["--manifest", str(good), "--out", str(tmp_path / "seed.db")])
    assert rc == 0
