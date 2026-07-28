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
    """The curated manifest: every entry resolves to a valid ComputeRequest --
    inline entries directly, stored bundles via their request.json -- and the
    set spans quiver and family algebra kinds."""
    import json
    seed = _load()
    manifest = seed.load_manifest(_MANIFEST)
    assert 3 <= len(manifest) <= 10
    kinds = set()
    for entry in manifest:
        if set(entry) == {"stored"}:
            bundle = pathlib.Path(entry["stored"])
            assert bundle.is_dir(), f"stored bundle missing: {bundle}"
            spec = json.loads((bundle / "request.json").read_text(encoding="utf-8"))
            assert (bundle / "result.json").exists(), f"{bundle} has no result.json"
        else:
            spec = entry
        ComputeRequest.model_validate(spec)      # raises on any invalid example
        kinds.add(spec["algebra"]["kind"])
    assert "quiver" in kinds and "family" in kinds


def test_stored_bundle_results_match_their_requests():
    """Each stored result.json was produced for exactly the request beside it
    (same algebra spec) and carries a non-empty results block for every
    requested compute kind."""
    import json
    seed = _load()
    for entry in seed.load_manifest(_MANIFEST):
        if set(entry) != {"stored"}:
            continue
        bundle = pathlib.Path(entry["stored"])
        spec = json.loads((bundle / "request.json").read_text(encoding="utf-8"))
        result = json.loads((bundle / "result.json").read_text(encoding="utf-8"))
        assert result["algebra"] == spec["algebra"], bundle
        from quiverlab.hpc.spec import parse_compute_item
        want = {parse_compute_item(c).kind for c in spec["compute"]}
        assert set(result["results"]) == want, bundle


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


def test_seed_from_stored_bundle_copies_without_recompute(tmp_path):
    """A {stored: dir} entry seeds by COPYING the bundle: the cache row appears
    under the request's canonical key and the artifacts are the bundle's files
    verbatim -- no compute dispatch runs (the result carries a sentinel no
    runner would produce, and it survives)."""
    import json
    seed = _load()
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    (bundle / "request.json").write_text(json.dumps(LOOP_CARTAN))
    sentinel = {"algebra": LOOP_CARTAN["algebra"],
                "results": {"cartan": {"matrix": [[3]], "sentinel": "stored"}}}
    (bundle / "result.json").write_text(json.dumps(sentinel))
    (bundle / "report.html").write_text("<html>stored report</html>")
    out = tmp_path / "seed" / "seed-cache.db"
    ok, skipped, total = seed.seed([{"stored": str(bundle)}], out)
    assert (ok, skipped, total) == (1, 0, 1)
    store = JobStore(out)
    row = store.cache_row(_key(LOOP_CARTAN))
    assert row is not None
    art = out.parent / "artifacts" / row["job_id"]
    assert json.loads((art / "result.json").read_text())["results"]["cartan"][
        "sentinel"] == "stored"
    assert (art / "report.html").read_text() == "<html>stored report</html>"
    assert not (art / "request.json").exists()   # the request is not an artifact


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
