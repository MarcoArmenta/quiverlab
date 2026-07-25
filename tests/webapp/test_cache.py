"""Plan 25 -- result-cache unit tests.

Two layers, both pure/fast:

  * the canonicalizer (``webapp/server/cache.py``): a deterministic key over the
    versioned compute request. Semantically identical requests MUST collide (dict
    key order is irrelevant; a JSON round-trip that turns arrow tuples into lists
    must not change the key); different mathematics -- or a different library
    version -- MUST NOT.
  * the store's ``result_cache`` table (``webapp/server/store.py``): put/get,
    LRU + version sweep, the retention "pin" (a cached job survives the ordinary
    job-retention purge), and the no-PII column contract.
"""
import sqlite3

from webapp.server.cache import canonical_blob, canonical_key
from webapp.server.schema import ComputeRequest
from webapp.server.store import JobStore


# --------------------------------------------------------------------------- #
# Canonicalizer
# --------------------------------------------------------------------------- #

_FAMILY = {"schema": 1,
           "algebra": {"kind": "family", "family": "QuantumCI",
                       "params": {"q": 1}, "field": {"kind": "GF", "p": 2, "n": 1}},
           "compute": ["hh_cohomology:0..3"],
           "artifacts": {"pdf": False, "tikz": False}}

_V = "0.1.0.dev0"


def test_key_is_deterministic():
    assert canonical_key(_FAMILY, _V) == canonical_key(_FAMILY, _V)
    assert len(canonical_key(_FAMILY, _V)) == 64          # sha256 hex


def test_dict_key_order_is_irrelevant():
    # params, field, and the top-level dict spelled in a different key order is the
    # SAME mathematics -> the SAME key (canonicalization sorts keys before hashing).
    reordered = {"artifacts": {"tikz": False, "pdf": False},
                 "compute": ["hh_cohomology:0..3"],
                 "algebra": {"field": {"n": 1, "p": 2, "kind": "GF"},
                             "params": {"q": 1}, "kind": "family",
                             "family": "QuantumCI"},
                 "schema": 1}
    assert canonical_key(reordered, _V) == canonical_key(_FAMILY, _V)


def test_tuple_and_list_arrows_collide():
    # QuiverAlgebraSpec.model_dump() yields arrow *tuples*; a JSON round-trip
    # through the store turns them into *lists*. The API checks the cache with the
    # in-memory (tuple) form and the worker records it from the reloaded (list)
    # form -- they MUST produce one key, so canonicalization goes through JSON.
    tup = {"schema": 1, "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False},
           "algebra": {"kind": "quiver", "vertices": [1, 2],
                       "arrows": {"a": (1, 2)}, "relations": [],
                       "field": {"kind": "GF", "p": 3, "n": 1}}}
    lst = {"schema": 1, "compute": ["cartan"], "artifacts": {"pdf": False, "tikz": False},
           "algebra": {"kind": "quiver", "vertices": [1, 2],
                       "arrows": {"a": [1, 2]}, "relations": [],
                       "field": {"kind": "GF", "p": 3, "n": 1}}}
    assert canonical_key(tup, _V) == canonical_key(lst, _V)


def test_library_version_is_part_of_the_key():
    # A version bump changes the key -> natural cache invalidation.
    assert canonical_key(_FAMILY, "0.1.0") != canonical_key(_FAMILY, "0.2.0")


def test_different_mathematics_do_not_collide():
    base = canonical_key(_FAMILY, _V)

    def mut(**over):
        d = {k: (dict(v) if isinstance(v, dict) else list(v) if isinstance(v, list) else v)
             for k, v in _FAMILY.items()}
        d.update(over)
        return canonical_key(d, _V)

    # different field prime
    alg = dict(_FAMILY["algebra"]); alg["field"] = {"kind": "GF", "p": 3, "n": 1}
    assert mut(algebra=alg) != base
    # different params
    alg2 = dict(_FAMILY["algebra"]); alg2["params"] = {"q": 2}
    assert mut(algebra=alg2) != base
    # different invariant / degree range
    assert mut(compute=["hh_cohomology:0..4"]) != base
    assert mut(compute=["cartan"]) != base
    # different artifact request (pdf) -> different produced files -> different key
    assert mut(artifacts={"pdf": True, "tikz": False}) != base
    # different schema version
    assert mut(schema=2) != base


def test_compute_order_is_significant():
    # Conservative decision (see cache.py docstring): the compute LIST order is
    # part of the key because it changes the produced artifacts (result.json key
    # order, the reproduce snippet, and -- with two HH kinds -- which one's PDF is
    # rendered). So a permuted compute list is a distinct cache entry, never a
    # subtly-different replay. Dict key order (above) still does not matter.
    a = dict(_FAMILY); a["compute"] = ["cartan", "dimension"]
    b = dict(_FAMILY); b["compute"] = ["dimension", "cartan"]
    assert canonical_key(a, _V) != canonical_key(b, _V)


def test_canonical_blob_is_sorted_json():
    blob = canonical_blob(_FAMILY, _V)
    assert blob.startswith('{"lib":')          # "lib" sorts before "req"
    assert " " not in blob                       # compact separators


# --------------------------------------------------------------------------- #
# Module requests (Plan 26): the module block canonicalizes deterministically
# THROUGH canonical_key, exactly like family/quiver requests. Specs go through
# the schema's model_dump (the normalized form the API and the worker both hash).
# --------------------------------------------------------------------------- #

_ALG = {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
        "relations": ["x*x*x"], "field": {"kind": "GF", "p": 2, "n": 1}}


def _mod_spec(**module):
    body = {"schema": 2, "algebra": _ALG, "compute": ["dimension_vector"],
            "artifacts": {"pdf": False, "tikz": False}, "module": module}
    return ComputeRequest.model_validate(body).model_dump(by_alias=True)


def test_same_module_typed_twice_collides():
    a = _mod_spec(dims={"1": 2}, maps={"x": [[0, 0], [1, 0]]})
    b = _mod_spec(dims={"1": 2}, maps={"x": [[0, 0], [1, 0]]})
    assert canonical_key(a, _V) == canonical_key(b, _V)


def test_module_side_omitted_equals_explicit_right():
    # The load-bearing side-default invariant: omitted side and side="right" are ONE
    # cache key (the schema always emits side="right"), so a user who never touches
    # the side toggle hits the same entry as one who set it to the default.
    omitted = _mod_spec(dims={"1": 2}, maps={"x": [[0, 0], [1, 0]]})
    explicit = _mod_spec(dims={"1": 2}, maps={"x": [[0, 0], [1, 0]]}, side="right")
    assert canonical_key(omitted, _V) == canonical_key(explicit, _V)


def test_module_dict_key_order_is_irrelevant():
    a = _mod_spec(dims={"1": 1, "2": 0}, maps={"x": [[0]]})
    # a spec spelled with the module dict keys in a different order is the SAME key
    b = _mod_spec(maps={"x": [[0]]}, side="right", dims={"2": 0, "1": 1})
    assert canonical_key(a, _V) == canonical_key(b, _V)


def test_builtin_nested_side_collides_with_lifted_side():
    # The builtin's nested `side` is lifted to the canonical top-level `side`, so the
    # task's `{"builtin": {..., "side": "left"}}` and the lifted form are one key.
    a = _mod_spec(builtin={"kind": "simple", "vertex": 1, "side": "left"})
    b = _mod_spec(builtin={"kind": "simple", "vertex": 1}, side="left")
    assert canonical_key(a, _V) == canonical_key(b, _V)


def test_different_modules_do_not_collide():
    base = canonical_key(_mod_spec(dims={"1": 2}, maps={"x": [[0, 0], [1, 0]]}), _V)
    assert canonical_key(_mod_spec(dims={"1": 3},
                                   maps={"x": [[0, 0, 0], [1, 0, 0], [0, 1, 0]]}), _V) != base
    assert canonical_key(_mod_spec(dims={"1": 2}, maps={"x": [[0, 0], [0, 0]]}), _V) != base
    assert canonical_key(_mod_spec(dims={"1": 2}, maps={"x": [[0, 0], [1, 0]]},
                                   side="left"), _V) != base
    assert canonical_key(_mod_spec(builtin={"kind": "simple", "vertex": 1}), _V) != base


def test_module_request_key_differs_from_same_algebra_without_module():
    with_mod = canonical_key(_mod_spec(dims={"1": 2}, maps={"x": [[0, 0], [1, 0]]}), _V)
    plain = ComputeRequest.model_validate(
        {"schema": 1, "algebra": _ALG, "compute": ["cartan"],
         "artifacts": {"pdf": False, "tikz": False}}).model_dump(by_alias=True)
    assert "module" not in plain                       # non-module dump is unchanged
    assert canonical_key(plain, _V) != with_mod


# --------------------------------------------------------------------------- #
# Store: result_cache table
# --------------------------------------------------------------------------- #

def _store(tmp_path):
    s = JobStore(tmp_path / "jobs.sqlite3")
    s.init_schema()
    return s


def _done_job(s, spec=None, ip="1.2.3.4"):
    jid = s.create_job(spec or {"compute": ["cartan"]}, ip=ip)
    s.claim_next()
    s.mark_done(jid, artifact_dir="/x/" + jid)
    return jid


def test_put_get_roundtrip(tmp_path):
    s = _store(tmp_path)
    jid = _done_job(s)
    assert s.cache_get("k1", "2026-01-01T00:00:00Z") is None      # miss
    s.cache_put("k1", jid, "v1", "2026-01-01T00:00:00Z")
    assert s.cache_get("k1", "2026-01-02T00:00:00Z") == jid       # hit


def test_get_bumps_hits_and_recency(tmp_path):
    s = _store(tmp_path)
    jid = _done_job(s)
    s.cache_put("k", jid, "v1", "2026-01-01T00:00:00Z")
    assert s.cache_row("k")["hits"] == 0
    s.cache_get("k", "2026-01-05T00:00:00Z")
    s.cache_get("k", "2026-01-06T00:00:00Z")
    row = s.cache_row("k")
    assert row["hits"] == 2
    assert row["last_hit_at"] == "2026-01-06T00:00:00Z"          # recency advanced


def test_put_is_idempotent_first_writer_wins(tmp_path):
    # The benign identical-request race: two jobs compute the same math and both
    # record the same key. The first put pins its job; the second is a no-op (never
    # a crash, never a repoint) -- so exactly one cache row, stably pinned.
    s = _store(tmp_path)
    j1 = _done_job(s)
    j2 = _done_job(s)
    s.cache_put("k", j1, "v1", "2026-01-01T00:00:00Z")
    s.cache_put("k", j2, "v1", "2026-01-02T00:00:00Z")           # no crash, no repoint
    assert s.cache_get("k", "2026-01-03T00:00:00Z") == j1


def test_get_drops_dangling_row(tmp_path):
    # Defensive: if the referenced job is somehow gone/not-done, get returns a miss
    # AND self-heals by dropping the stale cache row.
    s = _store(tmp_path)
    jid = s.create_job({"compute": ["cartan"]}, ip="i")          # pending, not done
    s.cache_put("k", jid, "v1", "2026-01-01T00:00:00Z")
    assert s.cache_get("k", "2026-01-02T00:00:00Z") is None      # not 'done' -> miss
    assert s.cache_row("k") is None                              # dropped


def test_cached_job_survives_retention_pin(tmp_path):
    # The lifetime rule: the ordinary job-retention purge must NOT delete a job
    # that a live cache row references -- its artifacts must survive to be replayed.
    s = _store(tmp_path)
    cached = _done_job(s, ip="a")
    plain = _done_job(s, ip="b")
    s.cache_put("k", cached, "v1", "2026-01-01T00:00:00Z")
    removed = s.purge_older_than("2999-01-01T00:00:00Z")         # everything is "old"
    assert plain in removed                                      # unpinned -> reclaimed
    assert cached not in removed                                 # pinned -> survives
    assert s.get_job(cached) is not None
    assert s.get_job(plain) is None


def test_cache_sweep_lru_keeps_most_recent(tmp_path):
    s = _store(tmp_path)
    for i in range(5):
        jid = _done_job(s)
        # distinct, increasing last_hit_at so recency order is unambiguous
        s.cache_put(f"k{i}", jid, "v1", f"2026-01-0{i + 1}T00:00:00Z")
    removed = s.cache_sweep(max_entries=3, current_version="v1")
    assert removed == 2
    # The two least-recently-hit (k0, k1) are evicted; the three newest remain.
    assert s.cache_row("k0") is None and s.cache_row("k1") is None
    assert all(s.cache_row(f"k{i}") is not None for i in (2, 3, 4))


def test_cache_sweep_purges_other_versions(tmp_path):
    s = _store(tmp_path)
    old = _done_job(s)
    new = _done_job(s)
    s.cache_put("old", old, "0.1.0", "2026-01-01T00:00:00Z")
    s.cache_put("new", new, "0.2.0", "2026-01-02T00:00:00Z")
    removed = s.cache_sweep(max_entries=100, current_version="0.2.0")
    assert removed == 1
    assert s.cache_row("old") is None                           # stale version purged
    assert s.cache_row("new") is not None


def test_cache_row_is_math_only_no_pii(tmp_path):
    # Privacy: the cache row carries mathematics only -- key, job ref, version,
    # timestamps, hit counter. No email/ip/token columns, and no value that could
    # identify the requester.
    s = _store(tmp_path)
    jid = _done_job(s)
    s.cache_put("k", jid, "v1", "2026-01-01T00:00:00Z")
    conn = sqlite3.connect(s.db_path)
    try:
        cols = {r[1] for r in conn.execute("PRAGMA table_info(result_cache)").fetchall()}
    finally:
        conn.close()
    forbidden = {"email", "email_hash", "ip", "token", "contact", "lang"}
    assert cols & forbidden == set(), f"cache table leaks PII columns: {cols & forbidden}"
    assert cols == {"key", "job_id", "quiverlab_version", "created_at",
                    "last_hit_at", "hits"}
