#!/usr/bin/env python
"""Plan 28 -- build-time seeding of the offline GUI's result cache.

Reads a manifest of compute requests (``webapp/precomputed/manifest.yaml``), runs
each one through the SAME dispatch the webapp uses (``webapp.server.runner.run_spec``),
and records the finished result in a FRESH SQLite result cache at ``--out``. The
image ships that DB (plus its sibling ``artifacts/`` tree); ``webapp.server.offline``
copies the bundle into the user data dir on first run, so the offline GUI replays
the examples instantly (Plan-25 cache) with zero recompute.

Invariants:
  * math-only cache rows -- no email/ip/token ever enters ``result_cache``
    (Plan 25); the backing job rows carry an EMPTY ip, never a real address.
  * idempotent -- an already-seeded request (its canonical key present) is skipped,
    so rerunning into the same DB is a no-op.
  * resilient -- a per-example failure is warned and skipped; the script exits
    non-zero only if EVERY example failed (or the manifest was empty).

Runnable standalone::

    python container/seed_cache.py --manifest webapp/precomputed/manifest.yaml \\
        --out /build/seed/seed-cache.db

The seed bundle is ``<out>`` (the DB) with artifacts under ``<out>.parent/artifacts/``.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Runnable standalone from anywhere: put the repo root (this file lives in
# ``container/``) on sys.path so ``import webapp...`` resolves without a PYTHONPATH.
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

_log = logging.getLogger("quiverlab_seed_cache")


def load_manifest(path) -> list[dict]:
    """Parse the YAML manifest into a list of entries. Accepts either a
    top-level ``{examples: [...]}`` mapping or a bare list. An entry is EITHER a
    full compute request (computed at seed time) OR a stored-bundle reference
    ``{stored: <dir>}`` (relative to the manifest file) whose directory carries a
    precomputed ``request.json`` + ``result.json`` (+ extra artifacts) -- the
    zero-recompute path for the big curated examples."""
    import yaml

    path = Path(path)
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if isinstance(data, dict):
        data = data.get("examples", [])
    if not isinstance(data, list):
        raise ValueError(f"manifest {path!r} must be a list (or an 'examples:' list)")
    out = []
    for entry in data:
        if isinstance(entry, dict) and set(entry) == {"stored"}:
            # resolve now, against the manifest's own directory
            entry = {"stored": str((path.parent / entry["stored"]).resolve())}
        out.append(entry)
    return out


def _now_iso() -> str:
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def seed(manifest: list[dict], out_db) -> tuple[int, int, int]:
    """Seed ``out_db`` from ``manifest``. Returns ``(ok, skipped, total)``.

    A fresh (or existing) SQLite store is opened at ``out_db``; artifacts land in
    the sibling ``out_db.parent/artifacts/<job_id>/`` tree -- exactly the bundle
    shape ``webapp.server.offline.seed_first_run`` copies. Each request is
    validated, computed via ``run_spec``, marked done, and pinned in the cache
    under its canonical key."""
    # Lazy imports: keep --help import-light and pull no server (fastapi/uvicorn).
    from webapp.server.cache import canonical_key, library_version
    from webapp.server.runner import run_spec
    from webapp.server.schema import ComputeRequest
    from webapp.server.store import JobStore

    out_db = Path(out_db)
    out_db.parent.mkdir(parents=True, exist_ok=True)
    artifacts_dir = out_db.parent / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    store = JobStore(out_db)
    store.init_schema()
    version = library_version()

    ok = skipped = 0
    total = len(manifest)
    for i, spec in enumerate(manifest):
        label = f"example[{i}]"
        try:
            bundle = None
            if isinstance(spec, dict) and set(spec) == {"stored"}:
                # stored bundle: the request travels beside its precomputed result
                import json as _json
                bundle = Path(spec["stored"])
                label = f"example[{i}] ({bundle.name})"
                spec = _json.loads((bundle / "request.json").read_text(encoding="utf-8"))
            req = ComputeRequest.model_validate(spec)
            spec_dump = req.model_dump(by_alias=True)
            key = canonical_key(spec_dump, version)
            if store.cache_row(key) is not None:      # already seeded -> idempotent
                skipped += 1
                _log.info("%s: already cached, skipping", label)
                continue
            jid = store.create_job(spec_dump, ip="")  # no PII: empty ip
            store.mark_running(jid)
            art = artifacts_dir / jid
            art.mkdir(parents=True, exist_ok=True)
            if bundle is not None:
                # copy the precomputed artifacts (everything but the request);
                # then a loud sanity parse: the stored result must be a result.
                import json as _json
                import shutil as _shutil
                for f in sorted(bundle.iterdir()):
                    if f.name == "request.json" or f.is_dir():
                        continue
                    _shutil.copy2(f, art / f.name)
                stored_result = _json.loads(
                    (art / "result.json").read_text(encoding="utf-8"))
                if "results" not in stored_result:
                    raise ValueError(f"stored bundle {bundle} result.json carries "
                                     "no 'results' block")
            else:
                run_spec(req, art)                    # writes result.json + artifacts
            store.mark_done(jid, str(art))
            store.cache_put(key, jid, version, _now_iso())
            ok += 1
            _log.info("%s: seeded (job %s)", label, jid)
        except Exception as exc:                      # skip + warn; never abort the batch
            _log.warning("%s: failed (%s: %s) -- skipping", label,
                         type(exc).__name__, exc)
    # Collapse the WAL into the single .db file so the shipped seed is a complete,
    # self-contained file: the offline first-run copy takes ONLY the .db (not the
    # -wal/-shm sidecars), so any WAL-buffered rows must be folded in first.
    import sqlite3
    conn = sqlite3.connect(out_db)
    try:
        conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
    finally:
        conn.close()
    return ok, skipped, total


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description="Seed the offline GUI result cache.")
    ap.add_argument("--manifest", required=True, help="path to manifest.yaml")
    ap.add_argument("--out", required=True, help="output SQLite seed cache path")
    args = ap.parse_args(argv)

    manifest = load_manifest(args.manifest)
    ok, skipped, total = seed(manifest, args.out)
    _log.info("seed complete: %d ok, %d skipped, %d total -> %s",
              ok, skipped, total, args.out)
    # Exit non-zero only if NOTHING was seeded AND nothing was already present
    # (every example failed, or the manifest was empty).
    if ok == 0 and skipped == 0:
        _log.error("no examples seeded")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
