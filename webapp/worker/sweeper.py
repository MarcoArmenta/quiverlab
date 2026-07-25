"""Retention sweeper: delete finished jobs and their artifact dirs older than
``cfg.retention_days``. ``now_iso`` is passed in explicitly (never an implicit
wall-clock read) so the cutoff is deterministic and testable.

Alongside it (Plan 25) runs the result-cache sweep: a library-version purge plus
an LRU size cap over ``result_cache``. The two interact by design -- the cache
"pins" a finished job against retention (``purge_older_than`` skips cache-referenced
jobs), so the cache sweep runs FIRST: it lifts pins from stale-version / LRU-evicted
entries, and the retention pass in the same tick reclaims any now-unpinned job that
is also past the cutoff."""
from __future__ import annotations

import shutil
from datetime import datetime, timedelta

from webapp.server import cache
from webapp.server.config import Config
from webapp.server.store import JobStore


def _cutoff(now_iso: str, days: int) -> str:
    now = datetime.strptime(now_iso, "%Y-%m-%dT%H:%M:%SZ")
    return (now - timedelta(days=days)).strftime("%Y-%m-%dT%H:%M:%SZ")


def sweep_cache_once(store: JobStore, cfg: Config) -> int:
    """Version-purge + LRU size-cap the result cache; return rows evicted. A no-op
    when the cache is disabled. Eviction only removes fast-lookup rows and their
    retention pins -- artifacts are reclaimed by the retention pass, not here."""
    if not cfg.cache_enabled:
        return 0
    return store.cache_sweep(cfg.cache_max_entries, cache.library_version())


def sweep_once(store: JobStore, cfg: Config, now_iso: str) -> list[str]:
    """Purge DB rows finished before the retention cutoff and delete their
    artifact dirs; return the purged job ids. Cache-pinned jobs are left in place
    by ``purge_older_than`` (their artifacts back cached replays)."""
    cutoff = _cutoff(now_iso, cfg.retention_days)
    removed = store.purge_older_than(cutoff)
    for jid in removed:
        shutil.rmtree(cfg.artifacts_dir / jid, ignore_errors=True)
    return removed
