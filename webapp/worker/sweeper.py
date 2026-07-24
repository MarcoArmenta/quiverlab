"""Retention sweeper: delete finished jobs and their artifact dirs older than
``cfg.retention_days``. ``now_iso`` is passed in explicitly (never an implicit
wall-clock read) so the cutoff is deterministic and testable."""
from __future__ import annotations

import shutil
from datetime import datetime, timedelta

from webapp.server.config import Config
from webapp.server.store import JobStore


def _cutoff(now_iso: str, days: int) -> str:
    now = datetime.strptime(now_iso, "%Y-%m-%dT%H:%M:%SZ")
    return (now - timedelta(days=days)).strftime("%Y-%m-%dT%H:%M:%SZ")


def sweep_once(store: JobStore, cfg: Config, now_iso: str) -> list[str]:
    """Purge DB rows finished before the retention cutoff and delete their
    artifact dirs; return the purged job ids."""
    cutoff = _cutoff(now_iso, cfg.retention_days)
    removed = store.purge_older_than(cutoff)
    for jid in removed:
        shutil.rmtree(cfg.artifacts_dir / jid, ignore_errors=True)
    return removed
