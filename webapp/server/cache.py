"""Plan 25 -- the result cache: a deterministic key over a compute request, plus
the thin lookup/record helpers the three tiers share.

Exact computations are DETERMINISTIC: the same canonical request always yields the
same result. So a finished result can be replayed for any later identical request
-- across users, across tiers -- with zero recompute and zero loss of correctness.

The key is ``sha256`` of a canonical JSON encoding of ``{"lib": <library version>,
"req": <the versioned compute-request spec>}``. Canonicalisation:

  * **Sorted keys** -- dict key order (``params``, ``arrows``, ``field`` blocks) is
    irrelevant; the same mathematics always hashes identically.
  * **Through JSON** -- ``model_dump()`` yields arrow *tuples*, a store round-trip
    yields *lists*; ``json.dumps`` erases the distinction, so the API (in-memory
    tuple form) and the worker (reloaded list form) compute ONE key.
  * **Library version in the key** -- a version bump changes every key, so old
    entries can never be replayed for a new build (natural invalidation); the
    schema version rides along inside the request spec, so it invalidates too.
  * **Compute-list order is significant** (a deliberately conservative choice): the
    order changes the produced artifacts (``result.json`` key order, the
    ``reproduce`` snippet, and -- with two HH kinds -- which one's worked-steps PDF
    is rendered). Treating a permuted list as a distinct entry guarantees a cache
    hit replays EXACTLY what this request would have produced, never a subtly
    different artifact. Dict key order is the only thing normalised away.
  * **Artifact flags are in the key** -- a ``pdf: true`` request produces a
    ``trace.pdf`` a ``pdf: false`` run does not, so they are distinct entries.

Nothing user-identifying enters the key: it is a hash of mathematics + version
only. Big-job requests key on the spec with ``email``/``lang`` stripped (the same
stripped spec the store persists), so a big example and an anonymous queued
example of the same mathematics share one cache entry.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone


def canonical_blob(spec: dict, library_version: str) -> str:
    """The canonical JSON string hashed into the cache key. Pure and deterministic:
    sorted keys, compact separators, JSON-normalised (tuples == lists)."""
    payload = {"lib": library_version, "req": spec}
    return json.dumps(payload, separators=(",", ":"), sort_keys=True, default=str)


def canonical_key(spec: dict, library_version: str) -> str:
    """The cache key: ``sha256`` hex of :func:`canonical_blob`. Semantically
    identical requests collide; different mathematics (or library versions) do not."""
    return hashlib.sha256(canonical_blob(spec, library_version).encode("utf-8")).hexdigest()


def library_version() -> str:
    """The running quiverlab version, part of every key (imported lazily so this
    module stays importable without quiverlab present, mirroring ``store.py``)."""
    import quiverlab
    return getattr(quiverlab, "__version__", "unknown")


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# --------------------------------------------------------------------------- #
# Tier-shared helpers (all three compute tiers check the cache FIRST; the worker
# records the key on completion). Config gates the whole feature.
# --------------------------------------------------------------------------- #

def lookup(store, cfg, spec: dict, now_iso: str | None = None) -> str | None:
    """Return the cached done-job id for this request spec, or None on a miss (or
    when the cache is disabled). Bumps the entry's LRU recency + hit counter."""
    if not getattr(cfg, "cache_enabled", True):
        return None
    key = canonical_key(spec, library_version())
    return store.cache_get(key, now_iso or _now_iso())


def record(store, cfg, job, now_iso: str | None = None) -> None:
    """Record a freshly-finished job in the cache so any later identical request
    replays it. No-op when the cache is disabled. Called by the worker AFTER the
    job is marked done -- a failure here must never lose the result, so callers
    wrap this in a try/except."""
    if not getattr(cfg, "cache_enabled", True):
        return
    v = library_version()
    store.cache_put(canonical_key(job.spec, v), job.id, v, now_iso or _now_iso())
