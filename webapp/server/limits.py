"""Rate-limit helpers over the salted-IP-hash counters in the store. Each helper
returns a human refusal message (surfaced as a 429) or None when the request is
within limits. The ``ip_hash`` passed in is already the salted SHA-256 -- the raw
address never reaches this layer or the store."""
from __future__ import annotations

from webapp.server.config import Config
from webapp.server.store import JobStore


def check_can_queue(store: JobStore, cfg: Config, ip_hash: str, now_iso: str) -> str | None:
    """None if a new anonymous job may be queued, else a refusal message."""
    if store.count_pending() >= cfg.global_queue_max:
        return ("The queue is full right now. Please retry shortly, or run "
                "locally: pip install quiverlab")
    if store.count_running_for_ip(ip_hash) >= cfg.per_ip_running_max:
        return "You already have a job running. Please wait for it to finish."
    if store.count_created_today_for_ip(ip_hash, now_iso) >= cfg.per_ip_daily_max:
        return ("Daily job limit reached for your address. For heavier use, run "
                "locally: pip install quiverlab")
    return None


def check_feedback_allowed(store: JobStore, cfg: Config, ip_hash: str,
                           now_iso: str) -> str | None:
    """None if this address may submit feedback today, else a refusal message."""
    if store.count_feedback_today_for_ip(ip_hash, now_iso) >= cfg.feedback_daily_max:
        return "Daily feedback limit reached for your address. Please try again tomorrow."
    return None
