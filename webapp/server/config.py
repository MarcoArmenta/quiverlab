"""Environment-driven configuration for quiverlab-web. All limits are
overridable via QLWEB_* env vars so the DRAC deployment can be tuned without
code changes. Values here are wall-clock/byte/count knobs, not algebra."""
from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Mapping


def _int(env: Mapping[str, str], key: str, default: int) -> int:
    raw = env.get(key)
    return int(raw) if raw is not None and raw != "" else default


@dataclass(frozen=True)
class Config:
    data_dir: Path
    db_path: Path
    artifacts_dir: Path
    instant_wall_seconds: int
    instant_ops_threshold: int
    instant_max_degree: int
    job_wall_seconds: int
    job_mem_bytes: int
    result_max_bytes: int
    worker_processes: int
    per_ip_running_max: int
    per_ip_daily_max: int
    global_queue_max: int
    big_queue_max: int
    per_email_running_max: int
    per_email_weekly_max: int
    retention_days: int
    ip_hash_salt: str
    feedback_daily_max: int
    admin_token: str
    queued_ops_threshold: int
    queued_max_degree: int
    big_ops_threshold: int
    big_max_degree: int
    big_job_wall_seconds: int
    big_job_mem_bytes: int
    smtp_host: str
    smtp_port: int
    smtp_user: str
    smtp_pass: str
    smtp_from: str
    token_secret: str
    email_hash_salt: str
    big_token_ttl_seconds: int
    public_base_url: str
    docs_url: str
    cache_enabled: bool
    cache_max_entries: int

    @property
    def big_jobs_enabled(self) -> bool:
        # The big-job tier requires an outbound SMTP relay (spec §17); with no
        # relay configured the tier is disabled and the app says "run locally".
        return bool(self.smtp_host and self.smtp_from)

    @classmethod
    def from_env(cls, env: Mapping[str, str] | None = None) -> "Config":
        env = os.environ if env is None else env
        data_dir = Path(env.get("QLWEB_DATA_DIR", "/data"))
        cpu = os.cpu_count() or 2
        return cls(
            data_dir=data_dir,
            db_path=data_dir / "quiverlab_web.sqlite3",
            artifacts_dir=data_dir / "artifacts",
            instant_wall_seconds=_int(env, "QLWEB_INSTANT_WALL_SECONDS", 5),
            instant_ops_threshold=_int(env, "QLWEB_INSTANT_OPS_THRESHOLD", 2_000_000),
            instant_max_degree=_int(env, "QLWEB_INSTANT_MAX_DEGREE", 8),
            job_wall_seconds=_int(env, "QLWEB_JOB_WALL_SECONDS", 900),
            job_mem_bytes=_int(env, "QLWEB_JOB_MEM_BYTES", 4 * 1024 ** 3),
            result_max_bytes=_int(env, "QLWEB_RESULT_MAX_BYTES", 32 * 1024 ** 2),
            worker_processes=_int(env, "QLWEB_WORKER_PROCESSES", max(1, cpu - 2)),
            per_ip_running_max=_int(env, "QLWEB_PER_IP_RUNNING_MAX", 1),
            per_ip_daily_max=_int(env, "QLWEB_PER_IP_DAILY_MAX", 100),
            global_queue_max=_int(env, "QLWEB_GLOBAL_QUEUE_MAX", 200),
            big_queue_max=_int(env, "QLWEB_BIG_QUEUE_MAX", 20),
            per_email_running_max=_int(env, "QLWEB_PER_EMAIL_RUNNING_MAX", 1),
            per_email_weekly_max=_int(env, "QLWEB_PER_EMAIL_WEEKLY_MAX", 5),
            retention_days=_int(env, "QLWEB_RETENTION_DAYS", 90),
            # Salt for hashing client IPs before storage; MUST be set in prod
            # (a fixed default only keeps dev/tests deterministic).
            ip_hash_salt=env.get("QLWEB_IP_HASH_SALT", "quiverlab-dev-salt"),
            feedback_daily_max=_int(env, "QLWEB_FEEDBACK_DAILY_MAX", 5),
            # Empty by default: the /admin/feedback route is not registered
            # unless a token is set (see Task 12).
            admin_token=env.get("QLWEB_ADMIN_TOKEN", ""),
            # Tier thresholds (heuristic ops estimates; config-overridable).
            queued_ops_threshold=_int(env, "QLWEB_QUEUED_OPS_THRESHOLD", 500_000_000),
            queued_max_degree=_int(env, "QLWEB_QUEUED_MAX_DEGREE", 20),
            big_ops_threshold=_int(env, "QLWEB_BIG_OPS_THRESHOLD", 50_000_000_000),
            big_max_degree=_int(env, "QLWEB_BIG_MAX_DEGREE", 40),
            big_job_wall_seconds=_int(env, "QLWEB_BIG_JOB_WALL_SECONDS", 4 * 3600),
            big_job_mem_bytes=_int(env, "QLWEB_BIG_JOB_MEM_BYTES", 16 * 1024 ** 3),
            # Outbound SMTP relay for verification/completion mail (spec §17).
            # Unset (default) disables the big-job tier.
            smtp_host=env.get("QLWEB_SMTP_HOST", ""),
            smtp_port=_int(env, "QLWEB_SMTP_PORT", 587),
            smtp_user=env.get("QLWEB_SMTP_USER", ""),
            smtp_pass=env.get("QLWEB_SMTP_PASS", ""),
            smtp_from=env.get("QLWEB_SMTP_FROM", ""),
            token_secret=env.get("QLWEB_TOKEN_SECRET", "quiverlab-dev-token-secret"),
            # Salt for the per-email rate-limit hash. Optional: defaults to
            # token_secret (the original choice) so nothing breaks if it is unset;
            # set it to decouple the two secrets and to rotate rate-limit buckets
            # independently of the magic-link signing key.
            email_hash_salt=(env.get("QLWEB_EMAIL_HASH_SALT")
                             or env.get("QLWEB_TOKEN_SECRET", "quiverlab-dev-token-secret")),
            big_token_ttl_seconds=_int(env, "QLWEB_BIG_TOKEN_TTL_SECONDS", 3600),
            public_base_url=env.get("QLWEB_PUBLIC_BASE_URL", "http://127.0.0.1:8000"),
            # Empty by default: the Docs nav link is absent (no dead href) unless
            # a docs site URL is configured (spec §3 docs-link).
            docs_url=env.get("QLWEB_DOCS_URL", ""),
            # Result cache (Plan 25). On by default: a finished computation is
            # replayed for any later identical request (across users and tiers),
            # never recomputed. The LRU size cap bounds the pinned artifacts.
            cache_enabled=_int(env, "QLWEB_CACHE_ENABLED", 1) != 0,
            cache_max_entries=_int(env, "QLWEB_CACHE_MAX_ENTRIES", 1000),
        )


@lru_cache(maxsize=1)
def get_config() -> Config:
    return Config.from_env()
