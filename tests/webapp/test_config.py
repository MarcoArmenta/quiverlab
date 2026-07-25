from pathlib import Path

from webapp.server.config import Config


def test_from_env_defaults(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    assert cfg.data_dir == tmp_path
    assert cfg.db_path == tmp_path / "quiverlab_web.sqlite3"
    assert cfg.artifacts_dir == tmp_path / "artifacts"
    assert cfg.instant_wall_seconds == 5
    assert cfg.instant_ops_threshold == 2_000_000
    assert cfg.job_wall_seconds == 900
    assert cfg.job_mem_bytes == 4 * 1024 ** 3
    assert cfg.result_max_bytes == 32 * 1024 ** 2
    assert cfg.retention_days == 90
    assert cfg.ip_hash_salt  # non-empty default so dev/tests are deterministic
    assert cfg.feedback_daily_max == 5
    assert cfg.admin_token == ""      # admin route disabled unless a token is set
    assert cfg.big_job_wall_seconds == 4 * 3600
    assert cfg.big_job_mem_bytes == 16 * 1024 ** 3
    assert cfg.per_email_weekly_max == 5
    assert cfg.big_jobs_enabled is False   # no SMTP configured by default
    assert cfg.docs_url == ""              # Docs nav link absent by default


def test_env_overrides(tmp_path):
    cfg = Config.from_env({
        "QLWEB_DATA_DIR": str(tmp_path),
        "QLWEB_JOB_WALL_SECONDS": "120",
        "QLWEB_INSTANT_OPS_THRESHOLD": "500000",
        "QLWEB_RESULT_MAX_BYTES": "1024",
        "QLWEB_IP_HASH_SALT": "prod-secret",
        "QLWEB_FEEDBACK_DAILY_MAX": "9",
        "QLWEB_ADMIN_TOKEN": "s3cret",
        "QLWEB_BIG_JOB_WALL_SECONDS": "7200",
    })
    assert cfg.job_wall_seconds == 120
    assert cfg.instant_ops_threshold == 500_000
    assert cfg.result_max_bytes == 1024
    assert cfg.ip_hash_salt == "prod-secret"
    assert cfg.feedback_daily_max == 9
    assert cfg.admin_token == "s3cret"
    assert cfg.big_job_wall_seconds == 7200


def test_smtp_enables_big_jobs(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_SMTP_HOST": "smtp.relay.example",
                           "QLWEB_SMTP_FROM": "quiverlab@example.org"})
    assert cfg.big_jobs_enabled is True
