from pathlib import Path

import pytest

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


def test_instant_rate_limit_defaults(tmp_path):
    # Correction #2: the instant-tier flood-throttle knobs exist with sane defaults.
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    assert cfg.instant_rate_max == 60
    assert cfg.instant_rate_window_seconds == 60
    over = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                            "QLWEB_INSTANT_RATE_MAX": "3",
                            "QLWEB_INSTANT_RATE_WINDOW_SECONDS": "10"})
    assert over.instant_rate_max == 3 and over.instant_rate_window_seconds == 10


# --------------------------------------------------------------------------- #
# Correction #4 -- the production-secret guard runs on EVERY create_app boot path
# (not only the env-config boot), with an explicit `enforce_secrets=False` opt-out.
# --------------------------------------------------------------------------- #

def _smtp_default_secret_cfg(tmp_path):
    # Big-job tier enabled (SMTP set) while the signing secret is still the public
    # repo default -> the insecure production config the guard must refuse.
    return Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                            "QLWEB_SMTP_HOST": "relay", "QLWEB_SMTP_FROM": "q@e.org"})


def test_explicit_config_production_secret_is_enforced(tmp_path):
    # The defect: a deploy that constructs Config EXPLICITLY (cfg is not None) used to
    # bypass the guard entirely. Now an explicit insecure cfg is refused at boot.
    from webapp.server.app import create_app
    cfg = _smtp_default_secret_cfg(tmp_path)
    assert cfg.big_jobs_enabled is True
    with pytest.raises(RuntimeError, match="insecure configuration"):
        create_app(cfg)


def test_enforce_secrets_false_opts_out(tmp_path):
    # The documented opt-out (tests / dev harnesses): the same insecure cfg boots
    # when enforce_secrets=False is passed explicitly.
    from webapp.server.app import create_app
    app = create_app(_smtp_default_secret_cfg(tmp_path), enforce_secrets=False)
    assert app is not None


def test_real_secret_boots_without_opt_out(tmp_path):
    # A real (non-default) signing secret satisfies the guard on the default path.
    from webapp.server.app import create_app
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path),
                           "QLWEB_SMTP_HOST": "relay", "QLWEB_SMTP_FROM": "q@e.org",
                           "QLWEB_TOKEN_SECRET": "a-real-non-default-secret"})
    assert create_app(cfg) is not None


def test_smtp_off_config_boots_by_default(tmp_path):
    # The offline / SMTP-off case: big-job tier disabled -> the guard is a no-op even
    # with the default secret, so the offline app (which never sets SMTP) is unaffected.
    from webapp.server.app import create_app
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    assert cfg.big_jobs_enabled is False
    assert create_app(cfg) is not None
