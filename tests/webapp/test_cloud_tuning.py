"""Cloud capacity tuning (webapp/deploy) contract tests.

The deployed website ships a tuned cloud profile in ``docker-compose.yml`` for
the DRAC persistent instance (~16 vCPU / 50 GB RAM / 200 GB /data). These tests
prove that profile:

- parses cleanly through ``Config.from_env`` and yields the intended values,
- stays within the RAM arithmetic (per-job RLIMIT_AS x worst-case concurrent
  workers must fit RAM with margin -- the shared FIFO queue means the worst case
  is every worker running a big job),
- leaves the INSTANT tier at the library code defaults (the anonymous DoS
  surface is off-limits), and
- keeps the graceful-drain / requeue relationship coherent with the longer walls.

Nothing here touches ``src/`` config defaults: the profile lives entirely in the
compose ``environment:`` values, so ``tests/webapp/test_config.py`` (which pins
the code defaults) stays green.
"""
import re
from pathlib import Path

import yaml

from webapp.server.config import Config

DEPLOY = Path("webapp/deploy")
GiB = 1024 ** 3

# ---- Minimal compose ${VAR}-interpolation, matching docker compose semantics. --
# ${VAR}          -> override or ""              (unused here)
# ${VAR:-default} -> override or default
# ${VAR:?message} -> override or a dummy (a required var; dummy lets parse proceed)
_VAR = re.compile(r"\$\{([A-Z_][A-Z0-9_]*)(?::-(?P<d>[^}]*)|:\?(?P<m>[^}]*))?\}")


def _interp(value: str, overrides: dict[str, str]) -> str:
    def sub(m: "re.Match") -> str:
        name, dflt = m.group(1), m.group("d")
        if name in overrides:
            return overrides[name]
        if dflt is not None:
            return dflt
        return f"dummy-{name.lower()}"          # required ${VAR:?...}, no value given
    return _VAR.sub(sub, str(value))


def _load_compose() -> dict:
    return yaml.safe_load((DEPLOY / "docker-compose.yml").read_text(encoding="utf-8"))


def _config_from_compose(overrides: dict[str, str] | None = None) -> Config:
    """Resolve the app service's environment the way `docker compose` would, then
    feed it to the real config loader -- the mandate's 'compose env values parse
    through Config' check."""
    overrides = overrides or {}
    env = _load_compose()["services"]["app"]["environment"]
    resolved = {k: _interp(v, overrides) for k, v in env.items()}
    return Config.from_env(resolved)


# --------------------------------------------------------------------------- #
# 1. The compose profile parses through Config and yields the tuned values.
# --------------------------------------------------------------------------- #

def test_compose_tuning_parses_through_config():
    cfg = _config_from_compose()
    # Queued (anonymous large) tier raised from a laptop-grade cap.
    assert cfg.job_wall_seconds == 3600            # 15 min -> 1 h
    assert cfg.job_mem_bytes == 8 * GiB            # 4 -> 8 GiB
    # Big (email-verified HPC) tier.
    assert cfg.big_job_wall_seconds == 86400       # 4 -> 24 h
    assert cfg.big_job_mem_bytes == 20 * GiB       # 16 -> 20 GiB
    # Fleet size (RAM-bound).
    assert cfg.worker_processes == 2
    # Tier-routing thresholds accept what a laptop can't.
    assert cfg.queued_ops_threshold == 20_000_000_000
    assert cfg.queued_max_degree == 30
    assert cfg.big_ops_threshold == 500_000_000_000
    assert cfg.big_max_degree == 60
    # Cache + retention sized for the 200 GB volume.
    assert cfg.cache_max_entries == 10_000
    assert cfg.retention_days == 365
    # Backlog + fairness budgets.
    assert cfg.global_queue_max == 1000
    assert cfg.big_queue_max == 50
    assert cfg.per_email_weekly_max == 10
    assert cfg.per_ip_daily_max == 300


def test_app_and_worker_share_one_env_anchor():
    # A single YAML anchor feeds both services, so the app's estimator and the
    # worker's per-job caps can never drift apart.
    compose = _load_compose()
    app_env = compose["services"]["app"]["environment"]
    wrk_env = compose["services"]["worker"]["environment"]
    assert app_env == wrk_env
    # Both still carry the required-secret guards (existing deploy test relies on
    # this too, but assert here since the anchor restructured the file).
    for env in (app_env, wrk_env):
        assert "${QLWEB_IP_HASH_SALT:?" in env["QLWEB_IP_HASH_SALT"]
        assert "${QLWEB_TOKEN_SECRET:?" in env["QLWEB_TOKEN_SECRET"]


# --------------------------------------------------------------------------- #
# 2. Every knob is overridable from .env (the ${VAR:-default} form).
# --------------------------------------------------------------------------- #

def test_compose_tuning_is_overridable():
    cfg = _config_from_compose({
        "QLWEB_JOB_WALL_SECONDS": "7200",
        "QLWEB_BIG_JOB_MEM_BYTES": str(40 * GiB),
        "QLWEB_WORKER_PROCESSES": "1",
    })
    assert cfg.job_wall_seconds == 7200
    assert cfg.big_job_mem_bytes == 40 * GiB
    assert cfg.worker_processes == 1
    # The documented single-huge-job campaign profile (1 worker x 40 GiB) is
    # itself RAM-safe on a 50 GiB box.
    assert cfg.worker_processes * cfg.big_job_mem_bytes <= 50 * GiB


# --------------------------------------------------------------------------- #
# 3. RAM arithmetic: worst-case (all workers big) fits with margin.
# --------------------------------------------------------------------------- #

def test_ram_arithmetic_fits_with_margin():
    cfg = _config_from_compose()
    total_ram = 50 * GiB                                   # the persistent-instance box
    # The queue is a single shared FIFO with no per-tier priority, so the worst
    # case is EVERY worker running a big job at once.
    fleet_worst = cfg.worker_processes * cfg.big_job_mem_bytes
    assert fleet_worst <= 40 * GiB                          # documented usable budget
    headroom = total_ram - fleet_worst
    assert headroom >= 8 * GiB      # OS + app + caddy + SQLite + instant children
    # Queued-only worst case is smaller still.
    assert cfg.worker_processes * cfg.job_mem_bytes <= fleet_worst


# --------------------------------------------------------------------------- #
# 4. The INSTANT tier is UNCHANGED -- compose must NOT override its knobs, so
#    Config falls back to the library code defaults (the DoS surface stays put).
# --------------------------------------------------------------------------- #

def test_instant_tier_left_at_code_defaults():
    env = _load_compose()["services"]["app"]["environment"]
    for knob in ("QLWEB_INSTANT_WALL_SECONDS", "QLWEB_INSTANT_OPS_THRESHOLD",
                 "QLWEB_INSTANT_MAX_DEGREE", "QLWEB_INSTANT_RATE_MAX",
                 "QLWEB_INSTANT_RATE_WINDOW_SECONDS"):
        assert knob not in env, f"cloud profile must not tune the instant knob {knob}"
    cfg = _config_from_compose()
    assert cfg.instant_wall_seconds == 5
    assert cfg.instant_ops_threshold == 2_000_000
    assert cfg.instant_max_degree == 8
    assert cfg.instant_rate_max == 60
    assert cfg.instant_rate_window_seconds == 60
    # Per-identity RUNNING caps stay 1 (the concurrency + RAM bound).
    assert cfg.per_ip_running_max == 1
    assert cfg.per_email_running_max == 1


# --------------------------------------------------------------------------- #
# 5. Graceful-stop / requeue coherence with the longer queued wall.
# --------------------------------------------------------------------------- #

def test_stop_grace_covers_queued_wall_plus_join_slack():
    # run_loop.main() joins its loops within (job_wall + 20s); compose's
    # stop_grace_period must sit at/above that so main reaps cleanly before
    # Docker's SIGKILL. Big jobs (24 h) are intentionally not covered -- they are
    # requeued at next startup.
    compose = _load_compose()
    cfg = _config_from_compose()
    grace = int(str(compose["services"]["worker"]["stop_grace_period"]).rstrip("s"))
    assert grace >= cfg.job_wall_seconds + 20
    # But the grace must NOT balloon to the big wall (a deploy can't block a day).
    assert grace < cfg.big_job_wall_seconds


# --------------------------------------------------------------------------- #
# 6. The profile is documented (arithmetic, burst story, honest won't-fit note).
# --------------------------------------------------------------------------- #

def test_provisioning_documents_cloud_tuning():
    text = (DEPLOY / "PROVISIONING.md").read_text(encoding="utf-8")
    low = text.lower()
    assert "cloud capacity tuning" in low
    # The RAM arithmetic is stated.
    assert "QLWEB_WORKER_PROCESSES" in text and "QLWEB_BIG_JOB_MEM_BYTES" in text
    # The burst-instance story points at the HPC batch CLI.
    assert "80 vcpu" in low and "quiverlab-hpc" in low
    # Honest scope: the dim-220 product examples still won't fit.
    assert "220" in text


def test_env_example_documents_overrides():
    text = (DEPLOY / ".env.example").read_text(encoding="utf-8")
    assert "CLOUD CAPACITY TUNING" in text
    # A few representative override knobs are shown.
    for knob in ("QLWEB_WORKER_PROCESSES", "QLWEB_BIG_JOB_MEM_BYTES",
                 "QLWEB_JOB_WALL_SECONDS", "QLWEB_RETENTION_DAYS"):
        assert knob in text
    # The RAM rule is spelled out for the operator.
    assert "QLWEB_WORKER_PROCESSES * QLWEB_BIG_JOB_MEM_BYTES" in text
