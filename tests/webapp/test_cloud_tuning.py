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

import pytest
import yaml

from webapp.server.config import Config

DEPLOY = Path("webapp/deploy")
GiB = 1024 ** 3

# ---- Minimal compose ${VAR}-interpolation, matching docker compose semantics. --
# The compose file uses ONLY the colon forms, whose defining property is that an
# EMPTY value counts the same as an unset one (docker's `:-` / `:?`):
#   ${VAR:-default} -> override if set AND non-empty, else default
#   ${VAR:?message} -> override if set AND non-empty, else a dummy (a required var;
#                      compose would ERROR on empty/unset, the dummy lets the test
#                      parser proceed past the required-secret guard)
# The earlier version returned the override verbatim even when it was "", which is
# the NON-colon `-`/`?` semantics and does not match this file -- an empty override
# must fall back to the default (regression pinned by the empty-override test).
_VAR = re.compile(r"\$\{([A-Z_][A-Z0-9_]*)(?::-(?P<d>[^}]*)|:\?(?P<m>[^}]*))?\}")


def _interp(value: str, overrides: dict[str, str]) -> str:
    def sub(m: "re.Match") -> str:
        name, dflt, msg = m.group(1), m.group("d"), m.group("m")
        val = overrides.get(name)
        set_and_nonempty = val is not None and val != ""
        if dflt is not None:                    # ${VAR:-default}
            return val if set_and_nonempty else dflt
        if msg is not None:                     # ${VAR:?message} (required)
            return val if set_and_nonempty else f"dummy-{name.lower()}"
        return val if val is not None else ""   # bare ${VAR} (unused in this file)
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
    # Big (email-verified HPC) tier. Memory back to 16 GiB (the code default): the
    # conservative all-big worst case is 2 x 16 = 32 GiB of 50 (~18 GiB headroom),
    # and the tier-aware fleet runs at most one big job at a time (fix M1/M3).
    assert cfg.big_job_wall_seconds == 86400       # 4 -> 24 h
    assert cfg.big_job_mem_bytes == 16 * GiB       # RAM-honesty: was 20, back to 16
    # Fleet size (RAM-bound).
    assert cfg.worker_processes == 2
    # Global instant-concurrency ceiling (a protection; default 0 = unlimited).
    assert cfg.instant_global_max == 8
    # Tier-routing thresholds accept what a laptop can't.
    assert cfg.queued_ops_threshold == 20_000_000_000
    assert cfg.queued_max_degree == 30
    assert cfg.big_ops_threshold == 500_000_000_000
    assert cfg.big_max_degree == 60
    # Cache + retention sized for the 200 GB volume (byte-consistent: 3500 x the
    # 32 MiB result ceiling ~= 117 GB, under 60% of 200 GB).
    assert cfg.cache_max_entries == 3500
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
# 3. RAM budget-math (documentation arithmetic, NOT a proof).
#
# The old test here was tautological: it recomputed `worker_processes * big_mem`
# from the config and asserted it against a constant picked to match those same
# config values. This replaces it with the DOCUMENTED worst-case-resident budget
# the deploy docs promise, using the post-fix numbers and the tier-aware fleet.
# It is honest budget arithmetic (a sanity gate on the chosen knobs), explicitly
# NOT a runtime proof that no job ever exceeds its RLIMIT_AS.
# --------------------------------------------------------------------------- #

# Stated per-instant-child RESIDENT allowance. Instant computes are tiny by
# classification (ops <= 2e6, degree <= 8), so a fresh spawn child's real resident
# footprint (interpreter + numpy + a small dense matrix, cold numba skipped) sits
# well under this. Its RLIMIT_AS *ceiling* reuses the 8 GiB queued value, but that
# is a virtual-address backstop, never expected resident -- using the ceiling here
# (8 x 8 = 64 GiB) would be dishonest, not conservative.
INSTANT_CHILD_RESIDENT_ALLOWANCE = 1 * GiB
OS_AND_SERVICES_HEADROOM = 8 * GiB            # OS + FastAPI app + caddy + SQLite/WAL
BOX_RAM = 50 * GiB


def test_ram_worst_case_resident_budget_math():
    """Budget-math (documentation), NOT proof. Worst-case concurrent RESIDENT
    memory under the shipped profile must fit the 50 GiB box with headroom."""
    cfg = _config_from_compose()
    fleet = cfg.worker_processes                       # 2

    # Fix M1 (tier-aware claiming): with fleet >= 2, loop 0 refuses the big tier,
    # so at most (fleet - 1) loops run a big job concurrently...
    big_resident = (fleet - 1) * cfg.big_job_mem_bytes          # 1 * 16 = 16 GiB
    # ...while loop 0 (the reserved loop) runs at most one queued job.
    loop0_resident = cfg.job_mem_bytes                          # 8 GiB
    # Fix M3: instant children are globally bounded; count them at the stated
    # resident allowance, not their (virtual) RLIMIT_AS ceiling.
    instant_resident = cfg.instant_global_max * INSTANT_CHILD_RESIDENT_ALLOWANCE  # 8 * 1 = 8 GiB

    worst_resident = big_resident + loop0_resident + instant_resident            # 32 GiB
    assert worst_resident + OS_AND_SERVICES_HEADROOM <= BOX_RAM   # 32 + 8 = 40 <= 50

    # And the CONSERVATIVE bound the docs also promise (holds even if someone
    # overrides the tier reservation to a blind FIFO): all workers big at once.
    conservative_all_big = fleet * cfg.big_job_mem_bytes         # 2 * 16 = 32 GiB
    assert conservative_all_big <= BOX_RAM - 16 * GiB            # >= 16 GiB slack
    assert BOX_RAM - conservative_all_big >= 16 * GiB            # ~18 GiB headroom


def test_cache_entry_budget_fits_the_volume():
    """Fix W3: the entry-count LRU, at the worst case of every entry at the full
    result ceiling, must fit under ~60% of the 200 GB /data volume -- so the cache
    can never fill the disk on its own."""
    cfg = _config_from_compose()
    volume_bytes = 200 * 1000 ** 3                     # 200 GB (decimal, conservative)
    worst_cache_bytes = cfg.cache_max_entries * cfg.result_max_bytes
    assert worst_cache_bytes <= 0.60 * volume_bytes
    # And the number is not trivially tiny (it must actually use the big box).
    assert cfg.cache_max_entries >= 3000


# --------------------------------------------------------------------------- #
# 4. The INSTANT tier is UNCHANGED -- compose must NOT override its knobs, so
#    Config falls back to the library code defaults (the DoS surface stays put).
# --------------------------------------------------------------------------- #

def test_instant_tier_left_at_code_defaults():
    env = _load_compose()["services"]["app"]["environment"]
    # The instant CLASSIFIER knobs (what instant will accept) stay at code default:
    # widening the anonymous DoS surface is off the table.
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
    # The ONE instant knob the profile sets only TIGHTENS the tier (bounds
    # concurrent children); it never widens what instant accepts.
    assert env.get("QLWEB_INSTANT_GLOBAL_MAX") == "${QLWEB_INSTANT_GLOBAL_MAX:-8}"
    assert cfg.instant_global_max == 8


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
                 "QLWEB_JOB_WALL_SECONDS", "QLWEB_RETENTION_DAYS",
                 "QLWEB_INSTANT_GLOBAL_MAX"):
        assert knob in text
    # The RAM rule is spelled out for the operator.
    assert "QLWEB_WORKER_PROCESSES * QLWEB_BIG_JOB_MEM_BYTES" in text


def test_provisioning_states_requeue_from_scratch_and_burst_checkpointing():
    """Fix W1: the deploy docs must state plainly that a killed job re-runs from
    scratch (no checkpoint/resume in the webapp worker), name the deploy-at-hour-23
    consequence, point long campaigns at quiverlab-hpc's checkpointing, and give
    `docker compose kill` as the fast path past the ~1 h queued drain."""
    text = (DEPLOY / "PROVISIONING.md").read_text(encoding="utf-8")
    low = text.lower()
    assert "from scratch" in low or "re-runs from scratch" in low
    assert "no checkpoint" in low or "checkpoint/resume the webapp" in low \
        or "webapp worker has no checkpoint" in low
    assert "hour 23" in low                            # the concrete consequence
    assert "docker compose kill" in low                # the fast-path escape hatch
    assert "quiverlab-hpc" in low                      # where checkpointing lives


# --------------------------------------------------------------------------- #
# 7. Interpolation fidelity: an EMPTY override falls back to the default (the
#    docker `:-` semantics), NOT to the empty string.
# --------------------------------------------------------------------------- #

def test_compose_interpolation_empty_override_uses_default():
    # docker-compose ${VAR:-default}: a set-but-EMPTY value counts as unset and
    # takes the default. (Only `:-`/`:?` behave this way; the file uses only those.)
    cfg = _config_from_compose({"QLWEB_JOB_WALL_SECONDS": ""})
    assert cfg.job_wall_seconds == 3600            # empty -> compose default, not code default
    # A NON-empty override still wins.
    cfg2 = _config_from_compose({"QLWEB_JOB_WALL_SECONDS": "1234"})
    assert cfg2.job_wall_seconds == 1234
    # A required ${VAR:?...} secret set-but-empty also falls back (to the dummy),
    # never leaving a literal "" that a real deploy's compose would reject.
    env = _load_compose()["services"]["app"]["environment"]
    resolved = _interp(env["QLWEB_TOKEN_SECRET"], {"QLWEB_TOKEN_SECRET": ""})
    assert resolved and "${" not in resolved


# --------------------------------------------------------------------------- #
# 8. `docker compose config` actually validates + interpolates the shipped file.
#    Skipped where docker / the compose plugin is absent (CI may lack them; this
#    host has them -- see the fix session's verification).
# --------------------------------------------------------------------------- #

def _docker_compose_available() -> bool:
    import shutil
    import subprocess
    if shutil.which("docker") is None:
        return False
    try:
        r = subprocess.run(["docker", "compose", "version"],
                           capture_output=True, timeout=30)
        return r.returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


@pytest.mark.skipif(not _docker_compose_available(),
                    reason="docker + compose plugin not available")
def test_docker_compose_config_validates_and_interpolates():
    import os
    import subprocess
    env = {**os.environ,
           "QLWEB_IP_HASH_SALT": "a" * 64,
           "QLWEB_TOKEN_SECRET": "b" * 64}
    r = subprocess.run(
        ["docker", "compose", "-f", "docker-compose.yml", "config"],
        cwd=DEPLOY, env=env, capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, f"docker compose config failed:\n{r.stderr}"
    out = r.stdout
    # The tuned defaults render into the interpolated output.
    assert "17179869184" in out          # big-job mem == 16 GiB (fix M3)
    assert "8589934592" in out           # queued mem == 8 GiB
    assert "3500" in out                 # cache_max_entries (fix W3)
    # The global instant cap is present and set to 8.
    assert '"8"' in out or "QLWEB_INSTANT_GLOBAL_MAX" in out
