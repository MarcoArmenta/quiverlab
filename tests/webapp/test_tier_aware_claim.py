"""Fix M1 -- tier-aware claiming so a run of big jobs can never freeze the
anonymous instant/queued tier.

``JobStore.claim_next`` gains an ``exclude_tiers`` filter (a pure narrowing of the
eligible rows: same BEGIN IMMEDIATE atomicity, same FIFO order). The worker fleet
(``webapp/worker/run_loop.py``) reserves loop 0 for the instant/queued tier by
excluding ``('big',)`` -- but only when the fleet has >= 2 loops, so a single-loop
fleet (the offline embedded worker) stays byte-identical and still serves big jobs.

These tests drive the store + the pure fleet-decision helper directly; they do not
spawn worker processes (that path is exercised by test_deploy_assets.py)."""
import pytest

from webapp.server.config import Config
from webapp.server.store import JobStore
from webapp.worker.run_loop import exclude_tiers_for

pytestmark = pytest.mark.fast


def _store(tmp_path):
    cfg = Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})
    s = JobStore(cfg.db_path)
    s.init_schema()
    return s


# --------------------------------------------------------------------------- #
# claim_next(exclude_tiers=...) narrows the eligible rows, keeps FIFO + atomicity.
# --------------------------------------------------------------------------- #

def test_claim_next_default_claims_any_tier_fifo(tmp_path):
    s = _store(tmp_path)
    first = s.create_job({"schema": 1}, ip="h", tier="big")
    s.create_job({"schema": 1}, ip="h", tier="anonymous")
    # No exclusion -> oldest row wins, byte-identical to before.
    assert s.claim_next().id == first


def test_claim_next_excluding_big_skips_big_rows(tmp_path):
    s = _store(tmp_path)
    # Big is oldest; the anonymous (queued) job is newest.
    s.create_job({"schema": 1}, ip="h", tier="big")
    queued = s.create_job({"schema": 1}, ip="h", tier="anonymous")
    # Excluding big claims the queued job even though it is NOT the oldest overall.
    claimed = s.claim_next(exclude_tiers=("big",))
    assert claimed is not None and claimed.id == queued
    assert claimed.status == "running"


def test_claim_next_excluding_big_returns_none_when_only_big_pending(tmp_path):
    s = _store(tmp_path)
    s.create_job({"schema": 1}, ip="h", tier="big")
    # Loop 0 finds nothing eligible and idles -- it does NOT claim the big job.
    assert s.claim_next(exclude_tiers=("big",)) is None
    # ...but an unrestricted loop still claims it.
    assert s.claim_next() is not None


def test_claim_next_exclude_preserves_fifo_within_eligible(tmp_path):
    s = _store(tmp_path)
    a1 = s.create_job({"schema": 1}, ip="h", tier="anonymous")
    s.create_job({"schema": 1}, ip="h", tier="big")
    s.create_job({"schema": 1}, ip="h", tier="anonymous")   # newer anonymous
    # Among the eligible (non-big) rows the OLDEST is still chosen.
    assert s.claim_next(exclude_tiers=("big",)).id == a1


# --------------------------------------------------------------------------- #
# The mandate scenario: [big, big, queued] with a fleet of 2.
# --------------------------------------------------------------------------- #

def test_fleet_two_serves_queued_while_one_big_waits(tmp_path):
    s = _store(tmp_path)
    big1 = s.create_job({"schema": 1}, ip="h", tier="big")
    big2 = s.create_job({"schema": 1}, ip="h", tier="big")
    queued = s.create_job({"schema": 1}, ip="h", tier="anonymous")

    # Loop 0 (reserved: excludes big) and loop 1 (claims anything) each take one.
    loop0 = s.claim_next(exclude_tiers=exclude_tiers_for(0, 2))
    loop1 = s.claim_next(exclude_tiers=exclude_tiers_for(1, 2))

    claimed_ids = {loop0.id, loop1.id}
    # The anonymous queued job is served immediately (by loop 0)...
    assert loop0.id == queued
    # ...loop 1 takes the oldest big...
    assert loop1.id == big1
    # ...and exactly one big is still waiting -- the queue is NOT frozen.
    assert big2 not in claimed_ids
    assert s.get_job(big2).status == "pending"
    assert s.count_pending() == 1


def test_fleet_one_claims_everything_including_big(tmp_path):
    s = _store(tmp_path)
    big1 = s.create_job({"schema": 1}, ip="h", tier="big")
    s.create_job({"schema": 1}, ip="h", tier="big")
    # A single-loop fleet must still serve big jobs (the offline app depends on it).
    assert exclude_tiers_for(0, 1) is None
    assert s.claim_next(exclude_tiers=exclude_tiers_for(0, 1)).id == big1


# --------------------------------------------------------------------------- #
# The pure fleet-reservation decision.
# --------------------------------------------------------------------------- #

def test_exclude_tiers_for_reservation_rule():
    # Loop 0 reserves itself for instant/queued ONLY when a second loop exists.
    assert exclude_tiers_for(0, 2) == ("big",)
    assert exclude_tiers_for(0, 3) == ("big",)
    # A single-loop fleet excludes nothing (byte-compatible; offline app).
    assert exclude_tiers_for(0, 1) is None
    # Every non-zero loop claims any tier.
    assert exclude_tiers_for(1, 2) is None
    assert exclude_tiers_for(2, 3) is None
