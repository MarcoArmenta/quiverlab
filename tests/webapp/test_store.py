from webapp.server.store import JobStore


def _store(tmp_path):
    s = JobStore(tmp_path / "jobs.sqlite3")
    s.init_schema()
    return s


def test_create_and_get(tmp_path):
    s = _store(tmp_path)
    jid = s.create_job({"compute": ["cartan"]}, ip="1.2.3.4")
    job = s.get_job(jid)
    assert job.status == "pending"
    assert job.spec["compute"] == ["cartan"]
    assert job.ip == "1.2.3.4"


def test_claim_transitions_to_running(tmp_path):
    s = _store(tmp_path)
    jid = s.create_job({"compute": ["cartan"]}, ip="1.2.3.4")
    claimed = s.claim_next()
    assert claimed.id == jid
    assert s.get_job(jid).status == "running"
    assert s.claim_next() is None  # nothing left pending


def test_requeue_stale_running_makes_claimable(tmp_path):
    s = _store(tmp_path)
    jid = s.create_job({"compute": ["cartan"]}, ip="1.2.3.4")
    assert s.claim_next().id == jid          # flipped to running
    assert s.get_job(jid).status == "running"
    assert s.claim_next() is None            # worker died: row stranded in running
    # Startup requeue adopts the orphan back into the queue and clears started_at.
    assert s.requeue_stale_running() == [jid]
    job = s.get_job(jid)
    assert job.status == "pending"
    assert job.started_at is None
    assert s.claim_next().id == jid          # claimable again
    # A quiescent fleet (nothing running) requeues nothing.
    s.mark_done(jid, artifact_dir="/x")
    assert s.requeue_stale_running() == []


def test_ip_counts(tmp_path):
    s = _store(tmp_path)
    s.create_job({}, ip="9.9.9.9")
    s.claim_next()
    assert s.count_running_for_ip("9.9.9.9") == 1
    assert s.count_running_for_ip("0.0.0.0") == 0


def test_mark_done_and_progress(tmp_path):
    s = _store(tmp_path)
    jid = s.create_job({}, ip="1.1.1.1")
    s.claim_next()
    s.update_progress(jid, {"degree": 3, "of": 6})
    s.mark_done(jid, artifact_dir="/data/artifacts/" + jid)
    job = s.get_job(jid)
    assert job.status == "done"
    assert job.progress == {"degree": 3, "of": 6}


def test_feedback_roundtrip_and_count(tmp_path):
    s = _store(tmp_path)
    fid = s.create_feedback("problem", "HH^3 looks off for A5", "me@example.org",
                            ip="hashy", job_ref="01AN4Z07BY79KA1307SR9X4MV3")
    rows = s.list_feedback()
    assert rows and rows[0]["id"] == fid
    assert rows[0]["category"] == "problem"
    assert rows[0]["job_ref"] == "01AN4Z07BY79KA1307SR9X4MV3"
    assert rows[0]["extra"] is None            # non-structured category
    assert s.count_feedback_today_for_ip("hashy", "2099-01-01T00:00:00Z") == 0
    # created_at defaults to now; count for its own day is 1
    today = rows[0]["created_at"]
    assert s.count_feedback_today_for_ip("hashy", today) == 1
    assert s.count_feedback_today_for_ip("other", today) == 0


def test_feedback_extra_json_roundtrip(tmp_path):
    s = _store(tmp_path)
    import json
    extra = json.dumps({"reference": "arXiv:1406.2300",
                        "why_relevant": "Chouhy-Solotar resolution used here."})
    s.create_feedback("literature", "Please cite this.", None,
                      ip="hashy", job_ref=None, extra=extra)
    row = s.list_feedback()[0]
    assert json.loads(row["extra"])["reference"] == "arXiv:1406.2300"


def test_pending_big_is_single_use(tmp_path):
    s = _store(tmp_path)
    pid = s.create_pending_big({"compute": ["hh_cohomology:0..30"]},
                               email="a@b.c", email_hash="eh", lang="es")
    got = s.consume_pending_big(pid)
    assert got is not None
    assert got["email"] == "a@b.c" and got["lang"] == "es"
    assert s.consume_pending_big(pid) is None      # second use finds nothing


def test_big_job_caps_land_in_row_and_email_clears(tmp_path):
    s = _store(tmp_path)
    jid = s.create_job({"compute": ["cartan"]}, ip="1.2.3.4", tier="big",
                       email="a@b.c", email_hash="eh",
                       wall_seconds=14400, mem_bytes=17179869184, lang="en")
    job = s.get_job(jid)
    assert job.tier == "big" and job.wall_seconds == 14400 and job.mem_bytes == 17179869184
    assert job.email == "a@b.c"
    s.clear_email(jid)
    cleared = s.get_job(jid)
    assert cleared.email is None                   # plaintext gone...
    assert cleared.email_hash == "eh"              # ...hash kept for rate-limiting


def test_big_email_rate_counts(tmp_path):
    s = _store(tmp_path)
    s.create_job({}, ip="i", tier="big", email="a@b.c", email_hash="eh")
    assert s.count_big_running_for_email_hash("eh") == 1
    assert s.count_big_running_for_email_hash("other") == 0
    s.create_job({}, ip="i", tier="big", email_hash="eh")   # second big job, same hash
    assert s.count_big_since_for_email_hash("eh", "1970-01-01T00:00:00Z") == 2
    assert s.count_big_pending() == 2
