#!/usr/bin/env python
"""Full-stack local smoke for quiverlab-web (Task 15, acceptance).

Starts the REAL stack as OS subprocesses -- uvicorn running the app factory
(``webapp.server.app:create_app``) and ONE worker poll loop
(``python -m webapp.worker.run_loop``) -- in a throwaway data dir, then drives it
over HTTP with ``httpx`` exactly as a browser/client would:

  1. GET /            -- English landing page renders.
  2. GET /es/         -- Spanish landing page renders (localized string present).
  3. GET /api/catalog -- the introspected family catalog.
  4. POST /api/compute (instant GF(p))          -- runs synchronously, dims inline.
  5. POST /api/compute (queued tier)             -- worker runs it; poll
     /api/jobs/{id} to `done`, then download result.json + the worked-steps
     artifact.
  6. POST /api/feedback                          -- round-trips to a reference id.
  7. POST /api/jobs/big (SMTP off)               -- clean 503 "disabled" response.

Then both subprocesses are shut down gracefully and the script VERIFIES no
process is left alive in either subprocess group (no orphans).

This is intentionally a script, not a pytest: it allocates a real TCP port and
manages OS process teardown, which would make a CI unit test flaky. The
equivalent flow runs in-process (robust, no ports) in
``tests/webapp/test_acceptance.py``.

Run:  .venv/bin/python scripts/webapp_smoke.py
Exit code 0 == every step asserted OK and both process groups are clean.
"""
from __future__ import annotations

import json
import os
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import httpx

REPO_ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable

# QuantumCI(q=1): dim 4, the real reference family (the brief's
# `truncated_polynomial` is not in `quiverlab.families()`).
_QCI_GF2 = {"kind": "family", "family": "QuantumCI", "params": {"q": 1},
            "field": {"kind": "GF", "p": 2, "n": 1}}
_QCI_CC = {"kind": "family", "family": "QuantumCI", "params": {"q": 1},
           "field": {"kind": "CC"}}

_PASS = "\033[32mPASS\033[0m" if sys.stdout.isatty() else "PASS"
_FAIL = "\033[31mFAIL\033[0m" if sys.stdout.isatty() else "FAIL"

_step_n = 0


def check(label: str, ok: bool, detail: str = "") -> None:
    global _step_n
    _step_n += 1
    tag = _PASS if ok else _FAIL
    line = f"  [{_step_n:02d}] {tag}  {label}"
    if detail:
        line += f"  --  {detail}"
    print(line, flush=True)
    if not ok:
        raise AssertionError(label + (f" ({detail})" if detail else ""))


def free_port() -> int:
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


def wait_ready(base: str, timeout: float = 40.0) -> None:
    """Poll GET /api/catalog until the server answers 200 (or time out)."""
    deadline = time.monotonic() + timeout
    last = ""
    while time.monotonic() < deadline:
        try:
            r = httpx.get(base + "/api/catalog", timeout=2.0)
            if r.status_code == 200:
                return
            last = f"status {r.status_code}"
        except Exception as exc:                       # not up yet
            last = type(exc).__name__
        time.sleep(0.3)
    raise TimeoutError(f"server did not become ready in {timeout}s ({last})")


def group_alive(pgid: int) -> bool:
    """True iff any process remains in the group (signal 0 is an existence probe)."""
    try:
        os.killpg(pgid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:                            # exists but not ours to signal
        return True


def stop_group(proc: subprocess.Popen, name: str, grace: float = 25.0) -> None:
    """Graceful stop: SIGTERM the group leader (run_loop drains in-flight work and
    reaps its own loop children; uvicorn shuts down its server), wait, then SIGKILL
    the whole group if anything is still standing. Finally VERIFY the group empty."""
    pgid = proc.pid                                    # start_new_session => leader == pgid
    if proc.poll() is None:
        try:
            os.killpg(pgid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    try:
        proc.wait(timeout=grace)
    except subprocess.TimeoutExpired:
        pass
    # Orphan guard: nothing must survive in the group.
    deadline = time.monotonic() + 10.0
    while group_alive(pgid) and time.monotonic() < deadline:
        try:
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            break
        time.sleep(0.3)
    orphaned = group_alive(pgid)
    check(f"{name}: process group fully torn down (no orphans)", not orphaned,
          f"pgid={pgid} rc={proc.returncode}")


def main() -> int:
    data_dir = Path(tempfile.mkdtemp(prefix="qlweb-smoke-"))
    port = free_port()
    base = f"http://127.0.0.1:{port}"
    env = dict(os.environ)
    env.update({
        "QLWEB_DATA_DIR": str(data_dir),
        # Force a degree-4 request into the queued tier while a degree<=2 request
        # stays instant -- lets ONE server exercise both tiers.
        "QLWEB_INSTANT_MAX_DEGREE": "2",
        "QLWEB_WORKER_PROCESSES": "1",                 # exactly ONE worker loop
        "QLWEB_IP_HASH_SALT": "smoke-salt",
        # SMTP deliberately unset -> the big-job tier is disabled (503 leg).
        "NUMBA_NUM_THREADS": "2", "OMP_NUM_THREADS": "2",
        "PYTHONUNBUFFERED": "1",
    })

    print(f"quiverlab-web full-stack smoke", flush=True)
    print(f"  data_dir = {data_dir}", flush=True)
    print(f"  base_url = {base}\n", flush=True)

    server = worker = None
    server_log = open(data_dir / "uvicorn.log", "w")
    worker_log = open(data_dir / "worker.log", "w")
    try:
        server = subprocess.Popen(
            [PY, "-m", "uvicorn", "webapp.server.app:create_app", "--factory",
             "--host", "127.0.0.1", "--port", str(port), "--log-level", "warning"],
            cwd=REPO_ROOT, env=env, stdout=server_log, stderr=subprocess.STDOUT,
            start_new_session=True)
        worker = subprocess.Popen(
            [PY, "-m", "webapp.worker.run_loop"],
            cwd=REPO_ROOT, env=env, stdout=worker_log, stderr=subprocess.STDOUT,
            start_new_session=True)

        wait_ready(base)
        check("uvicorn + worker started; server answers /api/catalog", True,
              f"server pid={server.pid} worker pid={worker.pid}")

        with httpx.Client(base_url=base, timeout=30.0,
                          follow_redirects=True) as c:
            # 1. English landing page.
            r = c.get("/")
            check("GET / (EN) renders", r.status_code == 200 and
                  "Compute with quivers, exactly" in r.text)

            # 2. Spanish landing page (/es/ 307-redirects to /es).
            r = c.get("/es/")
            check("GET /es/ (ES) renders localized string",
                  r.status_code == 200 and "Polinomio de Coxeter" in r.text)

            # 3. Catalog.
            r = c.get("/api/catalog")
            fams = [f["name"] for f in r.json()["families"]]
            check("GET /api/catalog lists families",
                  r.status_code == 200 and "QuantumCI" in fams,
                  f"{len(fams)} families")

            # 4. Instant GF(p) compute (degree <= 2 stays instant): dims inline.
            r = c.post("/api/compute", json={
                "schema": 1, "algebra": _QCI_GF2,
                "compute": ["hh_cohomology:0..2"],
                "artifacts": {"pdf": False, "tikz": False}})
            j = r.json()
            dims = j.get("result", {}).get("results", {}).get(
                "hh_cohomology", {}).get("dims")
            check("POST /api/compute (instant) -> tier=instant, dims present",
                  r.status_code == 200 and j["tier"] == "instant" and dims == [4, 8, 12],
                  f"dims={dims}")

            # 5. Queued compute (degree 4 > instant cap) -> worker runs it.
            r = c.post("/api/compute", json={
                "schema": 1, "algebra": _QCI_GF2,
                "compute": ["hh_cohomology:0..4"],
                "artifacts": {"pdf": True, "tikz": False}})
            j = r.json()
            check("POST /api/compute (deep) -> 202 tier=queued",
                  r.status_code == 202 and j["tier"] == "queued", str(j))
            jid = j["job_id"]

            # Poll until the worker finishes the job.
            deadline = time.monotonic() + 60.0
            status = None
            while time.monotonic() < deadline:
                s = c.get(f"/api/jobs/{jid}").json()
                status = s["status"]
                if status in ("done", "failed"):
                    break
                time.sleep(0.5)
            check(f"poll /api/jobs/{{id}} -> done", status == "done",
                  f"status={status}")

            # Download the result artifact + assert dims.
            dl = c.get(f"/download/{jid}/result.json")
            res = json.loads(dl.content)
            qdims = res["results"]["hh_cohomology"]["dims"]
            check("GET /download/{id}/result.json -> 200, dims list",
                  dl.status_code == 200 and qdims == [4, 8, 12, 16, 20],
                  f"dims={qdims}")

            # Download the worked-steps artifact (PDF if a LaTeX toolchain is on
            # PATH, else the self-contained HTML fallback).
            art = c.get(f"/download/{jid}/trace.pdf")
            name = "trace.pdf"
            if art.status_code != 200:
                art = c.get(f"/download/{jid}/trace_steps.html")
                name = "trace_steps.html"
            check(f"GET /download/{{id}}/{name} (worked steps) -> 200",
                  art.status_code == 200, f"{len(art.content)} bytes")

            # 6. Feedback round-trip.
            r = c.post("/api/feedback", json={
                "category": "feature",
                "message": "Please add support for gentle algebras.",
                "contact": "", "job_ref": "", "website": ""})
            ref = r.json().get("reference")
            check("POST /api/feedback -> 201 with reference",
                  r.status_code == 201 and bool(ref), f"reference={ref}")

            # 7. Big-job tier with SMTP off -> clean disabled response.
            r = c.post("/api/jobs/big", json={
                "schema": 1, "algebra": _QCI_CC,
                "compute": ["hh_cohomology:0..30"],
                "artifacts": {"pdf": False, "tikz": False},
                "email": "user@example.org", "lang": "en"})
            check("POST /api/jobs/big (SMTP off) -> 503 BigJobsDisabled",
                  r.status_code == 503 and
                  r.json().get("error_type") == "BigJobsDisabled", str(r.json()))

        print("\n  All HTTP steps passed; shutting the stack down...\n", flush=True)
    finally:
        # Clean, verified teardown -- no orphaned processes.
        if worker is not None:
            stop_group(worker, "worker")
        if server is not None:
            stop_group(server, "uvicorn")
        server_log.close()
        worker_log.close()

    # Sanity: no ./quiverlab_traces at the repo root (worker child chdir's into
    # the job artifact dir; a leak here would be a regression).
    check("no quiverlab_traces/ leaked at repo root",
          not (REPO_ROOT / "quiverlab_traces").exists())

    shutil.rmtree(data_dir, ignore_errors=True)
    print("\nSMOKE OK -- full stack accepted.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
