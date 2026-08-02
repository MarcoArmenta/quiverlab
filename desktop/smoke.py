"""Cross-platform smoke for the frozen desktop app (CI + local).

    python desktop/smoke.py dist/QuiverLab        # or dist\\QuiverLab.exe

Starts the binary headless on a fixed port with an isolated data dir, then
asserts the two things a user needs: the GUI page answers, and a REAL exact
computation round-trips (quantum CI, q=1, GF(2): HH^0..2 = [4, 8, 12] -- a
frozen-bundle pin, the value itself is oracle-gated in the test suite). Serving
alone is NOT enough: the missing multiprocessing.freeze_support() bug served
pages fine and broke every queued job child.
"""
import json
import os
import subprocess
import sys
import tempfile
import time
import urllib.request

PORT = "8190"
BASE = f"http://127.0.0.1:{PORT}"
REQUEST = {
    "schema": 1,
    "algebra": {"kind": "family", "family": "QuantumCI", "params": {"q": 1},
                 "field": {"kind": "GF", "p": 2, "n": 1}},
    "compute": ["hh_cohomology:0..2"],
    "artifacts": {"pdf": False, "tikz": False},
}
EXPECTED_DIMS = [4, 8, 12]

# The COMPLETE-report check (Marco: "check the binaries for the app on all OS give
# a complete report with all product tables"). SAME tiny algebra, schema 2, but the
# whole HH product surface AND report artifacts (pdf=true, tikz=true). pdf=true forces
# the QUEUED tier -- the instant tier discards its artifact dir, so a report request
# served instantly returns no report at all (see estimator.classify: report_artifacts).
# Every kind computes within caps at tops 0..2 on this 4-dim algebra (verified before
# committing), so the report carries a real product Cayley table, not a vanishing note.
REPORT_REQUEST = {
    "schema": 2,
    "algebra": {"kind": "family", "family": "QuantumCI", "params": {"q": 1},
                 "field": {"kind": "GF", "p": 2, "n": 1}},
    "compute": ["hh_cohomology:0..2", "hh_homology:0..2", "cup:0..2", "cap:0..2",
                "bracket:0..2", "connes_b:0..2", "cyclic_homology:0..2"],
    "artifacts": {"pdf": True, "tikz": True},
}
# STRUCTURAL tokens from the renderers (quiverlab.trace.render_html / results_html),
# not prose that could be reworded: the emitted product-table element, the cup
# results-section anchor, the ordered-basis enumeration list, and the JSON-guide
# appendix anchor. All four are gone if the report is not the complete record.
REPORT_MARKERS = (
    '<table class="ql-matrix ql-cayley">',   # a product Cayley table (cup/cap/bracket)
    "id='cr-cup'",                           # the cup product results-section anchor
    "class='ql-enum'",                       # the ordered-basis enumeration list
    "id='json-guide'",                       # the "Reading the JSON record" appendix
)
# Every product/cyclic kind must be present AND non-refusing in results.
REPORT_PRODUCT_KINDS = ("cup", "cap", "bracket", "connes_b", "cyclic_homology")


def _find_dims(obj):
    """Recursively find the hh_cohomology dims list in a result payload."""
    if isinstance(obj, dict):
        if obj.get("kind") == "HH^" and "dims" in obj:
            return obj["dims"]
        for v in obj.values():
            found = _find_dims(v)
            if found is not None:
                return found
    elif isinstance(obj, list):
        for v in obj:
            found = _find_dims(v)
            if found is not None:
                return found
    return None


def _get_json(url, timeout=10):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return json.loads(r.read())


def _download_text(url, timeout=30):
    """GET an artifact download route; assert 200 and return the decoded body.
    urllib raises HTTPError on a non-2xx, so reaching the read means the artifact
    served -- we assert the status explicitly for a clear failure line."""
    with urllib.request.urlopen(url, timeout=timeout) as r:
        assert r.status == 200, f"{url} -> HTTP {r.status}"
        return r.read().decode("utf-8", "replace")


def _poll_done(jid):
    """Poll /api/jobs/{jid} to a terminal state (generous: the frozen binary runs
    the pure path from a cold start), asserting it finished 'done'."""
    for _ in range(300):                       # 300 x 2s = 10 min ceiling
        time.sleep(2)
        st = _get_json(f"{BASE}/api/jobs/{jid}")
        if st.get("status") in ("done", "failed", "error"):
            return st
    raise SystemExit(f"job {jid} never reached a terminal state")


def main() -> int:
    binary = os.path.abspath(sys.argv[1])
    env = dict(os.environ,
               QUIVERLAB_DESKTOP_NO_BROWSER="1",
               QUIVERLAB_DESKTOP_PORT=PORT,
               QUIVERLAB_DATA=tempfile.mkdtemp(prefix="ql-smoke-data-"))
    log = open(os.path.join(tempfile.gettempdir(), "ql-desktop-smoke.log"),
               "w", encoding="utf-8")
    proc = subprocess.Popen([binary], env=env, stdout=log,
                            stderr=subprocess.STDOUT)
    try:
        for _ in range(360):  # frozen onefile extracts on first start: be patient
            time.sleep(0.5)
            if proc.poll() is not None:
                raise SystemExit(f"binary exited early with {proc.returncode}")
            try:
                with urllib.request.urlopen(BASE + "/", timeout=2) as r:
                    page = r.read().decode("utf-8", "replace")
                    if r.status == 200:
                        break
            except OSError:
                pass
        else:
            raise SystemExit("server never answered")
        assert "quiverlab" in page.lower(), page[:200]
        print("smoke: GUI page OK")

        with urllib.request.urlopen(BASE + "/draw", timeout=10) as r:
            draw = r.read().decode("utf-8", "replace")
        assert r.status == 200 and 'id="qlgui"' in draw, draw[:200]
        with urllib.request.urlopen(BASE + "/gui/worker.js", timeout=10) as r:
            assert r.status == 200 and b"postMessage" in r.read()
        print("smoke: draw-a-quiver canvas OK")

        req = urllib.request.Request(
            BASE + "/api/compute", data=json.dumps(REQUEST).encode(),
            headers={"Content-Type": "application/json"})
        out = _get_json_from(req)
        if out.get("tier") == "queued":
            jid = out["job_id"]
            for _ in range(300):
                time.sleep(2)
                out = _get_json(f"{BASE}/api/jobs/{jid}")
                if out.get("status") in ("done", "failed", "error"):
                    break
            assert out.get("status") == "done", out
        dims = _find_dims(out)
        assert dims == EXPECTED_DIMS, f"dims {dims} != {EXPECTED_DIMS}"
        print(f"smoke: exact compute OK (HH^0..2 = {dims})")

        # ---- The COMPLETE worked-steps report with product Cayley tables. This is
        # what the binary must actually deliver end to end: a queued report job whose
        # downloadable HTML/JSON carry every product table.
        req2 = urllib.request.Request(
            BASE + "/api/compute", data=json.dumps(REPORT_REQUEST).encode(),
            headers={"Content-Type": "application/json"})
        out2 = _get_json_from(req2)
        # (1) pdf=true must NOT serve instant -- the instant tier keeps no artifacts.
        assert out2.get("tier") == "queued", f"pdf=true must queue, got {out2}"
        jid = out2["job_id"]
        st = _poll_done(jid)
        assert st.get("status") == "done", st
        print("smoke: report job queued + computed to done")

        # (2) the downloaded HTML report is the complete record; trace.json parses.
        report = _download_text(f"{BASE}/download/{jid}/trace_steps.html")
        missing = [m for m in REPORT_MARKERS if m not in report]
        assert not missing, f"report missing markers: {missing}"
        json.loads(_download_text(f"{BASE}/download/{jid}/trace.json"))
        print("smoke: report HTML complete (Cayley table + cup anchor + ordered "
              "basis + json-guide) and trace.json parses")

        # (3) result.json exposes the json guide and the cup product tables.
        result = json.loads(_download_text(f"{BASE}/download/{jid}/result.json"))
        assert result.get("json_guide"), "result.json json_guide missing/empty"
        results = result.get("results", {})
        assert results.get("cup", {}).get("tables"), "results.cup.tables missing/empty"
        print("smoke: result.json OK (json_guide + cup tables non-empty)")

        # (4) every product/cyclic kind computed within caps -- none refused.
        for kind in REPORT_PRODUCT_KINDS:
            block = results.get(kind)
            assert isinstance(block, dict), f"result missing {kind}: {block}"
            assert "error" not in block, f"{kind} refused: {block.get('error')}"
        print("smoke: all product tables present (%s)"
              % ", ".join(REPORT_PRODUCT_KINDS))
        return 0
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
        log.close()


def _get_json_from(req):
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read())


if __name__ == "__main__":
    sys.exit(main())
