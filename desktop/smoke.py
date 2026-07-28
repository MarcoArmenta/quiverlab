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
