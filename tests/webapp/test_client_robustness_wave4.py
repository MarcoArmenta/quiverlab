"""Wave 4 (v1.0.1, 2026-08-21): the compute-page client must survive an error
body of ANY shape.

app.js's compute submit handler printed ``data.error_type + ": " + data.message``
for every non-ok response. But ``/api/compute`` takes a Pydantic model and the
server registers no ``RequestValidationError`` handler, so a schema-invalid body
(an incomplete algebra spec — reachable from the real UI) comes back as FastAPI's
default ``{"detail": [...]}`` 422, where both fields are absent → the box read
"undefined: undefined". The draw-page twin ``webapp/static/gui/worker.js`` already
handles the ``detail`` shape (``protocolError``); these gates pin that the same
handling is ported into ``app.js`` and that a network / non-JSON failure no longer
turns a submit button into a silent no-op.
"""
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from webapp.server.app import create_app
from webapp.server.config import Config

ROOT = Path(__file__).resolve().parents[2]
APP_JS = ROOT / "webapp" / "static" / "app.js"
NODE = shutil.which("node")


def _src():
    return APP_JS.read_text(encoding="utf-8")


def _client(tmp_path):
    return TestClient(create_app(Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})))


def _extract_fn(src, name):
    """Return the source text of a top-level ``function <name>(...) {...}`` by
    brace matching (safe here: none of the target functions embed ``{``/``}`` in
    a string or template literal)."""
    i = src.find("function " + name + "(")
    if i < 0:
        return None
    j = src.find("{", i)
    depth = 0
    for k in range(j, len(src)):
        if src[k] == "{":
            depth += 1
        elif src[k] == "}":
            depth -= 1
            if depth == 0:
                return src[i:k + 1]
    return None


def _run_node(fn_src, calls):
    driver = (
        fn_src
        + "\nconst calls = " + json.dumps(calls) + ";\n"
        + "process.stdout.write(JSON.stringify(calls.map(function (c) {\n"
        + "  return protocolMessage(c[0], c[1], c[2]);\n"
        + "})));\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False,
                                     encoding="utf-8") as f:
        f.write(driver)
        path = f.name
    try:
        out = subprocess.run([NODE, path], capture_output=True, text=True,
                             timeout=30)
        assert out.returncode == 0, out.stderr
        return json.loads(out.stdout)
    finally:
        os.unlink(path)


# --------------------------------------------------------------------------- #
# Source-contract gates.
# --------------------------------------------------------------------------- #
def test_naive_error_type_shape_is_gone():
    # the exact line that typeset "undefined: undefined" on a schema-422.
    assert 'data.error_type + ": " + data.message' not in _src()


def test_protocolMessage_and_detail_handling_present():
    src = _src()
    assert "function protocolMessage" in src, "no protocolMessage helper in app.js"
    # it reads FastAPI's 422 {detail:[{loc,msg},…]} shape, like worker.js.
    assert ".detail" in src and ".loc" in src and ".msg" in src


def test_three_submit_handlers_are_try_wrapped():
    # compute + big-job + feedback each gain a try/catch (baseline was 2 in the
    # catalog/param helpers), so a rejected fetch / non-JSON body reaches errDiv.
    src = _src()
    assert src.count("try {") >= 5, src.count("try {")
    assert src.count("catch (") >= 5, src.count("catch (")
    # the network catch surfaces the failure via the error box, not silence.
    assert "errDiv(" in src


# --------------------------------------------------------------------------- #
# Behaviour gate (node): the helper never yields "undefined", whatever the shape.
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not NODE, reason="node not available")
def test_protocolMessage_handles_every_error_shape(tmp_path):
    fn = _extract_fn(_src(), "protocolMessage")
    assert fn, "protocolMessage not found / not brace-balanced in app.js"

    # the REAL FastAPI 422 the server emits for an incomplete algebra spec.
    r = _client(tmp_path).post("/api/compute", json={"schema": 1})
    assert r.status_code == 422
    detail_body = r.json()
    assert "detail" in detail_body and "error_type" not in detail_body, detail_body

    calls = [
        [422, detail_body, "compute failed"],                       # FastAPI schema 422
        [422, {"error_type": "TooLarge", "message": "too big"}, "x"],  # server _error_response
        [500, None, "server error"],                                # non-JSON 500 page
        [200, {"error": {"type": "SchemaError", "message": "bad relation"}}, "y"],  # probe shape
    ]
    outs = _run_node(fn, calls)
    for o in outs:
        assert isinstance(o, str) and o and "undefined" not in o, repr(o)
    assert "algebra" in outs[0], outs[0]          # the 422 surfaced the offending field
    assert "TooLarge" in outs[1], outs[1]
    assert "server error" in outs[2], outs[2]     # the fallback rode through
    assert "SchemaError" in outs[3] or "bad relation" in outs[3], outs[3]
