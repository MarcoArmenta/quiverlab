"""The draw-a-quiver canvas, ported from the retired docs landing page into the
webapp (served by the offline/desktop app and the deployed tier alike).

Contract under test:
  * ``webapp/static/gui/gui.js`` and ``gui.css`` are VENDORED BYTE-IDENTICAL
    copies of ``docs/gui/`` -- the canvas has one source of truth and the
    server-backed ``worker.js`` must keep speaking its exact message protocol.
  * ``/draw`` (and ``/es/draw``) render the ``#qlgui`` mount with no inline
    scripts (strict CSP), and every ``gui/…`` sibling the canvas fetches by
    page-relative URL is served under BOTH language prefixes.
  * ``/gui/manifest.json`` carries a truthful non-null wheel marker and the
    library version (gui.js refuses a null wheel).
  * ``POST /api/gui/probe`` is the server-backed ``run_build``: build-only,
    protocol payload always 200, honest error types, never a 500.
"""
import json
from pathlib import Path

from fastapi.testclient import TestClient

from webapp.server.app import create_app
from webapp.server.config import Config

ROOT = Path(__file__).resolve().parents[2]


def _client(tmp_path):
    return TestClient(create_app(Config.from_env({"QLWEB_DATA_DIR": str(tmp_path)})))


# --------------------------------------------------------------------------- #
# Vendoring parity: one canvas, two homes
# --------------------------------------------------------------------------- #
def test_gui_js_and_css_are_byte_identical_to_docs_gui():
    for name in ("gui.js", "gui.css"):
        docs = (ROOT / "docs" / "gui" / name).read_bytes()
        vendored = (ROOT / "webapp" / "static" / "gui" / name).read_bytes()
        assert docs == vendored, (
            f"{name}: webapp/static/gui/ has drifted from docs/gui/ -- "
            "edit docs/gui/ and re-copy (the canvas has ONE source of truth)")


def test_worker_shim_speaks_the_protocol_verbs():
    src = (ROOT / "webapp" / "static" / "gui" / "worker.js").read_text(encoding="utf-8")
    # The inbound commands (incl. the GitHub-#3 "random" verb) and the outbound
    # message types gui.js consumes.
    for token in ('"init"', '"run"', '"probe"', '"random"', '"calibrate"',
                  '"ready"', '"calibrated"', '"built"', '"result"',
                  '"trace"', '"artifacts"', '"done"', '"fatal"'):
        assert token in src, f"worker shim lost protocol token {token}"
    # Absolute API paths: a worker resolves relative URLs against its own
    # script URL, which would 404 under the /es mount.
    assert '"/api/compute"' in src and '"/api/gui/probe"' in src
    assert '"/api/gui/random-module"' in src


# --------------------------------------------------------------------------- #
# Pages + assets, both languages
# --------------------------------------------------------------------------- #
def test_draw_page_renders_in_both_languages(tmp_path):
    c = _client(tmp_path)
    for path in ("/draw", "/es/draw"):
        r = c.get(path)
        assert r.status_code == 200, path
        assert 'id="qlgui"' in r.text
        assert "gui/gui.js" in r.text
        # strict CSP: no inline scripts anywhere on the page
        for frag in r.text.split("<script"):
            if frag and not frag.startswith("<"):
                assert 'src="' in frag.split(">", 1)[0] or frag == r.text.split("<script")[0]


def test_nav_links_to_draw(tmp_path):
    c = _client(tmp_path)
    assert '/draw"' in c.get("/").text
    assert '/es/draw"' in c.get("/es/").text


def test_gui_assets_served_under_both_prefixes(tmp_path):
    c = _client(tmp_path)
    for prefix in ("", "/es"):
        for name, marker in (("gui.js", "qlgui"), ("gui.css", "qlgui"),
                             ("worker.js", "postMessage"),
                             ("mathjax-katex-shim.js", "renderMathInElement")):
            r = c.get(f"{prefix}/gui/{name}")
            assert r.status_code == 200, (prefix, name)
            assert marker in r.text
        presets = c.get(f"{prefix}/gui/presets.json")
        assert presets.status_code == 200
        assert isinstance(presets.json(), list) and presets.json()


def test_gui_asset_whitelist_404s_unknown_names(tmp_path):
    c = _client(tmp_path)
    assert c.get("/gui/secrets.txt").status_code == 404
    assert c.get("/gui/..%2f..%2fserver%2fapp.py").status_code == 404


def test_gui_manifest_versions_the_server_engine(tmp_path):
    import quiverlab
    m = _client(tmp_path).get("/gui/manifest.json").json()
    assert m["wheel"]                     # non-null: gui.js gates on it
    assert m["quiverlab_version"] == getattr(quiverlab, "__version__", "unknown")


# --------------------------------------------------------------------------- #
# The probe endpoint (server-backed run_build)
# --------------------------------------------------------------------------- #
def test_probe_builds_a_family_and_reports_dim(tmp_path):
    c = _client(tmp_path)
    r = c.post("/api/gui/probe", json={
        "schema": 1,
        "algebra": {"kind": "family", "family": "QuantumCI", "params": {"q": 1},
                     "field": {"kind": "GF", "p": 2, "n": 1}},
        "compute": ["hh_cohomology:0..2"]})
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is True
    assert body["dim"] > 0
    assert body["algebra"]                # human label for the "built" block


def test_probe_reports_honest_error_not_500(tmp_path):
    c = _client(tmp_path)
    r = c.post("/api/gui/probe", json={
        "schema": 1,
        "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                     "relations": [], "field": {"kind": "GF", "p": 2, "n": 1}},
        "compute": ["hh_cohomology:0..2"]})
    assert r.status_code == 200
    body = r.json()
    # kQ with a free loop and no relations is infinite-dimensional: the probe
    # must surface the library's loud refusal as a typed protocol error.
    assert body["ok"] is False
    assert body["error"]["type"]
    assert body["error"]["message"]


def test_probe_accepts_schema2_module_requests(tmp_path):
    """The live 422 regression: the canvas tags module-carrying requests as
    schema 2; the probe endpoint reads the RAW body (never ComputeRequest), so
    a module block -- or any half-finished editor state -- can never 422."""
    c = _client(tmp_path)
    r = c.post("/api/gui/probe", json={
        "schema": 2,
        "algebra": {"kind": "quiver", "vertices": [1, 2],
                     "arrows": {"a": [1, 2]}, "relations": [],
                     "field": {"kind": "GF", "p": 2, "n": 1}},
        "compute": ["dimension_vector"],
        "module": {"builtin": {"kind": "simple", "vertex": 1, "side": "right"}}})
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is True and body["dim"] == 3


def test_probe_never_422s_on_malformed_bodies(tmp_path):
    c = _client(tmp_path)
    for payload in ({"schema": 3, "algebra": {}},
                    {"algebra": {"kind": "family"}},
                    {"schema": 1},
                    ["not", "a", "dict"]):
        r = c.post("/api/gui/probe", json=payload)
        assert r.status_code == 200, payload
        body = r.json()
        assert body["ok"] is False and body["error"]["type"], payload
