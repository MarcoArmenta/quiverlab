"""Cross-runner + app-shell pieces of Marco's 2026-08-03 report pass.

  * both runners (``quiverlab.hpc.spec`` and the Pyodide twin
    ``docs/gui/runner.py``) stamp the SAME ``resolved`` provenance on Ext/Tor
    blocks -- which module was resolved, by which resolution;
  * the app lands on the draw-a-quiver page right away: the desktop launcher
    opens ``/draw``, and the offline banner points there.
"""
import json
import tempfile

# NO oracle-class marker: tests/webapp/ collects only with the [web] extra, and the
# Plan-32 audit requires the audited class counts to be environment-independent.

_REQ = {
    "schema": 2,
    "algebra": {"kind": "quiver", "vertices": [1], "arrows": {"x": [1, 1]},
                "relations": ["x*x*x"], "field": {"kind": "GF", "p": 5, "n": 1}},
    "compute": ["ext:0..2", "tor:0..2"],
    "module": {"builtin": {"kind": "simple", "vertex": 1}, "side": "right"},
    "ext_target": {"builtin": {"kind": "simple", "vertex": 1}, "side": "right"},
    "tor_target": {"builtin": {"kind": "simple", "vertex": 1}, "side": "left"},
}


def _server_results(req):
    from quiverlab.hpc import spec
    with tempfile.TemporaryDirectory() as d:
        return spec.run(req, d)["results"]


def _pyodide_results(req):
    import importlib.util
    import pathlib
    path = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"
    s = importlib.util.spec_from_file_location("_gui_runner_m0803", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)
    assert json.loads(mod.run_build(json.dumps(req)))["ok"]
    out = {}
    for item in req["compute"]:
        r = json.loads(mod.compute_one(item))
        assert r["ok"], r
        out[item.split(":")[0]] = r["block"]
    return out


def test_both_runners_stamp_resolved_provenance_identically():
    server = _server_results(_REQ)
    gui = _pyodide_results(_REQ)
    for kind in ("ext", "tor"):
        assert server[kind]["resolved"] == gui[kind]["resolved"]
        r = server[kind]["resolved"]
        assert r["module"] == "M" and r["side"] == "right"
        assert "minimal projective" in r["resolution"]
        # ... and the SAME displayed resolution of M (Marco 2026-08-03 pass 2)
        assert server[kind]["resolution"] == gui[kind]["resolution"]
        assert server[kind]["resolution"]["summands"]


def test_offline_banner_points_at_draw(tmp_path):
    from webapp.server.offline import (_banner_lines, build_offline_config,
                                       detect_resources, runtime_caps)
    res = detect_resources()
    cfg = build_offline_config(tmp_path, res)
    caps = runtime_caps(cfg, res)
    lines = _banner_lines(8000, cfg, caps, open_hint=True)
    assert any("http://localhost:8000/draw" in ln for ln in lines)


def test_desktop_launcher_opens_draw(monkeypatch):
    import importlib.util
    import pathlib
    path = pathlib.Path(__file__).resolve().parents[2] / "desktop" / "launcher.py"
    s = importlib.util.spec_from_file_location("_ql_launcher_m0803", path)
    mod = importlib.util.module_from_spec(s)
    s.loader.exec_module(mod)

    opened = []

    class _Conn:
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    monkeypatch.setattr(mod.socket, "create_connection", lambda *a, **k: _Conn())
    monkeypatch.setattr(mod.webbrowser, "open", lambda url: opened.append(url))
    mod._open_when_ready(4321, timeout_s=2)
    assert opened == ["http://localhost:4321/draw"]
