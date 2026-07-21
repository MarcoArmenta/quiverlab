import importlib.util
import json
import pathlib

import pytest

HOOK_PATH = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "gui_build_hook.py"


def _hook():
    spec = importlib.util.spec_from_file_location("gui_build_hook", HOOK_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_presets_round_trip_through_runner(runner):
    presets = _hook().generate_presets()
    assert len(presets) >= 4
    labels = [p["label"] for p in presets]
    assert any("Kronecker" in l for l in labels)
    for p in presets:
        req = {"schema": 1,
               "algebra": {"kind": "quiver", "vertices": p["vertices"],
                           "arrows": p["arrows"], "relations": p["relations"],
                           "field": p["field"]},
               "compute": []}
        out = json.loads(runner.run_build(json.dumps(req)))
        assert out["ok"], (p["label"], out)


@pytest.mark.deep     # shells out to a full PEP-517 wheel build (~10s+)
def test_build_wheel_and_manifest(tmp_path):
    hook = _hook()
    gui = tmp_path / "site" / "gui"
    gui.mkdir(parents=True)
    name = hook.build_wheel(gui)
    assert name.startswith("quiverlab-") and name.endswith("-py3-none-any.whl")
    assert (gui / name).exists()
    hook.on_post_build({"site_dir": str(tmp_path / "site")})
    manifest = json.loads((gui / "manifest.json").read_text())
    assert manifest["schema"] == 1 and manifest["wheel"] == name
    assert manifest["quiverlab_version"]
    presets = json.loads((gui / "presets.json").read_text())
    assert isinstance(presets, list) and presets


def test_skip_wheel_env(tmp_path, monkeypatch):
    hook = _hook()
    monkeypatch.setenv("QLGUI_SKIP_WHEEL", "1")
    hook.on_post_build({"site_dir": str(tmp_path / "site")})
    manifest = json.loads((tmp_path / "site" / "gui" / "manifest.json").read_text())
    assert manifest["wheel"] is None
    assert (tmp_path / "site" / "gui" / "presets.json").exists()
