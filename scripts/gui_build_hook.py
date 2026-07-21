"""mkdocs build hook (Plan 10): package the GUI engine payload into the site.

on_post_build writes, under <site>/gui/:
  quiverlab-<version>-py3-none-any.whl   pip wheel of THIS checkout (the GUI
                                         always runs the exact code being documented)
  manifest.json                          {"schema": 1, "wheel": ..., "quiverlab_version": ...}
  presets.json                           curated examples EXTRACTED via the library
                                         (A.quiver / A.relations), never hand-written

QLGUI_SKIP_WHEEL=1 skips only the (slow) wheel build for fast `mkdocs serve`
iterations; manifest.json then carries {"wheel": null} and the GUI shows an
"engine payload not built" chip instead of loading Pyodide."""
import json
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

REPO = pathlib.Path(__file__).resolve().parents[1]


def _preset_algebras():
    import quiverlab as ql
    ql.verbose = False
    entries = [
        ("Kronecker quiver (CC)",
         ql.Quiver(vertices=[1, 2], arrows={"a": (1, 2), "b": (1, 2)})
           .algebra(field=ql.CC),
         {"kind": "CC"}),
        ("A3 path, zero relation a*b (CC)",
         ql.Quiver(vertices=[1, 2, 3], arrows={"a": (1, 2), "b": (2, 3)})
           .algebra(relations=["a*b"], field=ql.CC),
         {"kind": "CC"}),
        ("Truncated polynomial GF(2)[x]/(x^3)",
         ql.truncated_polynomial(n=3, field=ql.GF(2)),
         {"kind": "GF", "p": 2, "n": 1}),
    ]
    for i, A in enumerate(ql.zoo(dim_max=12)):
        if i >= 2:
            break
        entries.append(("Exact zoo #%d — dim %d (CC)" % (i + 1, A.dim), A,
                        {"kind": "CC"}))
    return entries


def generate_presets():
    return [{"label": label,
             "vertices": list(A.quiver.vertices),
             "arrows": {k: list(v) for k, v in A.quiver.arrows.items()},
             "relations": [str(r) for r in A.relations],
             "field": field_spec}
            for label, A, field_spec in _preset_algebras()]


def build_wheel(gui_dir):
    """pip-wheel this checkout into gui_dir; return the wheel filename."""
    gui_dir = pathlib.Path(gui_dir)
    # pip builds in-tree: a stale build/ from a previous run collides
    # ([Errno 17] on dist-info), so every rebuild must start clean.
    shutil.rmtree(REPO / "build", ignore_errors=True)
    with tempfile.TemporaryDirectory() as td:
        proc = subprocess.run(
            [sys.executable, "-m", "pip", "wheel", "--no-deps", "-w", td, str(REPO)],
            capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError("GUI wheel build failed:\n%s" % proc.stderr)
        wheels = list(pathlib.Path(td).glob("quiverlab-*.whl"))
        if len(wheels) != 1:
            raise RuntimeError("expected exactly one quiverlab wheel, got %r"
                               % ([w.name for w in wheels],))
        shutil.copy2(wheels[0], gui_dir / wheels[0].name)
        return wheels[0].name


def on_post_build(config):
    import quiverlab
    gui = pathlib.Path(config["site_dir"]) / "gui"
    gui.mkdir(parents=True, exist_ok=True)
    wheel = None
    if not os.environ.get("QLGUI_SKIP_WHEEL"):
        wheel = build_wheel(gui)
    (gui / "manifest.json").write_text(json.dumps(
        {"schema": 1, "wheel": wheel, "quiverlab_version": quiverlab.__version__}))
    (gui / "presets.json").write_text(json.dumps(generate_presets(), indent=1))
