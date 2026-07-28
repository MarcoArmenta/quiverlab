# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for the quiverlab desktop app (download -> double-click ->
# the offline GUI opens in the browser). Build from the repo root:
#
#   pip install -e ".[fast,web,hpc]" pyinstaller
#   pyinstaller desktop/quiverlab.spec
#
# Design: `webapp/` is EXCLUDED from module analysis and shipped verbatim as
# data files; the launcher puts the bundle root on sys.path so `import webapp`
# resolves to real files (its __file__-relative templates/static/precomputed
# paths then work unchanged -- the same trick as the container's
# PYTHONPATH=/app). The library and the webapp's third-party deps are pulled
# in explicitly since the excluded webapp is never traced.
import pathlib

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

ROOT = pathlib.Path(SPECPATH).resolve().parent      # repo root (spec lives in desktop/)

datas = [(str(ROOT / "webapp"), "webapp")]
datas += collect_data_files("quiverlab")            # packaged bib / catalogs

hiddenimports = (
    collect_submodules("quiverlab")                 # lazy imports everywhere
    + collect_submodules("uvicorn")                 # dynamic loop/protocol picks
    + collect_submodules("fastapi")                 # e.g. fastapi.staticfiles
    + collect_submodules("starlette")               # fastapi's runtime backend
    + ["jinja2", "pydantic", "ulid", "yaml", "sqlite3",
       "email.mime.text", "multiprocessing"]
)

a = Analysis(
    [str(ROOT / "desktop" / "launcher.py")],
    pathex=[str(ROOT)],
    datas=datas,
    hiddenimports=hiddenimports,
    excludes=["webapp",                             # shipped as data, see above
              "matplotlib", "tkinter", "IPython", "pytest"],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name="quiverlab",
    debug=False,
    strip=False,
    upx=False,
    console=True,          # keep the banner/log visible; closing it stops the app
    onefile=True,
)
