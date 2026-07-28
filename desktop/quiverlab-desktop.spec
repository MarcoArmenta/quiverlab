# PyInstaller spec for the QuiverLab desktop app (one double-clickable file).
#
# Build (from the repo root, with `pip install -e ".[hpc,web]" pyinstaller`):
#     pyinstaller --noconfirm desktop/quiverlab-desktop.spec
# Output: dist/QuiverLab (dist/QuiverLab.exe on Windows). CI: .github/workflows/desktop.yml.
#
# Design (mirrors desktop/launcher.py):
# * numba/llvmlite are EXCLUDED -- the launcher pins QUIVERLAB_NO_NUMBA=1 and the
#   pure-exact kernel path is parity-gated to identical results.
# * the whole webapp/ tree ships as DATA (templates, static, vendored KaTeX,
#   i18n, precomputed manifest ride along by construction); the launcher puts the
#   bundle root on sys.path.
# * collect_submodules('quiverlab'): the engine/dispatch layers import lazily, so
#   static analysis alone would miss them.
# * an optional desktop/seed-cache.db (built by container/seed_cache.py) lands in
#   the bundle root; the launcher points QUIVERLAB_SEED_CACHE at it.

import os

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

spec_dir = os.path.dirname(os.path.abspath(SPEC))
root = os.path.dirname(spec_dir)

datas = [(os.path.join(root, "webapp"), "webapp")]
datas += collect_data_files("quiverlab")  # references.bib, zoo_catalog.json, ...
seed = os.path.join(spec_dir, "seed-cache.db")
if os.path.exists(seed):
    datas.append((seed, "."))

hiddenimports = (
    collect_submodules("quiverlab")
    + collect_submodules("webapp")
    + collect_submodules("uvicorn")   # uvicorn resolves loops/protocols dynamically
    + ["yaml", "ulid"]                # lazy imports (hpc.spec / clusterconfig; store)
)

a = Analysis(
    [os.path.join(spec_dir, "launcher.py")],
    pathex=[root],
    datas=datas,
    hiddenimports=hiddenimports,
    excludes=[
        "numba", "llvmlite",          # pure-exact path by design (see header)
        "tkinter",                    # matplotlib is Agg-only here
        "pytest", "IPython", "PyQt5", "PyQt6", "PySide2", "PySide6",
    ],
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    name="QuiverLab",
    console=True,   # the server banner IS the UI feedback; browser opens on top
    upx=False,
    strip=False,
)
