"""QuiverLab desktop launcher -- the PyInstaller onefile entry point.

One double-clickable file: start the offline GUI (webapp/server/offline.py --
embedded worker, no network, vendored KaTeX) and open the browser on it. This
file only wires the frozen environment; ALL behavior lives in the library:

* pure-exact kernels: ``QUIVERLAB_NO_NUMBA=1`` -- numba is deliberately NOT
  bundled. Engine parity between the numba and pure paths is test-gated
  (results identical, the pure path is just slower), and leaving LLVM out keeps
  the binary small and the freeze reliable on every OS.
* ``webapp/`` ships inside the bundle (PyInstaller ``--add-data``); the bundle
  root goes on ``sys.path`` before the import below.
* port: 8000 when free, else an ephemeral one -- two running copies must not
  collide. ``QUIVERLAB_DESKTOP_PORT`` pins it (CI smoke).
* the browser opens once the server actually answers;
  ``QUIVERLAB_DESKTOP_NO_BROWSER=1`` disables that (CI smoke, headless use).

Not a module of the library: lives in ``desktop/``, never imported by
``quiverlab`` or ``webapp``.
"""
import multiprocessing
import os
import socket
import sys
import threading
import time
import webbrowser

# Before ANY quiverlab import: the desktop app runs the pure-exact kernel path.
os.environ.setdefault("QUIVERLAB_NO_NUMBA", "1")
# Headless matplotlib (quiverlab.viz imports it; there is no display loop here).
os.environ.setdefault("MPLBACKEND", "Agg")


def _bundle_root() -> str:
    """The directory holding the bundled ``webapp/`` tree: PyInstaller's
    ``_MEIPASS`` when frozen, the repo root when run from a source checkout."""
    frozen = getattr(sys, "_MEIPASS", None)
    if frozen:
        return frozen
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


_ROOT = _bundle_root()
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

_SEED = os.path.join(_ROOT, "seed-cache.db")
if os.path.exists(_SEED):
    os.environ.setdefault("QUIVERLAB_SEED_CACHE", _SEED)

from webapp.server.offline import serve_offline  # noqa: E402  (path set above)


def _pick_port() -> int:
    env = os.environ.get("QUIVERLAB_DESKTOP_PORT")
    if env:
        return int(env)
    for candidate in (8000, 0):  # 0 -> OS-assigned free port
        s = socket.socket()
        try:
            s.bind(("127.0.0.1", candidate))
            port = s.getsockname()[1]
            s.close()
            return port
        except OSError:
            s.close()
    raise RuntimeError("no free TCP port")


def _open_when_ready(port: int, timeout_s: int = 120) -> None:
    """Poll until the server answers, then open the browser tab exactly once."""
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=1):
                break
        except OSError:
            time.sleep(0.3)
    else:
        return
    # Land on the draw-a-quiver-and-compute page right away (Marco 2026-08-03);
    # the catalog stays one click away at "/".
    webbrowser.open(f"http://localhost:{port}/draw")


def main() -> int:
    port = _pick_port()
    if os.environ.get("QUIVERLAB_DESKTOP_NO_BROWSER") != "1":
        threading.Thread(target=_open_when_ready, args=(port,),
                         daemon=True, name="ql-desktop-browser").start()
    try:
        serve_offline(port=port, open_hint=True)
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    # The worker runs each job in a multiprocessing *spawn* child, which in a
    # frozen app re-executes this very binary -- without freeze_support() every
    # job child would boot a second server instead of running the job.
    multiprocessing.freeze_support()
    sys.exit(main())
