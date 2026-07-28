"""quiverlab desktop launcher (Plan: download button -> double-click -> GUI).

The PyInstaller entry point for the standalone desktop app: start the offline
webapp (the same `quiverlab-hpc gui` server -- embedded worker, loopback-only,
no network) and open the user's browser at it once it answers. Frozen-bundle
layout: the whole `webapp/` source tree ships as DATA files next to the
executable payload and `sys.path` gains that directory, so `import webapp`
resolves to real files and its template/static paths work unchanged.

Not part of the wheel -- built by desktop/quiverlab.spec via
.github/workflows/desktop.yml into per-OS binaries attached to releases.
"""
import multiprocessing

# MUST run before anything else: the offline worker computes each job in a
# resource-capped ``spawn`` child, and in a frozen bundle spawn re-executes
# THIS executable -- without freeze_support() every child would boot the
# launcher (and its server) again instead of the multiprocessing bootstrap,
# forking a process storm. freeze_support() intercepts the child argv.
multiprocessing.freeze_support()

import os
import socket
import sys
import threading
import time
import urllib.request
import webbrowser


def _bundle_root():
    # onefile: payload extracted to _MEIPASS; onedir: the app directory.
    return getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))


def _free_port(preferred: int = 8000) -> int:
    for port in (preferred, 0):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(("127.0.0.1", port))
                return s.getsockname()[1]
        except OSError:
            continue
    return preferred


def _open_when_ready(url: str, timeout_s: int = 90) -> None:
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        try:
            with urllib.request.urlopen(url, timeout=2):
                webbrowser.open(url)
                return
        except Exception:
            time.sleep(0.5)
    # server never answered -- the banner in the console still shows the URL.


def main() -> int:
    root = _bundle_root()
    if root not in sys.path:
        sys.path.insert(0, root)
    try:
        from webapp.server.offline import serve_offline
    except ImportError as exc:
        print("quiverlab desktop: bundled webapp not found (%s)" % exc,
              file=sys.stderr)
        return 70
    port = _free_port(8000)
    url = f"http://127.0.0.1:{port}/"
    threading.Thread(target=_open_when_ready, args=(url,), daemon=True).start()
    try:
        serve_offline(port=port, open_hint=True)
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
