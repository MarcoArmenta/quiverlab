"""The ``gui`` verb must work no matter where it is invoked from.

``webapp`` is deliberately NOT packaged into the wheel; it lives next to
``src/`` in a source checkout and next to the editable install inside the
container image (``/app``). A console script never puts the cwd on ``sys.path``
and the PEP-660 editable hook maps only ``quiverlab`` -- so ``_cmd_gui`` locates
the checkout root itself (``_ensure_webapp_on_path``). This is the regression
that shipped in the Plan-28 image: ``quiverlab-hpc gui`` died with
``No module named 'webapp'`` from every cwd, including ``WORKDIR /app``.

The container additionally needs a non-loopback bind (``docker run -p`` cannot
reach 127.0.0.1 inside the container's netns), hence the ``--host`` flag with
the ``QUIVERLAB_GUI_HOST`` env default the image sets to 0.0.0.0.
"""
import os
import subprocess
import sys

import pytest


def _run(code: str, cwd, env_extra=None) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    if env_extra:
        env.update(env_extra)
    return subprocess.run([sys.executable, "-c", code], cwd=str(cwd),
                          env=env, capture_output=True, text=True)


def test_ensure_webapp_on_path_from_foreign_cwd(tmp_path):
    """From a cwd with no ``webapp`` in sight, the helper finds the checkout
    root via ``quiverlab.__file__`` and makes ``import webapp`` succeed."""
    proc = _run(
        "from quiverlab.hpc.cli import _ensure_webapp_on_path\n"
        "_ensure_webapp_on_path()\n"
        "import webapp\n"
        "print('ok')\n",
        cwd=tmp_path,
    )
    assert proc.returncode == 0, proc.stderr
    assert "ok" in proc.stdout


def test_ensure_webapp_on_path_is_idempotent(tmp_path):
    proc = _run(
        "import sys\n"
        "from quiverlab.hpc.cli import _ensure_webapp_on_path\n"
        "_ensure_webapp_on_path()\n"
        "_ensure_webapp_on_path()\n"
        "import webapp\n"
        "n = sum(1 for p in sys.path if p == sys.path[0])\n"
        "assert n == 1, sys.path\n"
        "print('ok')\n",
        cwd=tmp_path,
    )
    assert proc.returncode == 0, proc.stderr
    assert "ok" in proc.stdout


def test_gui_parser_host_default_and_flag():
    from quiverlab.hpc.cli import _build_parser
    args = _build_parser().parse_args(["gui"])
    assert args.host == "127.0.0.1"
    args = _build_parser().parse_args(["gui", "--host", "0.0.0.0"])
    assert args.host == "0.0.0.0"


def test_gui_parser_host_env_default(tmp_path):
    """The container image sets QUIVERLAB_GUI_HOST=0.0.0.0 so a plain
    ``docker run -p 8000:8000 ... gui`` is reachable without extra flags."""
    proc = _run(
        "from quiverlab.hpc.cli import _build_parser\n"
        "args = _build_parser().parse_args(['gui'])\n"
        "assert args.host == '0.0.0.0', args.host\n"
        "print('ok')\n",
        cwd=tmp_path,
        env_extra={"QUIVERLAB_GUI_HOST": "0.0.0.0"},
    )
    assert proc.returncode == 0, proc.stderr
    assert "ok" in proc.stdout


def test_serve_offline_accepts_host_kwarg():
    pytest.importorskip("fastapi")
    import inspect

    from quiverlab.hpc.cli import _ensure_webapp_on_path
    _ensure_webapp_on_path()
    from webapp.server.offline import serve_offline
    assert "host" in inspect.signature(serve_offline).parameters
