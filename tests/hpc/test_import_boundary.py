"""Import-boundary contract (Plan 28): importing ``quiverlab.hpc`` (and its CLI /
report / resources / deepen_runner submodules) must NOT pull the web stack
(fastapi/uvicorn/jinja2/pydantic/ulid) or the ``webapp`` package, and must NOT
eagerly import the sanctioned engine internal ``quiverlab.engine.deepen`` (that
import is lazy, inside ``run_deepen``). Base ``import quiverlab`` also stays clean.
Run in a subprocess so the assertion sees a pristine ``sys.modules``."""
import subprocess
import sys

_PROBE = r"""
import sys
BANNED = {"fastapi", "uvicorn", "jinja2", "pydantic", "ulid", "webapp"}

import quiverlab.hpc
import quiverlab.hpc.cli
import quiverlab.hpc.report
import quiverlab.hpc.resources
import quiverlab.hpc.deepen_runner
import quiverlab.hpc.spec

tops = {m.split(".")[0] for m in sys.modules}
leaked = BANNED & tops
assert not leaked, f"quiverlab.hpc pulled banned modules: {leaked}"
assert "quiverlab.engine.deepen" not in sys.modules, \
    "engine.deepen imported eagerly (must stay lazy in run_deepen)"
print("ok")
"""

_BASE_PROBE = r"""
import sys, importlib.abc
BANNED = {"fastapi", "uvicorn", "jinja2", "pydantic", "ulid", "webapp"}

class _Blocker(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path, target=None):
        if name.split(".")[0] in BANNED:
            raise ModuleNotFoundError(f"blocked {name}")
        return None

sys.meta_path.insert(0, _Blocker())
import quiverlab
import quiverlab.hpc          # must import with the web stack unavailable
assert not (set(sys.modules) & BANNED)
print("ok")
"""


def test_hpc_import_pulls_no_web_or_engine_deepen():
    proc = subprocess.run([sys.executable, "-c", _PROBE],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip().endswith("ok")


def test_hpc_imports_with_web_stack_uninstalled():
    proc = subprocess.run([sys.executable, "-c", _BASE_PROBE],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip().endswith("ok")
