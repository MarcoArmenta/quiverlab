import subprocess
import sys

_PROBE = r"""
import sys, importlib.abc, importlib.machinery
BANNED = {"fastapi", "uvicorn", "jinja2", "ulid", "pydantic"}

class _Blocker(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path, target=None):
        top = name.split(".")[0]
        if top in BANNED:
            raise ModuleNotFoundError(f"blocked {name} (simulated uninstall)")
        return None

sys.meta_path.insert(0, _Blocker())
import quiverlab  # must succeed with the web stack unavailable
assert not (set(sys.modules) & BANNED), f"library imported web deps: {set(sys.modules) & BANNED}"
print("ok")
"""


def test_library_imports_without_web_stack():
    proc = subprocess.run([sys.executable, "-c", _PROBE],
                          capture_output=True, text=True, encoding="utf-8")
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip().endswith("ok")
