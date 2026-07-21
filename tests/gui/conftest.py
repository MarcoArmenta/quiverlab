"""Load docs/gui/runner.py as a FRESH module per test (clean _state), exactly
the file the browser gets — no copies, no sys.path games."""
import importlib.util
import pathlib

import pytest

RUNNER_PATH = pathlib.Path(__file__).resolve().parents[2] / "docs" / "gui" / "runner.py"


@pytest.fixture()
def runner():
    spec = importlib.util.spec_from_file_location("gui_runner", RUNNER_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
