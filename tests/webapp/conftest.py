"""Collection guard: the webapp suite needs the [web] extra.

Jobs that install it (the CI fast matrix and deep legs) run these tests for
real; any context without it (the QPA cross-check job, a bare `pip install -e
".[dev]"` checkout) must SKIP the directory at collection instead of erroring
with ModuleNotFoundError -- pytest imports every test module during collection
even when -m deselects them all.
"""
import pytest

pytest.importorskip("fastapi", reason="webapp tests need the [web] extra")
