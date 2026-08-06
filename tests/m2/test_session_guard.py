"""Session guards for the Macaulay2 bridge -- run everywhere, no M2 needed."""
import os
import pytest

from quiverlab.errors import M2UnavailableError
from quiverlab.m2 import session

pytestmark = pytest.mark.fast


def test_unavailable_raises_with_hint(monkeypatch):
    monkeypatch.setattr(session, "_which_m2", lambda: None)
    session.m2_available.cache_clear()
    with pytest.raises(M2UnavailableError) as e:
        session.require_m2()
    assert "brew" in str(e.value) or "apt" in str(e.value)
    session.m2_available.cache_clear()


def test_skip_predicate_env_override(monkeypatch):
    monkeypatch.setattr(session, "_which_m2", lambda: None)
    session.m2_available.cache_clear()
    monkeypatch.delenv("QUIVERLAB_REQUIRE_M2", raising=False)
    assert session.should_skip_m2() is True
    monkeypatch.setenv("QUIVERLAB_REQUIRE_M2", "1")
    assert session.should_skip_m2() is False   # CI: absence must FAIL, not skip
    session.m2_available.cache_clear()


def test_import_does_not_probe(monkeypatch):
    # importing quiverlab.m2 must not shell out (mirror qpa laziness)
    import importlib
    import quiverlab.m2
    importlib.reload(quiverlab.m2)   # no error even with no M2 on PATH
