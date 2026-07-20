"""The [qpa] guard behaves without GAP installed (runs on every cell).

These are marked `fast` explicitly so they run in the normal matrix even though
they live under tests/qpa/ (the path default would be `qpa`)."""
import pytest

from quiverlab import Quiver, GF, QpaUnavailableError
from quiverlab.qpa import gap_available

pytestmark = pytest.mark.fast


def test_gap_available_is_boolean_and_cheap():
    assert isinstance(gap_available(), bool)      # no exception when GAP absent


def test_crosscheck_raises_cleanly_without_backend():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(2))
    if gap_available():
        pytest.skip("GAP present; guard-absence path not exercised here")
    with pytest.raises(QpaUnavailableError) as e:
        A.crosscheck("hochschild", 2)
    assert "quiverlab[qpa]" in str(e.value) or "WSL" in str(e.value)
