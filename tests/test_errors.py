import pytest
import quiverlab
from quiverlab.errors import (
    QuiverlabError, ExactnessError, FieldError, RelationError,
    AdmissibilityError, NotFiniteDimensionalError, DepthLimitError,
)


def test_hierarchy():
    for cls in (ExactnessError, FieldError, RelationError, AdmissibilityError,
                NotFiniteDimensionalError, DepthLimitError):
        assert issubclass(cls, QuiverlabError)


def test_reexported_from_package():
    assert quiverlab.ExactnessError is ExactnessError


def test_messages_carry_hint():
    err = ExactnessError("0.5 is a float", hint="write '1/2' instead")
    assert "0.5 is a float" in str(err) and "write '1/2' instead" in str(err)
