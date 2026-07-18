"""quiverlab exceptions (spec §7). Every message states the problem and a fix-it hint."""


class QuiverlabError(Exception):
    def __init__(self, message: str, hint: str | None = None):
        self.hint = hint
        super().__init__(message if hint is None else f"{message}  [hint: {hint}]")


class ExactnessError(QuiverlabError):
    """A float (or other non-exact input) tried to enter quiverlab."""


class FieldError(QuiverlabError):
    """Unsupported field, or an entry that does not live in the stated field."""


class RelationError(QuiverlabError):
    """A relation string is malformed, non-composable, or not parallel."""


class AdmissibilityError(QuiverlabError):
    """The ideal cannot be certified admissible."""


class NotFiniteDimensionalError(QuiverlabError):
    """The presented algebra is (or cannot be certified) finite-dimensional."""


class DepthLimitError(QuiverlabError):
    """A guard stopped a computation; the certified range is stated in the message."""
