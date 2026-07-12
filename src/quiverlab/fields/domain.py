"""The Domain protocol: exact field arithmetic behind one interface (spec §5, component 1)."""
from fractions import Fraction

from quiverlab.errors import ExactnessError

_FLOAT_HINT = "quiverlab is exact-only: write '1/3' or Fraction(1, 3), never 0.333"


def reject_inexact(x):
    """Loud gate for obviously non-exact Python inputs (spec D3)."""
    if isinstance(x, bool):
        raise ExactnessError(f"{x!r} is a bool, not a scalar", hint="use 0 or 1")
    if isinstance(x, (float, complex)):
        raise ExactnessError(f"{x!r} is a float", hint=_FLOAT_HINT)
    if isinstance(x, str) and any(
        ch == "." and (i + 1 < len(x) and x[i + 1].isdigit() or i > 0 and x[i - 1].isdigit())
        for i, ch in enumerate(x)
    ):
        raise ExactnessError(f"decimal literal in {x!r}", hint=_FLOAT_HINT)
    return x


class Domain:
    """Abstract exact field. Elements are plain Python objects; ops go through the domain."""

    name: str = "?"
    characteristic: int = 0

    # -- construction hooks used by Algebra builders ------------------------
    def parse_entry(self, x):
        """Pre-parse a raw user entry (overridden by CC). Default: pass through."""
        return reject_inexact(x)

    def make_domain(self, entries):
        """Return the concrete Domain for these entries (overridden by CC)."""
        return self

    # -- required field interface ------------------------------------------
    def coerce(self, x):
        raise NotImplementedError

    def zero(self):
        raise NotImplementedError

    def one(self):
        raise NotImplementedError

    def add(self, a, b):
        raise NotImplementedError

    def neg(self, a):
        raise NotImplementedError

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        raise NotImplementedError

    def inv(self, a):
        raise NotImplementedError

    def is_zero(self, a) -> bool:
        raise NotImplementedError

    def eq(self, a, b) -> bool:
        return self.is_zero(self.sub(a, b))

    def to_str(self, a) -> str:
        return str(a)

    def __repr__(self):
        return self.name


def parse_rational(x) -> Fraction:
    """int | Fraction | 'a/b' string -> Fraction, loudly exact."""
    reject_inexact(x)
    if isinstance(x, int):
        return Fraction(x)
    if isinstance(x, Fraction):
        return x
    if isinstance(x, str):
        try:
            return Fraction(x.strip())
        except (ValueError, ZeroDivisionError) as exc:
            raise ExactnessError(
                f"cannot read {x!r} as an exact rational", hint="use forms like '2' or '-1/3'"
            ) from exc
    raise ExactnessError(f"cannot read {type(x).__name__} {x!r} as an exact scalar",
                         hint=_FLOAT_HINT)
