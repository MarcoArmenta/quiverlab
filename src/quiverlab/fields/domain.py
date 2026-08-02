"""The Domain protocol: exact field arithmetic behind one interface (spec §5, component 1)."""
import re
from fractions import Fraction

from quiverlab.errors import ExactnessError

_FLOAT_HINT = "quiverlab is exact-only: write '1/3' or Fraction(1, 3), never 0.333"
# digit-E-digit run: catches exponent notation like "15e-1" that Fraction() silently parses.
_SCI_NOTATION = re.compile(r"\d[eE][+-]?\d")


def reject_inexact(x):
    """Loud gate for obviously non-exact Python inputs (spec D3)."""
    if isinstance(x, bool):
        raise ExactnessError(f"{x!r} is a bool, not a scalar", hint="use 0 or 1")
    if isinstance(x, (float, complex)):
        raise ExactnessError(f"{x!r} is a float", hint=_FLOAT_HINT)
    if isinstance(x, str):
        if any(
            ch == "." and (i + 1 < len(x) and x[i + 1].isdigit() or i > 0 and x[i - 1].isdigit())
            for i, ch in enumerate(x)
        ):
            raise ExactnessError(f"decimal literal in {x!r}", hint=_FLOAT_HINT)
        if _SCI_NOTATION.search(x):
            raise ExactnessError(f"scientific notation in {x!r}", hint=_FLOAT_HINT)
    return x


def exact_rational(x):
    """If ``x`` carries the exact rational protocol -- an integer ``.numerator``
    and ``.denominator`` -- return the pair ``(p, q)`` as Python ints; otherwise
    ``None``. This admits values by PROTOCOL, not by a name-list: sympy's own QQ
    elements (``sympy.external.pythonmpq.PythonMPQ`` in a pure install, ``gmpy2``
    ``mpq``/``mpz`` when the C backend is present), ``fractions.Fraction``,
    ``int``, and ``sympy`` ``Rational``/``Integer`` all qualify. Such a value is
    an EXACT rational -- refusing it as "not an exact scalar" is a false refusal
    (these leak out of sympy Matrix/domain arithmetic over QQ and re-enter the
    exact-scalar readers).

    Float/complex/Decimal have no integer num/den pair and miss. ``bool`` IS a
    ``numbers.Rational`` and would qualify, so callers MUST reject it upstream
    via :func:`reject_inexact` (house style) BEFORE consulting this helper."""
    num = getattr(x, "numerator", None)
    den = getattr(x, "denominator", None)
    if num is None or den is None:
        return None
    if isinstance(num, (float, complex)) or isinstance(den, (float, complex)):
        return None
    try:
        p, q = int(num), int(den)
    except (TypeError, ValueError):
        return None
    if p != num or q != den:          # a non-integral num/den that int() truncated
        return None
    return (p, q)


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
    """int | Fraction | 'a/b' string -> Fraction, loudly exact. Also accepts any
    exact rational carried by the numerator/denominator protocol (PythonMPQ,
    gmpy2 mpq/mpz, sympy Rational/Integer) -- how sympy's QQ domain represents its
    own elements, which leak back into the readers through Matrix arithmetic."""
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
    pq = exact_rational(x)
    if pq is not None:
        return Fraction(*pq)
    raise ExactnessError(f"cannot read {type(x).__name__} {x!r} as an exact scalar",
                         hint=_FLOAT_HINT)
