"""GF(p): exact prime-field arithmetic on plain ints (spec §4.3)."""
import sympy

from quiverlab.errors import FieldError
from quiverlab.fields.domain import Domain, parse_rational


class PrimeField(Domain):
    characteristic: int

    def __init__(self, p: int):
        if not (isinstance(p, int) and not isinstance(p, bool) and sympy.isprime(p)):
            raise FieldError(f"GF({p!r}): {p!r} is not a prime",
                             hint="use GF(p) with p prime, or GF(p**n) for prime powers")
        self.p = p
        self.characteristic = p
        self.name = f"GF({p})"

    def coerce(self, x):
        if isinstance(x, int) and not isinstance(x, bool):
            return x % self.p
        q = parse_rational(x)
        if q.denominator % self.p == 0:
            raise FieldError(
                f"{q} has denominator divisible by {self.p}, so it is not in GF({self.p})",
                hint="clear denominators in your presentation or choose another characteristic",
            )
        return (q.numerator * pow(q.denominator, -1, self.p)) % self.p

    def zero(self):
        return 0

    def one(self):
        return 1 % self.p

    def add(self, a, b):
        return (a + b) % self.p

    def neg(self, a):
        return (-a) % self.p

    def mul(self, a, b):
        return (a * b) % self.p

    def inv(self, a):
        if a % self.p == 0:
            raise ZeroDivisionError(f"0 has no inverse in GF({self.p})")
        return pow(a, -1, self.p)

    def is_zero(self, a):
        return a % self.p == 0
