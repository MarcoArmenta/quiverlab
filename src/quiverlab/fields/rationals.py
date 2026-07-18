"""Exact rationals: the char-0 workhorse domain (elements: fractions.Fraction)."""
from fractions import Fraction

from quiverlab.fields.domain import Domain, parse_rational


class RationalField(Domain):
    name = "QQ"
    characteristic = 0

    def coerce(self, x):
        if isinstance(x, Fraction) and not isinstance(x, bool):
            return x
        return parse_rational(x)

    def zero(self):
        return Fraction(0)

    def one(self):
        return Fraction(1)

    def add(self, a, b):
        return a + b

    def neg(self, a):
        return -a

    def mul(self, a, b):
        return a * b

    def inv(self, a):
        return Fraction(1) / a

    def is_zero(self, a):
        return a == 0


QQ = RationalField()
