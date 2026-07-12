"""GF(p^n): exact arithmetic on coefficient tuples modulo a monic irreducible
polynomial over GF(p) (spec §4.3). Elements: tuple[int, ...] of length n,
little-endian."""
import itertools

import sympy

from quiverlab.errors import FieldError
from quiverlab.fields.conway import CONWAY
from quiverlab.fields.domain import Domain, parse_rational
from quiverlab.fields.primefield import PrimeField


def _poly_trim(c):
    while c and c[-1] == 0:
        c = c[:-1]
    return c


def _poly_mod(a, mod, p):
    """a mod (monic) mod, coefficients mod p."""
    a = [x % p for x in a]
    dm = len(mod) - 1
    while len(_poly_trim(a)) - 1 >= dm:
        a = _poly_trim(a)
        d = len(a) - 1
        lead = a[-1]
        for i, c in enumerate(mod):
            a[d - dm + i] = (a[d - dm + i] - lead * c) % p
    return _poly_trim(a)


def _poly_mul(a, b, p):
    out = [0] * (len(a) + len(b) - 1) if a and b else []
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] = (out[i + j] + x * y) % p
    return _poly_trim(out)


def poly_is_irreducible(coeffs, p):
    """Trial division by every monic polynomial of degree <= deg/2 over GF(p).
    Fine for the bundled table sizes; loudly the wrong tool for huge q."""
    f = _poly_trim([c % p for c in coeffs])
    n = len(f) - 1
    if n < 1 or f[-1] != 1:
        return False
    for d in range(1, n // 2 + 1):
        for tail in itertools.product(range(p), repeat=d):
            g = list(tail) + [1]  # monic degree d
            if not _poly_mod(f, g, p):
                return False
    return True


class FiniteField(Domain):
    def __init__(self, p: int, n: int, modulus=None):
        if modulus is None:
            if (p, n) not in CONWAY:
                raise FieldError(
                    f"GF({p}^{n}) is beyond the bundled modulus table",
                    hint=f"pass an irreducible monic modulus: GF({p}**{n}, modulus=[c0, ..., 1])",
                )
            modulus = CONWAY[(p, n)]
        modulus = [c % p for c in modulus]
        if len(modulus) != n + 1 or modulus[-1] != 1:
            raise FieldError(f"modulus must be monic of degree {n} over GF({p})",
                             hint="little-endian coefficients [c0, ..., 1]")
        if not poly_is_irreducible(modulus, p):
            raise FieldError(f"modulus {modulus} is reducible over GF({p})",
                             hint="pick an irreducible polynomial (see Lübeck's Conway tables)")
        self.p, self.n, self.q = p, n, p**n
        self.modulus = list(modulus)
        self.characteristic = p
        self.name = f"GF({p}^{n})"
        self._prime = PrimeField(p)

    def _tup(self, coeffs):
        c = list(coeffs)[: self.n] + [0] * max(0, self.n - len(coeffs))
        return tuple(x % self.p for x in c)

    def gen(self):
        return self._tup([0, 1])

    def elements(self):
        return (self._tup(t) for t in itertools.product(range(self.p), repeat=self.n))

    def coerce(self, x):
        if isinstance(x, tuple) and len(x) == self.n:
            return self._tup(x)
        return self._tup([self._prime.coerce(x)])

    def zero(self):
        return self._tup([])

    def one(self):
        return self._tup([1])

    def add(self, a, b):
        return tuple((x + y) % self.p for x, y in zip(a, b))

    def neg(self, a):
        return tuple((-x) % self.p for x in a)

    def mul(self, a, b):
        return self._tup(_poly_mod(_poly_mul(list(a), list(b), self.p), self.modulus, self.p))

    def pow(self, a, k: int):
        out = self.one()
        base = a
        while k:
            if k & 1:
                out = self.mul(out, base)
            base = self.mul(base, base)
            k >>= 1
        return out

    def inv(self, a):
        if self.is_zero(a):
            raise ZeroDivisionError(f"0 has no inverse in {self.name}")
        return self.pow(a, self.q - 2)

    def is_zero(self, a):
        return all(x % self.p == 0 for x in a)

    def to_str(self, a):
        terms = []
        for i, c in enumerate(a):
            if c:
                terms.append(f"{c}" if i == 0 else (f"x^{i}" if c == 1 else f"{c}*x^{i}"))
        return " + ".join(terms) if terms else "0"


def GF(q, modulus=None):
    """The public finite-field constructor: GF(p) or GF(p^n) (spec §3.2)."""
    if not (isinstance(q, int) and not isinstance(q, bool) and q >= 2):
        raise FieldError(f"GF({q!r}): argument must be a prime power >= 2",
                         hint="examples: GF(2), GF(7), GF(4), GF(27)")
    fac = sympy.factorint(q)
    if len(fac) != 1:
        raise FieldError(f"GF({q}): {q} is not a prime power",
                         hint="finite fields exist only for prime powers")
    (p, n), = fac.items()
    if n == 1:
        if modulus is not None:
            raise FieldError(f"GF({p}) takes no modulus", hint="modulus is for GF(p^n), n >= 2")
        return PrimeField(p)
    return FiniteField(p, n, modulus)
