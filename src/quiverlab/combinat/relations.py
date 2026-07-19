"""Relation strings -> exact linear combinations of parallel paths (spec §3.3).
Grammar (this phase): terms joined by + and -; term = [coeff*] arrows;
'p^k' repeats an arrow k times. Paths read LEFT TO RIGHT.

A coefficient factor is either a plain rational ('2', '-1/3') or a sanctioned
exact non-rational scalar token -- the imaginary unit 'i', a root of unity
'E(n)', a radical 'sqrt(k)', or a parenthesised rational '(2)' -- parsed exactly
via sympy (Plan 06 Task 4, spec §3.3). Rational factors continue to yield
``Fraction`` (Plan-03/04 byte-compatible); non-rational scalars yield a sympy
``Expr``. The parser stays field-agnostic: it never picks a field. Floats fail
loudly (``ExactnessError``)."""
import re
from dataclasses import dataclass
from fractions import Fraction

import sympy

from quiverlab.errors import ExactnessError, RelationError
from quiverlab.fields.domain import parse_rational

_COEFF = re.compile(r"^[+-]?\d+(/\d+)?$")
# Arrow names start with a letter or underscore; a factor beginning with a
# digit, dot, or sign is numeric-intent and must go through the exact parser
# (so float-style input fails loudly instead of being mistaken for an arrow).
_NUMERIC_INTENT = re.compile(r"^[+-]?[.\d]")
_POW = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)\^(\d+)$")


def _E(n):
    """Primitive n-th root of unity exp(2*pi*i/n), exact (GAP's E(n) convention)."""
    return sympy.exp(2 * sympy.pi * sympy.I / int(n))


def _exact_scalar(tok):
    """Parse an exact non-rational coefficient token (i, E(n), sqrt(k), rationals).
    Returns Fraction for rationals, a sympy exact Expr otherwise. Loud on floats.
    Returns None when the token is not an exact scalar at all (e.g. an arrow name
    or an unknown symbol), leaving the caller to raise the usual RelationError."""
    try:
        expr = sympy.sympify(tok, locals={"i": sympy.I, "E": _E},
                             rational=False, evaluate=True)
    except (sympy.SympifyError, TypeError, SyntaxError, ValueError):
        return None
    if not getattr(expr, "is_number", False):
        return None
    if expr.atoms(sympy.Float):
        raise ExactnessError(f"coefficient {tok!r} contains a floating-point number",
                             hint="write exact scalars: '1/2', 'i', 'sqrt(2)', 'E(3)'")
    if expr.is_rational:
        return Fraction(int(sympy.numer(expr)), int(sympy.denom(expr)))
    return expr


@dataclass(frozen=True)
class Relation:
    terms: tuple  # tuple[tuple[Fraction | sympy.Expr, tuple[str, ...]], ...]
    source: object
    target: object

    @property
    def is_monomial(self) -> bool:
        return len(self.terms) == 1

    @property
    def max_length(self) -> int:
        return max(len(w) for _, w in self.terms)

    @property
    def min_length(self) -> int:
        return min(len(w) for _, w in self.terms)

    def __repr__(self):
        bits = []
        for c, w in self.terms:
            path = "*".join(w)
            if c == 1:
                bits.append(path)
            elif c == -1:
                bits.append(f"-{path}")
            else:
                bits.append(f"{c}*{path}")
        return " + ".join(bits).replace("+ -", "- ")


def _split_terms(s: str):
    """'a*b - 2*c + d' -> ['a*b', '-2*c', '+d'] (sign attached)."""
    s = s.replace(" ", "")
    if not s:
        raise RelationError("empty relation string", hint="e.g. 'a*b - c' or 'x^2'")
    out, cur = [], ""
    for ch in s:
        if ch in "+-" and cur:
            out.append(cur)
            cur = ch
        else:
            cur += ch
    out.append(cur)
    return out


def _parse_term(term: str, quiver):
    sign = Fraction(1)
    if term and term[0] in "+-":
        sign = Fraction(-1) if term[0] == "-" else Fraction(1)
        term = term[1:]
    factors = term.split("*")
    coeff = sign
    word: list[str] = []
    for f in factors:
        if not f:
            raise RelationError(f"malformed term {term!r}", hint="factors are joined by single *")
        if _COEFF.match(f):
            if word:
                raise RelationError(f"coefficient {f!r} appears after arrows in {term!r}",
                                    hint="write coefficients first: '2*a*b'")
            coeff = coeff * parse_rational(f)
            continue
        if _NUMERIC_INTENT.match(f):
            # Numeric-intent but not a plain rational: let the exact parser
            # decide. '0.5' raises ExactnessError (loud); an exotic exact form
            # returns a Fraction, handled just like the _COEFF branch above.
            if word:
                raise RelationError(f"coefficient {f!r} appears after arrows in {term!r}",
                                    hint="write coefficients first: '2*a*b'")
            coeff = coeff * parse_rational(f)
            continue
        m = _POW.match(f)
        if m:
            name, k = m.group(1), int(m.group(2))
            reps = [name] * k
        else:
            name, reps = f, [f]
        if name not in quiver.arrows:
            # Not an arrow: a factor may instead be a sanctioned exact non-rational
            # scalar (i, E(n), sqrt(k), a parenthesised rational). Try that before
            # failing; only a token that is neither an arrow nor an exact scalar is
            # a genuine RelationError (a float token still fails loudly, inside
            # _exact_scalar, as an ExactnessError).
            scalar = _exact_scalar(f)
            if scalar is not None:
                if word:
                    raise RelationError(
                        f"coefficient {f!r} appears after arrows in {term!r}",
                        hint="write coefficients first: '2*a*b'")
                coeff = coeff * scalar
                continue
            raise RelationError(f"unknown arrow {name!r} in relation term {term!r}",
                                hint=f"arrows are {sorted(quiver.arrows)}")
        word.extend(reps)
    if not word:
        raise RelationError(f"term {term!r} has no arrows",
                            hint="relations live in the arrow ideal")
    w = tuple(word)
    if not quiver.compose_ok(w):
        bad = next((i for i, (a, b) in enumerate(zip(w, w[1:]))
                    if quiver.target(a) != quiver.source(b)))
        raise RelationError(
            f"path {'*'.join(w)} is not composable: target({w[bad]}) = "
            f"{quiver.target(w[bad])} but source({w[bad + 1]}) = {quiver.source(w[bad + 1])}",
            hint="paths compose left to right: a*b needs target(a) == source(b)",
        )
    return coeff, w


def parse_relation(s: str, quiver) -> Relation:
    if not isinstance(s, str):
        raise RelationError(f"relations are strings, got {type(s).__name__}",
                            hint="e.g. 'a*b - c'")
    parsed = [_parse_term(t, quiver) for t in _split_terms(s)]
    # Combine like terms (sum coefficients per word, first-occurrence order),
    # then drop zero-sum words so full cancellation is caught as zero below.
    # Summing starts from Fraction(0): all-rational words stay Fraction (byte
    # compatible); any sympy Expr summand promotes that word's coefficient.
    combined: dict[tuple, "Fraction | sympy.Expr"] = {}
    for c, w in parsed:
        combined[w] = combined.get(w, Fraction(0)) + c
    terms = [(c, w) for w, c in combined.items() if c != 0]
    if not terms:
        raise RelationError(f"relation {s!r} is identically zero", hint="remove it")
    srcs = {quiver.word_source(w) for _, w in terms}
    tgts = {quiver.word_target(w) for _, w in terms}
    if len(srcs) > 1 or len(tgts) > 1:
        raise RelationError(
            f"terms of {s!r} are not parallel (sources {sorted(srcs, key=str)}, "
            f"targets {sorted(tgts, key=str)})",
            hint="every summand of a relation must share source and target",
        )
    return Relation(tuple(terms), srcs.pop(), tgts.pop())


def parse_relations(strings, quiver) -> list:
    return [parse_relation(s, quiver) for s in strings]
