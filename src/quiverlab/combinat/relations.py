"""Relation strings -> exact linear combinations of parallel paths (spec §3.3).
Grammar (this phase): terms joined by + and -; term = [rational*] arrows;
'p^k' repeats an arrow k times. Paths read LEFT TO RIGHT."""
import re
from dataclasses import dataclass
from fractions import Fraction

from quiverlab.errors import RelationError
from quiverlab.fields.domain import parse_rational

_COEFF = re.compile(r"^[+-]?\d+(/\d+)?$")
# Arrow names start with a letter or underscore; a factor beginning with a
# digit, dot, or sign is numeric-intent and must go through the exact parser
# (so float-style input fails loudly instead of being mistaken for an arrow).
_NUMERIC_INTENT = re.compile(r"^[+-]?[.\d]")
_POW = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)\^(\d+)$")


@dataclass(frozen=True)
class Relation:
    terms: tuple  # tuple[tuple[Fraction, tuple[str, ...]], ...]
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
    combined: dict[tuple, Fraction] = {}
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
