"""Quantum complete intersection A = k<x,y>/(x^2, y^2, xy + q yx) (spec §3.4).
General route; q may be non-rational (needs the Task-4 coefficient grammar).

The commutation coefficient q is emitted as a BARE grammar token (never a
parenthesised rational): a plain rational ('3', '1/2'), or a sanctioned exact
non-rational scalar ('i', 'E(n)', 'sqrt(k)'). Parenthesising a rational would
smuggle a '-' inside the parens, and _split_terms splits every '+'/'-' -- so
'x*y + (-1)*y*x' is rejected. A negative coefficient therefore folds its sign
OUT of the token, into the relation's own subtraction: 'x*y - 1*y*x'."""
from quiverlab.combinat.quiver import Quiver


def _q_token(q):
    """Stringify q into an exact grammar token. Strings pass through as the
    non-rational scalar ('i', 'E(3)', 'sqrt(2)'); ints/Fractions stringify to
    exact rationals ('3', '-1', '1/2') via str() (repr() would emit
    'Fraction(1, 2)', which the coefficient grammar does not read)."""
    if isinstance(q, str):
        return q
    return str(q)


def QuantumCI(q, field=None):
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    tok = _q_token(q)
    if tok.startswith("-"):                       # fold the sign into the relation's subtraction
        third = f"x*y - {tok[1:]}*y*x"
    else:
        third = f"x*y + {tok}*y*x"
    rels = ["x^2", "y^2", third]
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("quantum_ci", "qci_hh_oracle", "bardzell", "chouhy_solotar")
    return A
