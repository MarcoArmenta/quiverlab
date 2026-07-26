"""Quantum complete intersection A = k<x,y>/(x^a, y^b, xy + q yx) (spec §3.4).
General route; q may be non-rational (needs the Task-4 coefficient grammar).

The two nilpotency degrees a, b default to 2 (the classical dim-4 example); any
a, b >= 2 give the (a*b)-dimensional quantum complete intersection with basis
{ x^i y^j : 0 <= i < a, 0 <= j < b } (Bergh-Erdmann). Plan-33 D3 added a, b as a
BACKWARD-COMPATIBLE generalization: QuantumCI(q) is exactly QuantumCI(q, a=2, b=2)
and every field string 'x^2'/'y^2'/third is byte-identical to the pre-Plan-33 code.

COMMUTATION CONVENTION. The relation is `xy + q yx` (so yx = -xy/q), quiverlab's
established convention across the whole family and its test batteries -- NOT the
textbook `yx - q xy`. The two agree for GENERIC q (both are generic quantum CIs
with the same Hochschild (co)homology; verified: QCI(2,4) coh = [2,2,1,0,...] at a
generic q under either form). They differ only at the special/root-of-unity values
where the convention is load-bearing: under `xy + q yx` the anticommutative exterior
point is q = +1 (ExteriorAlgebra(2) = QuantumCI(1), byte-identical) and the
commutative point k[x,y]/(x^2,y^2) is q = -1. Keeping this convention is what makes
the a=b=2 output byte-identical to before.

The commutation coefficient q is emitted as a BARE grammar token (never a
parenthesised rational): a plain rational ('3', '1/2'), or a sanctioned exact
non-rational scalar ('i', 'E(n)', 'sqrt(k)'). Parenthesising a rational would
smuggle a '-' inside the parens, and _split_terms splits every '+'/'-' -- so
'x*y + (-1)*y*x' is rejected. A negative coefficient therefore folds its sign
OUT of the token, into the relation's own subtraction: 'x*y - 1*y*x'."""
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError


def _q_token(q):
    """Stringify q into an exact grammar token. Strings pass through as the
    non-rational scalar ('i', 'E(3)', 'sqrt(2)'); ints/Fractions stringify to
    exact rationals ('3', '-1', '1/2') via str() (repr() would emit
    'Fraction(1, 2)', which the coefficient grammar does not read)."""
    if isinstance(q, str):
        return q
    return str(q)


def _check_degree(name, d):
    """Nilpotency degree must be a plain int >= 2 (bool excluded: True == 1)."""
    if isinstance(d, bool) or not isinstance(d, int) or d < 2:
        raise QuiverlabError(
            f"QuantumCI: nilpotency degree {name}={d!r} must be an integer >= 2",
            hint="use e.g. QuantumCI(q, a=2, b=4) for k<x,y>/(x^2, y^4, xy + q yx)",
        )


def QuantumCI(q, a=2, b=2, field=None):
    """Quantum complete intersection k<x,y>/(x^a, y^b, xy + q yx).

    q: the commutation coefficient (see the module docstring for the token grammar
        and the `xy + q yx` convention). a, b: the nilpotency degrees of x, y
        (default 2, giving the classical dim-4 example). dim = a*b for any q != 0.
    """
    _check_degree("a", a)
    _check_degree("b", b)
    Q = Quiver([1], {"x": (1, 1), "y": (1, 1)})
    tok = _q_token(q)
    if tok.startswith("-"):                       # fold the sign into the relation's subtraction
        third = f"x*y - {tok[1:]}*y*x"
    else:
        third = f"x*y + {tok}*y*x"
    rels = [f"x^{a}", f"y^{b}", third]            # a=b=2 -> ['x^2', 'y^2', third] (byte-identical)
    A = Q.algebra(relations=rels, field=field)
    A._family_citations = ("quantum_ci", "qci_hh_oracle", "bardzell", "chouhy_solotar")
    return A
