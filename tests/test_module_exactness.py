"""The exactness gate at the MODULE surface (v1.0.1).

``Algebra.module`` reads a user-supplied action matrix per arrow. Every other
entry point into the library -- relations, ``QuantumCI``, ``from_structure_constants``
-- refuses an inexact entry loudly. The module surface did not: it coerced only
``int``/``Fraction`` and passed everything else through untouched, so a float was
silently computed with and returned WRONG dimensions (the rank of a float matrix
is not the rank of the exact matrix it approximates).

These tests pin the contract stated in the first line of the project README:
exact only -- floats fail loudly by design.
"""
from fractions import Fraction

import pytest

import quiverlab
from quiverlab import fields
from quiverlab.errors import ExactnessError, QuiverlabError


# --- the algebras under test -------------------------------------------------
# A_2 (1 -> 2), hereditary, one arrow: the smallest algebra whose module actions
# are a single non-trivial block, so an entry-level defect is visible in rad/top/soc.

def _a2(field):
    return quiverlab.linear_path_algebra(2, field=field)


def _action(entries):
    """The full vertex-ordered action of ``a1`` for dims {1:2, 2:2}: the 2x2 block
    ``entries`` placed at rows(vertex 2) x cols(vertex 1) of a 4x4 matrix."""
    mat = [[0] * 4 for _ in range(4)]
    for i in range(2):
        for j in range(2):
            mat[2 + i][j] = entries[i][j]
    return {"a1": mat}


# The rank-1 witness: row 2 is exactly one third of row 1, so the arrow has rank 1
# over QQ. In binary floating point the same matrix is numerically nonsingular.
_EXACT_RANK1 = [[Fraction(3, 5), Fraction(9, 10)],
                [Fraction(1, 5), Fraction(3, 10)]]
_FLOAT_RANK1 = [[0.6, 0.9], [0.2, 0.3]]


def test_float_action_entries_are_refused():
    """The headline: a float in an action matrix must fail loudly, not compute."""
    A = _a2(fields.QQ)
    with pytest.raises(ExactnessError):
        A.module({1: 2, 2: 2}, _action(_FLOAT_RANK1))


def test_float_action_entries_never_reach_a_dimension():
    """Regression for the v1.0.0 defect: the SAME module entered exactly and as
    decimals returned different rad/top/soc. Whatever happens now, the float form
    must not silently produce a dimension."""
    A = _a2(fields.QQ)
    exact = A.module({1: 2, 2: 2}, _action(_EXACT_RANK1))
    assert (exact.radical().dim, exact.top().dim, exact.socle().dim) == (1, 3, 3)

    with pytest.raises(QuiverlabError):
        A.module({1: 2, 2: 2}, _action(_FLOAT_RANK1)).radical()


@pytest.mark.parametrize("field_name", ["QQ", "CC", "GF(5)", "GF(9)"])
def test_float_refused_over_every_field_class(field_name):
    """The gate is a property of the module surface, not of one domain."""
    field = {"QQ": fields.QQ, "CC": quiverlab.CC,
             "GF(5)": quiverlab.GF(5), "GF(9)": quiverlab.GF(9)}[field_name]
    A = _a2(field)
    with pytest.raises(QuiverlabError):
        A.module({1: 1, 2: 1}, {"a1": [[0, 0], [0.5, 0]]})


def test_inexact_entry_types_are_all_refused():
    """float is not the only inexact carrier: numpy floats, sympy Floats, complex
    numbers and decimal strings must all be refused."""
    numpy = pytest.importorskip("numpy")
    import sympy

    A = _a2(fields.QQ)
    for bad in (0.5, numpy.float64(0.5), numpy.float32(0.5),
                sympy.Float(0.5), complex(1, 2), "0.5", "5e-1", True):
        with pytest.raises(QuiverlabError):
            A.module({1: 1, 2: 1}, {"a1": [[0, 0], [bad, 0]]})


def test_exact_fraction_strings_are_accepted():
    """Plan 26 documents ``'1/2'`` string entries for the no-code panel; the public
    ``Algebra.module`` must read them too, not raise a raw TypeError."""
    A = _a2(fields.QQ)
    M = A.module({1: 1, 2: 1}, {"a1": [[0, 0], ["1/2", 0]]})
    assert M.action["a1"][1][0] == Fraction(1, 2)
    assert M.radical().dim == 1


def test_exact_entries_still_construct_unchanged():
    """The gate must not disturb the exact inputs that already worked."""
    A = _a2(fields.QQ)
    M = A.module({1: 2, 2: 2}, _action(_EXACT_RANK1))
    assert M.action["a1"][2][0] == Fraction(3, 5)
    assert M.dim == 4


def test_native_domain_elements_pass_through():
    """An action whose entries are already elements of the algebra's domain -- the
    shape the library itself builds internally -- must still be accepted. Over an
    algebraic extension of CC these are ANP values that no parser can re-read, so
    the gate has to recognize them rather than try to coerce them."""
    ext = quiverlab.CC.make_domain([quiverlab.CC.parse_entry("sqrt(2)"),
                                    quiverlab.CC.parse_entry(1)])
    root2 = ext.coerce("sqrt(2)")
    from quiverlab.modules.module import _coerce_matrix
    assert _coerce_matrix([[root2]], ext)[0][0] == root2


def test_builtin_modules_round_trip_through_the_gate():
    """Every standard construction builds a Module through the same gate."""
    A = _a2(fields.QQ)
    for M in (A.simple(1), A.projective(1), A.injective(2)):
        assert M.dim >= 1


# --- native-element fast path must not become a bypass -----------------------

def test_invalid_gf_pn_tuples_are_refused():
    """`type(x) is type(dom.one())` is `tuple` over GF(p^n), so a naive native
    fast path waves through ANY tuple -- wrong length, out of range, or belonging
    to a different field -- reproducing the silently-wrong-dimension defect this
    module's gate exists to prevent. `dom.coerce` validates tuples, so they must
    never take the fallback."""
    A = _a2(quiverlab.GF(9))
    good = A.domain.coerce((1, 0))
    ok = A.module({1: 1, 2: 1}, {"a1": [[0, 0], [good, 0]]})
    assert ok.dim == 2                                   # the valid element still works
    for bad in ((1, 0, 0), (1,), ()):                    # wrong length: not an element
        with pytest.raises(QuiverlabError):
            A.module({1: 1, 2: 1}, {"a1": [[0, 0], [bad, 0]]})
    # An unreduced-but-valid element is NORMALIZED, not refused: 99 = 0 in GF(3),
    # so (99, 99) is a legitimate way to write (0, 0). The v1.0.0 defect was storing
    # it verbatim; coercing it is the correct behaviour.
    M = A.module({1: 1, 2: 1}, {"a1": [[0, 0], [(99, 99), 0]]})
    assert M.action["a1"][1][0] == (0, 0)


def test_element_of_a_different_finite_field_is_refused():
    """A GF(27) element must not be accepted into a GF(9) module just because both
    are represented as tuples."""
    A9 = _a2(quiverlab.GF(9))
    foreign = quiverlab.GF(27).coerce((0, 1, 0))
    with pytest.raises(QuiverlabError):
        A9.module({1: 1, 2: 1}, {"a1": [[0, 0], [foreign, 0]]})


def test_float_inside_a_tuple_is_refused_by_the_gate():
    A = _a2(quiverlab.GF(9))
    with pytest.raises(QuiverlabError):
        A.module({1: 1, 2: 1}, {"a1": [[0, 0], [(0.5, 0.5), 0]]})
