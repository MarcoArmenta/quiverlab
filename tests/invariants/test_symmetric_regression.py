"""Plan 29 Part 0 — the ``is_symmetric`` regression battery (QPA-verified live bug).

The former ``is_symmetric`` returned a SILENT WRONG ``False`` on provably
symmetric multi-vertex Nakayama (Brauer star) algebras ``kZ_n/J^L`` with
``n | (L - 1)``: over GF(p) it required the engine's Nakayama automorphism to be
the IDENTITY MATRIX, a sufficient-not-necessary test (``nu`` is defined only up to
inner automorphism). QPA's ``IsSymmetricAlgebra`` returns ``true`` for these.

The self-injective Nakayama ``kZ_n/J^L`` is weakly symmetric iff ``n | (L - 1)``,
and for Nakayama algebras weakly symmetric <=> symmetric (Skowronski–Yamagata,
*Frobenius Algebras I*, EMS 2011). Cross-checked against QPA in
``tests/qpa/test_symmetric_qpa.py``.
"""
import pytest

from quiverlab import CC, GF, Quiver, linear_path_algebra, truncated_polynomial
from quiverlab.families import NakayamaAlgebra, QuantumCI, TrivialExtension
from quiverlab.errors import FieldError
from quiverlab.fields import QQ

_FIELDS = [GF(32003), GF(2), QQ]                      # GF(p^n)/CC covered per-case below


def _star(n, ell, field):
    """The Brauer star kZ_n/J^ell = symmetric Nakayama when n | (ell - 1)."""
    return NakayamaAlgebra(n=n, l=ell, cyclic=True, field=field)


# -- the headline fix: multi-vertex symmetric Nakayama (Brauer stars) ----------
@pytest.mark.parametrize("n,ell", [(2, 3), (3, 4), (4, 5)])
@pytest.mark.parametrize("field", _FIELDS)
def test_brauer_stars_are_symmetric_and_weakly_symmetric(n, ell, field):
    # n | (ell - 1): 2|2, 3|3, 4|4  =>  symmetric (was WRONGLY False over GF(p))
    A = _star(n, ell, field)
    assert A.is_symmetric() is True
    assert A.is_weakly_symmetric() is True
    assert A.is_frobenius() is True


# -- the negatives: n does NOT divide (ell - 1) => not weakly symmetric --------
@pytest.mark.parametrize("n,ell", [(3, 3), (4, 7)])       # 3∤2, 4∤6
@pytest.mark.parametrize("field", _FIELDS)
def test_nakayama_non_dividing_are_not_symmetric(n, ell, field):
    A = _star(n, ell, field)
    assert A.is_weakly_symmetric() is False
    assert A.is_symmetric() is False


def test_self_injective_but_not_symmetric_member():
    # kZ_3/J^3: self-injective (all cyclic Nakayama are) yet NOT symmetric --
    # the Nakayama permutation is the nontrivial 3-cycle (3 ∤ 2).
    A = _star(3, 3, GF(32003))
    assert A.is_selfinjective() is True
    assert A.is_weakly_symmetric() is False
    assert A.is_symmetric() is False


# -- single-vertex positives stay byte-identical -------------------------------
@pytest.mark.parametrize("a", [2, 3, 4])
@pytest.mark.parametrize("field", _FIELDS + [CC])
def test_truncated_polynomial_symmetric_unchanged(a, field):
    # k[x]/(x^a) is commutative Frobenius => symmetric, over every field.
    A = truncated_polynomial(a, field=field)
    assert A.is_symmetric() is True
    assert A.is_weakly_symmetric() is True


# -- weak symmetry is STRICTLY weaker than symmetry ---------------------------
@pytest.mark.parametrize("q", [1, 2])                     # exterior (q=1), quantum CI (q=2)
@pytest.mark.parametrize("field", [QQ, GF(5)])
def test_weakly_symmetric_does_not_imply_symmetric(q, field):
    # Single vertex => identity Nakayama permutation => weakly symmetric, but the
    # Nakayama automorphism nu = diag(1, -q, -1/q, 1) is NOT inner => not symmetric.
    A = QuantumCI(q, field=field)
    assert A.is_weakly_symmetric() is True
    assert A.is_symmetric() is False


def test_commutative_ci_is_symmetric():
    assert QuantumCI(-1, field=QQ).is_symmetric() is True
    assert QuantumCI(-1, field=QQ).is_weakly_symmetric() is True


# -- honest refusal on a presentation-less algebra stays coherent --------------
def test_presentationless_symmetry_refuses_like_frobenius():
    """is_symmetric refuses (loudly) on a quiver-less algebra EXACTLY as
    is_frobenius does -- symmetric would imply Frobenius, so the two surfaces
    must agree. (The Plan-19 refusal contract; see test_refusal_surface.)"""
    from quiverlab.core.algebra import Algebra
    T = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]              # k[x]/(x^2), raw, no quiver
    A = Algebra.from_structure_constants(T, [1, 0], field=CC)
    with pytest.raises(FieldError):
        A.is_frobenius()
    with pytest.raises(FieldError):
        A.is_symmetric()


# -- Trivial extension T(A) = A |x D(A) is ALWAYS symmetric (Happel; ASS2006) --
# ROOT CAUSE (documented): T(A) is built from structure constants with NO quiver
# presentation, so the path-type-basis certifier cannot run and is_symmetric
# refuses loudly (honest -- consistent with is_frobenius also refusing). Turning
# this True requires a double-quiver presentation for TrivialExtension
# (families/trivial_extension.py), owned by a Plan-29 sibling; the trace-form
# certifier here then returns True unchanged (each T(A) already has a
# nondegenerate trace form -- verified during Plan 29 diagnosis). xfail until the
# presentation lands, at which point these flip to xpass.
def _kronecker(field):
    return Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(relations=[], field=field)


def _comm_square(field):
    return Quiver([1, 2, 3, 4],
                  {"a": (1, 2), "b": (1, 3), "c": (2, 4), "d": (3, 4)}
                  ).algebra(relations=["a*c - b*d"], field=field)


@pytest.mark.xfail(raises=FieldError, strict=False,
                   reason="TrivialExtension needs a double-quiver presentation "
                          "(Plan 29 families/ sibling) for the path-type certifier")
@pytest.mark.parametrize("build", [
    lambda: linear_path_algebra(2, field=QQ),            # kA_2
    lambda: linear_path_algebra(3, field=QQ),            # kA_3
    lambda: _kronecker(QQ),                              # 2-Kronecker
    lambda: _comm_square(QQ),                            # commutative square
])
def test_trivial_extension_is_symmetric(build):
    T = TrivialExtension(build())
    assert T.is_frobenius() is True
    assert T.is_symmetric() is True
