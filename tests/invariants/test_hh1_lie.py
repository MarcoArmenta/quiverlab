"""HH^1 as a Lie algebra = Der(A)/Inn(A) (Plan 70 / R11; RSS 1903.12145; Strametz 2006).
Field-general via the algebra's structure constants; char-0 classification behind a hard gate."""
import pytest

from quiverlab.fields import QQ, GF
from quiverlab.combinat.quiver import Quiver
from quiverlab.families import truncated_polynomial
from quiverlab.errors import QuiverlabError
from quiverlab.invariants.hh1_lie import (
    derivations, inner_derivations, hh1_lie_structure, rss_solvable_certificate,
    _lie_invariants_from_products, _field_general, _jacobi_holds,
    _matmul_commutator_flat, _span_basis, _rank, DEFAULT_MAXDIM,
)

xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _kron(field=QQ):
    return Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(field=field)   # kK2, dim 4


# -------------------------------------------------------------------- Task 1
@xeng
@pytest.mark.parametrize("n,expect", [(2, 1), (3, 2), (4, 3), (5, 4)])   # dim Der(k[x]/x^n)=n-1 (QQ)
def test_dim_hh1_matches_bar_engine_truncpoly(n, expect):
    A = truncated_polynomial(n, field=QQ)
    L = hh1_lie_structure(A)
    assert L.dim == expect
    assert L.dim == A.hochschild_cohomology(top=2).dims[1]               # cross-engine anchor
    assert L.dim_inn == 0                                                # commutative


@xeng
def test_dim_hh1_kronecker():
    A = _kron()
    L = hh1_lie_structure(A)
    assert L.dim == 3 == A.hochschild_cohomology(top=2).dims[1]


@selfcert
def test_inn_is_ideal_of_der():
    """[Der, Inn] subset Inn (so HH^1 = Der/Inn is a Lie algebra) -- CHECKED
    computationally: every commutator [D, ad] of a Der basis element with an Inn basis
    element lies in span(Inn) (adding it does not raise the Inn rank)."""
    A = _kron()
    L = hh1_lie_structure(A)
    assert L.dim_der - L.dim_inn == L.dim
    Der, Inn = derivations(A), inner_derivations(A)
    assert len(Der) == L.dim_der and len(Inn) == L.dim_inn
    dom = A.domain
    inn_basis = _span_basis(Inn, dom)
    inn_rank = _rank(inn_basis, dom)
    assert inn_rank == L.dim_inn > 0                     # Inn nontrivial here (kK2)
    for D in Der:
        for ad in inn_basis:
            br = _matmul_commutator_flat(D, ad, A.dim, dom)
            assert _rank(inn_basis + [br], dom) == inn_rank   # [D, ad] in span(Inn)


@selfcert
@pytest.mark.parametrize("factory", ["kron_qq", "witt_gf3", "x5_qq", "kron_sum_qq"])
def test_jacobi_identity_on_bracket_constants(factory):
    """The HH^1 bracket structure constants satisfy the Jacobi identity on the basis
    (an oracle_selfcert that the Der/Inn representative extraction is a genuine Lie
    algebra) -- verified across the pin zoo incl. the char-p Witt case and sl2 (+) sl2."""
    if factory == "kron_qq":
        A = _kron()
    elif factory == "witt_gf3":
        A = truncated_polynomial(3, field=GF(3))
    elif factory == "x5_qq":
        A = truncated_polynomial(5, field=QQ)
    else:                                                # sl2 (+) sl2, multi-factor
        A = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (1, 2), "c": (3, 4), "d": (3, 4)}).algebra(field=QQ)
    fg = _field_general(A, DEFAULT_MAXDIM)
    assert _jacobi_holds(fg["c"], fg["dim"], A.domain)


@selfcert
def test_der_inn_basis_provenance():
    A = _kron()
    L = hh1_lie_structure(A)
    assert L.basis == "der_inn"
    assert isinstance(L.constants, list) and len(L.constants) == L.dim


# -------------------------------------------------------------------- Task 2
@lit
@pytest.mark.parametrize("n,solv,nilp,ds", [
    (2, True, True,  [1, 0]),
    (3, True, False, [2, 1, 0]),
    (4, True, False, [3, 2, 0]),
    (5, True, False, [4, 3, 1, 0]),
])
def test_truncpoly_char0_solvable(n, solv, nilp, ds):
    A = truncated_polynomial(n, field=QQ)
    L = hh1_lie_structure(A)
    assert L.solvable is solv and L.nilpotent is nilp
    assert L.derived_series_dims == ds


@lit
@pytest.mark.parametrize("p,dim,solvable,perfect", [(2, 2, True, False), (3, 3, False, True),
                                                    (5, 5, False, True), (7, 7, False, True)])
def test_witt_char_p(p, dim, solvable, perfect):
    """n = char = p: HH^1 = W_1 (Jacobson-Witt), simple for p>=3 (perfect, not solvable),
    2-dim solvable for p=2. THE char-p side of the dichotomy."""
    A = truncated_polynomial(p, field=GF(p))
    L = hh1_lie_structure(A)                     # over char p: field-general block only
    assert L.dim == dim and L.solvable is solvable and L.perfect is perfect


@lit
@pytest.mark.parametrize("n,p,solvable", [(4, 2, True), (6, 2, True), (6, 3, False)])
def test_char_divides_n_but_not_prime(n, p, solvable):
    """Card refinement: char|n with n!=p does NOT force the simple Witt case."""
    A = truncated_polynomial(n, field=GF(p))
    assert hh1_lie_structure(A).solvable is solvable


@selfcert
def test_abelian_and_perfect_flags():
    assert hh1_lie_structure(truncated_polynomial(2, field=QQ)).abelian is True   # k[x]/x^2
    assert hh1_lie_structure(_kron(GF(3))).perfect is True                        # sl2 perfect


@selfcert
def test_is_solvable_nilpotent_thin_verdicts():
    from quiverlab.invariants.hh1_lie import is_solvable_hh1, is_nilpotent_hh1
    A = truncated_polynomial(3, field=QQ)
    assert is_solvable_hh1(A) is True and is_nilpotent_hh1(A) is False
    assert _kron().is_solvable_hh1() is False


@selfcert
def test_oversize_budget_refused_loudly():
    A = truncated_polynomial(3, field=QQ)
    with pytest.raises(QuiverlabError):
        hh1_lie_structure(A, budget=2)               # dim 3 > budget 2


# -------------------------------------------------------------------- Task 3
@lit
def test_kronecker_is_sl2_char0():
    A = _kron()                                                        # kK2
    L = hh1_lie_structure(A)
    assert L.dim == 3 and L.radical_dim == 0 and L.semisimple is True and L.simple is True
    assert L.sl2_count == 1 and L.toral_rank == 1 and L.levi_type == "A1"


@lit
def test_trivial_extension_kronecker_sl2_plus_radical():
    from quiverlab.families import TrivialExtension
    A = TrivialExtension(_kron())                                      # T(kK2), dim 8
    L = hh1_lie_structure(A)
    assert L.dim == 4 and L.radical_dim == 1 and L.levi_dim == 3
    assert L.sl2_count == 1 and L.simple is False    # k (x) sl2 : NOT sl2 itself (card refinement)


@lit
@pytest.mark.parametrize("n", [3, 5])
def test_truncpoly_char0_radical_is_whole(n):
    A = truncated_polynomial(n, field=QQ)
    L = hh1_lie_structure(A)
    assert L.radical_dim == L.dim and L.levi_dim == 0 and L.sl2_count == 0


@selfcert
def test_char_p_classification_refused_loudly():
    A = truncated_polynomial(3, field=GF(3))          # W_1 in char 3
    L = hh1_lie_structure(A)
    assert L.dim == 3 and L.solvable is False          # field-general part still computed
    assert L.radical_dim is None and L.sl2_count is None and L.char0_note   # loud gate
    with pytest.raises(QuiverlabError):
        hh1_lie_structure(A, require_char0=True)


@selfcert
def test_base_change_note_gated_on_arithmetic_not_formal_flag():
    L_qq = hh1_lie_structure(_kron(QQ))
    assert L_qq.algebraically_closed is False and L_qq.base_change_note
    # MAJOR fix: the DEFAULT field CC is FORMALLY closed (is_algebraically_closed True)
    # but its arithmetic is exact QQ -- the note MUST still be present.
    L_cc = hh1_lie_structure(Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra())   # default CC
    assert L_cc.algebraically_closed is True                 # formal flag
    assert L_cc.base_change_note and "QQ" in L_cc.base_change_note


@lit
def test_default_field_CC_equals_QQ_run():
    """MAJOR-fix oracle: a request with no field= (domain CC, char 0, formally closed,
    exact QQ arithmetic) COMPUTES and matches the explicit-QQ values."""
    A_cc = truncated_polynomial(3)                            # default CC
    A_qq = truncated_polynomial(3, field=QQ)
    L_cc, L_qq = hh1_lie_structure(A_cc), hh1_lie_structure(A_qq)
    assert L_cc.dim == L_qq.dim == 2 and L_cc.solvable is True
    assert L_cc.derived_series_dims == L_qq.derived_series_dims == [2, 1, 0]
    assert L_cc.radical_dim == L_qq.radical_dim and L_cc.sl2_count == L_qq.sl2_count == 0
    assert L_cc.base_change_note                              # present even though CC formally closed


@lit
def test_gf4_non_prime_field_arithmetic():
    """MINOR (b): exercise GF(p^n), n>1, in the Leibniz/commutator layer (Domain.inv).
    k[x]/x^3 over GF(4) -> dim 2 solvable (char 2 does not divide 3); ql HH^1 = 2."""
    A = truncated_polynomial(3, field=GF(4))
    L = hh1_lie_structure(A)
    assert L.dim == 2 == A.hochschild_cohomology(top=2).dims[1] and L.solvable is True


# -------------------------------------------------------------------- Task 4
from quiverlab.families import RadicalSquareZero


@lit
@pytest.mark.parametrize("Q,field", [
    ({"a": (1, 2), "b": (2, 1)}, QQ),               # k(1<->2)/rad^2 (no loops, no parallel)
])
def test_rss_criterion_positive(Q, field):
    A = RadicalSquareZero(Quiver([1, 2], Q), field=field)
    cert = rss_solvable_certificate(A)
    assert cert["solvable_by_criterion"] is True
    assert hh1_lie_structure(A).solvable is True    # criterion not vacuous: HH^1 != 0 here


@lit
def test_rss_criterion_3cycle_any_char():
    for field in (QQ, GF(3), GF(5)):
        A = RadicalSquareZero(Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (3, 1)}), field=field)
        assert rss_solvable_certificate(A)["solvable_by_criterion"] is True
        assert hh1_lie_structure(A).solvable is True


@selfcert
def test_rss_criterion_kronecker_not_covered():
    A = _kron()                                     # parallel arrows
    cert = rss_solvable_certificate(A)
    assert cert["solvable_by_criterion"] in (False, None) and cert["has_parallel"] is True
    assert hh1_lie_structure(A).solvable is False   # sl2 (correctly not solvable)


@selfcert
def test_rss_criterion_presentation_less_refused():
    from quiverlab.core.algebra import Algebra
    # a structure-constant-only algebra (no Gabriel quiver): loud refusal, but the
    # field-general HH^1-Lie block still computes.
    T = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]      # k[x]/x^2 raw structure constants
    A = Algebra.from_structure_constants(T, unit=[1, 0], field=QQ, check=True)
    assert A.quiver is None
    with pytest.raises(QuiverlabError):
        rss_solvable_certificate(A)
    assert hh1_lie_structure(A).dim == 1          # still computes field-generally


@xeng
@pytest.mark.parametrize("factory", ["truncpoly4_gf5", "kron_gf3"])
def test_gerstenhaber_bracket_agrees(factory):
    """The degree-(1,1) gerstenhaber_brackets route (Plan 35/P51, in-window over GF(p))
    agrees with the Der/Inn commutator on BASIS-INDEPENDENT data (dim, solvable,
    nilpotent, derived-series dims). Structure constants are NOT compared."""
    if factory == "truncpoly4_gf5":
        A = truncated_polynomial(4, field=GF(5))                  # dim HH^1 = 3, solvable
    else:
        A = _kron(GF(3))                                          # sl2
    L = hh1_lie_structure(A)
    g = A.gerstenhaber_brackets(top=2)                            # in-window degree (1,1)
    inv = _lie_invariants_from_products(g)
    assert inv["dim"] == L.dim and inv["solvable"] == L.solvable and inv["nilpotent"] == L.nilpotent
    assert inv["derived_series_dims"] == L.derived_series_dims
