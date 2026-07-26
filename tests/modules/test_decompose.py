"""Krull-Schmidt decomposition of A-modules (Plan 30 Part A, modules/decompose.py).

Theory oracles: simples / indecomposable projectives / injectives are indecomposable;
explicit direct sums decompose to the right summands with multiplicities; radicals match
known socle/top structure; Krull-Schmidt uniqueness (permuted constructions agree up to
iso); tau is additive over a direct sum. The honest budget refusal (char <= dim over
GF(p), where the trace-form radical is unreliable and no split is found) is asserted to
RAISE, not to guess. Field spread CC / GF(2) / GF(7) / GF(32003). QPA is the independent
oracle in tests/qpa/test_decompose_qpa.py.
"""
import pytest

from quiverlab import CC, GF, Quiver, QuantumCI, linear_path_algebra, truncated_polynomial
from quiverlab.errors import QuiverlabError
from quiverlab.modules import linalg_mod as lm
from quiverlab.modules.module import Module

pytestmark = [pytest.mark.oracle_literature, pytest.mark.oracle_selfcert]

# char > dim M for every module tested below -> the trace-form indecomposability
# certificate is rigorous (Dickson / Cohen-Ivanyos-Wales). GF(2) is handled separately.
_BIG_CHAR = [CC, GF(7), GF(32003)]
_ALL = [CC, GF(2), GF(7), GF(32003)]


def _comm_square(field):
    return Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3),
                                 "d": (3, 4)}).algebra(relations=["a*b - c*d"], field=field)


def _two_sink(field):
    # 1 -> 2, 1 -> 3 (no relations): rad P_1 = S_2 (+) S_3, a decomposable radical.
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (1, 3)}).algebra(field=field)


def _direct_sum(mods, name="DS"):
    """The external direct sum (+) mods as a single Module: every algebra basis-label
    acts block-diagonally on the concatenated basis (a genuine module in a non-vertex-
    block basis -- decomposition is basis-order blind)."""
    A = mods[0].algebra
    dom = A.domain
    labels = set().union(*(m.action.keys() for m in mods))
    dims = [m.dim for m in mods]
    N = sum(dims)
    offs, o = [], 0
    for d in dims:
        offs.append(o)
        o += d
    action = {}
    for lab in labels:
        Mb = lm.zeros(N, N, dom)
        for m, off in zip(mods, offs):
            blk = m.action.get(lab)
            if blk is None:
                continue
            for i in range(m.dim):
                for j in range(m.dim):
                    Mb[off + i][off + j] = blk[i][j]
        action[lab] = Mb
    return Module(A, N, action, name=name, side=mods[0].side)


def _reassemble(pairs):
    """(+) M_i^{m_i} rebuilt from a decompose() result, for the 'reassembles M' check."""
    mods = []
    for m, mult in pairs:
        mods.extend([m] * mult)
    return _direct_sum(mods, name="reassembled")


def _dimvec_multiset(pairs):
    return sorted(tuple(sorted(m.dimension_vector().items())) for m, mult in pairs
                  for _ in range(mult))


# --- simples / projectives / injectives are indecomposable (any char) -------
@pytest.mark.parametrize("field", _ALL)
def test_simples_indecomposable(field):
    for A in (linear_path_algebra(2, field=field), linear_path_algebra(3, field=field),
              _comm_square(field)):
        for v in A.quiver.vertices:
            S = A.simple(v)
            assert S.is_indecomposable()
            d = S.decompose()
            assert len(d) == 1 and d[0][1] == 1
            assert d[0][0].is_isomorphic(S)


@pytest.mark.parametrize("field", _ALL)   # P_v, I_v here all have dim End = 1 -> any char
def test_projectives_injectives_indecomposable(field):
    for A in (linear_path_algebra(3, field=field), _comm_square(field)):
        for v in A.quiver.vertices:
            for M in (A.projective(v), A.injective(v)):
                assert M.is_indecomposable(), M
                d = M.decompose()
                assert len(d) == 1 and d[0][1] == 1
                assert d[0][0].is_isomorphic(M)


# --- local-but-dim-End>1 indecomposables need char > dim (trace-form certificate) ---
@pytest.mark.parametrize("field", _BIG_CHAR)
def test_qci_and_truncated_polynomial_indecomposable(field):
    for M in (QuantumCI(-1, field=field).projective(1),
              QuantumCI(-1, field=field).projective(1).radical(),
              truncated_polynomial(3, field=field).projective(1),
              truncated_polynomial(3, field=field).projective(1).radical(),
              truncated_polynomial(2, field=field).projective(1)):
        assert M.is_indecomposable(), M
        assert M.decompose() == M.decompose()          # deterministic
        assert len(M.decompose()) == 1 and M.decompose()[0][1] == 1


# --- explicit direct sums decompose to the right summands + multiplicities --
@pytest.mark.parametrize("field", _ALL)
def test_direct_sum_of_simples(field):
    A = linear_path_algebra(2, field=field)
    S1S1 = _direct_sum([A.simple(1), A.simple(1)], "S1+S1")
    assert S1S1.check_module()[0]
    assert not S1S1.is_indecomposable()
    d = S1S1.decompose()
    assert len(d) == 1                                   # one iso class
    rep, mult = d[0]
    assert mult == 2 and rep.is_isomorphic(A.simple(1))


@pytest.mark.parametrize("field", _BIG_CHAR)
def test_direct_sum_mixed_multiplicities(field):
    # P_1 (+) S_2 (+) P_1 over kA_2: multiplicities (P_1, 2), (S_2, 1).
    A = linear_path_algebra(2, field=field)
    P1, S2 = A.projective(1), A.simple(2)
    M = _direct_sum([P1, S2, P1], "P1+S2+P1")
    assert M.check_module()[0]
    d = M.decompose()
    mult_by = {}
    for rep, mult in d:
        if rep.is_isomorphic(P1):
            mult_by["P1"] = mult
        elif rep.is_isomorphic(S2):
            mult_by["S2"] = mult
        else:
            pytest.fail(f"unexpected summand {rep}")
    assert mult_by == {"P1": 2, "S2": 1}
    assert _reassemble(d).is_isomorphic(M)               # (+) M_i^{m_i} ~ M


# --- radicals vs known structure --------------------------------------------
@pytest.mark.parametrize("field", _BIG_CHAR)
def test_radical_indecomposable_kA3(field):
    # kA_3 (1->2->3): rad P_1 = e_2 A = P_2, indecomposable.
    A = linear_path_algebra(3, field=field)
    radP1 = A.projective(1).radical()
    assert radP1.is_indecomposable()
    d = radP1.decompose()
    assert len(d) == 1 and d[0][1] == 1
    assert d[0][0].is_isomorphic(A.projective(2))


@pytest.mark.parametrize("field", _ALL)
def test_radical_decomposes_two_sink(field):
    # 1->2, 1->3 (no relations): rad P_1 = S_2 (+) S_3.
    A = _two_sink(field)
    radP1 = A.projective(1).radical()
    assert not radP1.is_indecomposable()
    d = radP1.decompose()
    reps = [rep for rep, _ in d]
    assert all(mult == 1 for _, mult in d)
    assert _dimvec_multiset(d) == sorted(
        [tuple(sorted(A.simple(2).dimension_vector().items())),
         tuple(sorted(A.simple(3).dimension_vector().items()))])
    assert any(r.is_isomorphic(A.simple(2)) for r in reps)
    assert any(r.is_isomorphic(A.simple(3)) for r in reps)


# --- Krull-Schmidt uniqueness (permuted / regrouped constructions agree) -----
@pytest.mark.parametrize("field", _BIG_CHAR)
def test_krull_schmidt_uniqueness(field):
    A = linear_path_algebra(2, field=field)
    P1, S2 = A.projective(1), A.simple(2)
    M1 = _direct_sum([P1, S2, P1], "order-a")
    M2 = _direct_sum([S2, P1, P1], "order-b")            # different ordering
    M3 = _direct_sum([P1, P1, S2], "order-c")            # yet another
    d1, d2, d3 = M1.decompose(), M2.decompose(), M3.decompose()
    assert _dimvec_multiset(d1) == _dimvec_multiset(d2) == _dimvec_multiset(d3)
    # multiplicities match up to iso, not merely dim vectors
    def _mult_of(d, X):
        return sum(mult for rep, mult in d if rep.is_isomorphic(X))
    for d in (d1, d2, d3):
        assert _mult_of(d, P1) == 2 and _mult_of(d, S2) == 1


# --- CC number-field factoring path (SympyExactDomain, sdom = QQ<i>) ---------
def test_decompose_over_number_field():
    # QuantumCI(i) computes in QQ<i>; the split search factors over that number field.
    A = QuantumCI("i", field=CC)
    S = A.simple(1)
    SS = _direct_sum([S, S], "S+S")
    assert not SS.is_indecomposable()
    d = SS.decompose()
    assert len(d) == 1 and d[0][1] == 2 and d[0][0].is_isomorphic(S)


# --- the honest budget refusal (char <= dim over GF(p), no split found) ------
def test_loud_budget_refusal_gf2():
    # k[x]/(x^2) as a right module over GF(2): End = k[x]/(x^2) is LOCAL (indecomposable)
    # but char 2 <= dim 2 makes the trace-form radical unreliable and there is no Fitting
    # split -> the engine must REFUSE loudly, never silently guess.
    Pg2 = truncated_polynomial(2, field=GF(2)).projective(1)
    with pytest.raises(QuiverlabError, match="characteristic 2"):
        Pg2.is_indecomposable()
    with pytest.raises(QuiverlabError, match="characteristic 2"):
        Pg2.decompose()
    # QuantumCI over GF(2) (dim 4) is the same honest-refusal situation.
    Q = QuantumCI(-1, field=GF(2)).projective(1)
    with pytest.raises(QuiverlabError):
        Q.decompose()


def test_gf2_dim_end_one_still_certifies():
    # Over GF(2) the dim-End = 1 fast path is rigorous in every characteristic, so
    # ordinary indecomposable projectives still decompose without refusing.
    A = linear_path_algebra(3, field=GF(2))
    for v in A.quiver.vertices:
        assert A.projective(v).is_indecomposable()


# --- tau-additivity smoke: tau(X (+) Y) ~ tau X (+) tau Y --------------------
@pytest.mark.parametrize("field", [CC, GF(7)])
def test_tau_additivity_smoke(field):
    A = linear_path_algebra(3, field=field)             # 1->2->3; S_1, S_2 non-projective
    S1, S2 = A.simple(1), A.simple(2)
    M = _direct_sum([S1, S2], "S1+S2")
    lhs = M.tau()                                       # tau of the direct sum
    rhs = _direct_sum([S1.tau(), S2.tau()], "tauS1+tauS2")
    assert lhs.is_isomorphic(rhs)


# --- the zero module is not indecomposable, and decomposes to nothing --------
def test_zero_module():
    A = linear_path_algebra(3)                          # kA_3 over CC
    z = A.projective(3).radical()                       # rad P_3 = rad S_3 = 0
    assert z.dim == 0
    assert not z.is_indecomposable()
    assert z.decompose() == []
