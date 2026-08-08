"""HH^* as a graded Lie module over HH^1 (Plan 71 / R12; Gerstenhaber deg-1 = Lie derivative).
Field-general on the normalized bar cochain complex; char-0 weights/decomposition behind P70's gate."""
import pytest
from quiverlab.fields import QQ, GF
from quiverlab.combinat.quiver import Quiver
from quiverlab.families import truncated_polynomial
from quiverlab.hochschild.lie_module import lie_module_action, lie_derivative_on_hh

xeng = pytest.mark.oracle_crossengine
selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


def _rad2loops(field):
    return Quiver([1], {"x": (1, 1), "y": (1, 1)}).algebra(
        relations=["x*x", "x*y", "y*x", "y*y"], field=field)


def _kron(field):
    return Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(field=field)


# ---------------------------------------------------------------------------
# Task 1: field-general action rho_n(D), dims, module axiom, inner-acts-zero
# ---------------------------------------------------------------------------
@xeng
@pytest.mark.parametrize("n,dims", [(2, [2, 1, 1]), (3, [3, 2, 2])])   # k[x]/x^n
def test_hh_dims_match_bar_engine(n, dims):
    A = truncated_polynomial(n, field=QQ)
    L = lie_module_action(A, top=2)
    assert L.hh_dims == dims == A.hochschild_cohomology(top=2).dims


@xeng
def test_kronecker_hh_dims():
    A = _kron(QQ)
    L = lie_module_action(A, top=2)
    assert L.hh_dims == [1, 3, 0] and L.hh1_dim == 3


@selfcert
@pytest.mark.parametrize("build", ["trunc3", "rad2loops", "kron"])
def test_module_axiom_and_inner_zero(build):
    A = ({"trunc3": lambda: truncated_polynomial(3, field=QQ),
          "rad2loops": lambda: _rad2loops(QQ),
          "kron": lambda: _kron(QQ)}[build])()
    L = lie_module_action(A, top=2)
    assert L.module_axiom_ok is True and L.inner_acts_zero is True


@selfcert
def test_lie_derivative_primitive_shape():
    """lie_derivative_on_hh returns the d_n x d_n action matrix on HH^n."""
    A = truncated_polynomial(3, field=QQ)
    from quiverlab.invariants.hh1_lie import derivations
    B = A.unit_adapted()
    Der = derivations(B)
    R = lie_derivative_on_hh(A, Der[0], 1)
    assert len(R) == 2 and all(len(row) == 2 for row in R)      # dim HH^1(k[x]/x^3) = 2


@selfcert
def test_basis_provenance_tag():
    A = truncated_polynomial(2, field=QQ)
    L = lie_module_action(A, top=2)
    assert L.basis == "der_inn/bar"


# ---------------------------------------------------------------------------
# Task 2: the in-window Gerstenhaber-bracket sign arbiter (cross-engine)
#
# Comparison is BASIS-INDEPENDENT (the Plan-35 rule): the two engines use
# different HH^1 and HH^n bases, so per-generator char-polys are NOT comparable;
# the conjugation/recombination invariants of the whole representation are (dim
# HH^n, dim of the associative envelope B_n = <rho_n(HH^1), I>, dim of the
# commutant End_{HH^1}(HH^n)). Structure constants are never compared.
# ---------------------------------------------------------------------------
def _envelope_dim(gens, dom, d):
    """dim of the unital associative subalgebra B = <gens, I> of M_d(k)."""
    from quiverlab.fields.linalg import rank
    from quiverlab.modules import linalg_mod as lm

    def flat(M):
        return [M[i][j] for i in range(d) for j in range(d)]
    basis_vecs = [flat(lm.identity(d, dom))]
    all_mats = [lm.identity(d, dom)]
    r = rank(basis_vecs, dom)
    i = 0
    while i < len(all_mats):
        M = all_mats[i]; i += 1
        for g in gens:
            P = lm.matmul(M, g, dom)
            fv = flat(P)
            rr = rank(basis_vecs + [fv], dom)
            if rr > r:
                basis_vecs.append(fv); r = rr; all_mats.append(P)
    return r


def _commutant_dim(gens, dom, d):
    """dim of {X in M_d(k) : X g = g X for every g in gens}."""
    from quiverlab.fields.linalg import nullspace
    rows = []
    for g in gens:
        for i in range(d):
            for l in range(d):
                row = [dom.zero()] * (d * d)
                for k in range(d):
                    row[i * d + k] = dom.add(row[i * d + k], g[k][l])
                    row[k * d + l] = dom.sub(row[k * d + l], g[i][k])
                if any(not dom.is_zero(x) for x in row):
                    rows.append(row)
    return len(nullspace(rows, dom)) if rows else d * d


def _action_invariants(L):
    """{n: (dim HH^n, dim B_n, dim End)} off our HHLieModule (GF(p))."""
    dom = GF(L.characteristic)
    out = {}
    for entry in L.action:
        n, dn = entry["n"], entry["dim"]
        if dn == 0:
            continue
        gens = [[[dom.coerce(int(g[i * dn + j])) for j in range(dn)] for i in range(dn)]
                for g in entry["gens"]]
        out[n] = (dn, _envelope_dim(gens, dom, dn), _commutant_dim(gens, dom, dn))
    return out


def _bracket_action_invariants(g):
    """{n: (dim HH^n, dim B_n, dim End)} off the Plan-35 degree-(1,n) bracket tables:
    the induced map [b_i, -] on HH^n is constants[k][i][j], one matrix per HH^1 gen."""
    import re
    dom = GF(int(re.search(r"GF\((\d+)\)", str(g.basis)).group(1)))
    out = {}
    for (p, n), t in g.tables.items():
        if p != 1:
            continue
        dl, dr, dout = t.dims
        if dr == 0:
            continue
        gens = [[[dom.coerce(int(t.constants[k][i][j])) for j in range(dr)]
                 for k in range(dout)] for i in range(dl)]
        out[n] = (dr, _envelope_dim(gens, dom, dr), _commutant_dim(gens, dom, dr))
    return out


@xeng
@pytest.mark.parametrize("build,p", [("trunc4", 5), ("kron", 3)])
def test_sign_arbiter_in_window(build, p):
    """The field-general L_D route agrees IN-WINDOW over GF(p) with gerstenhaber_brackets
    degree-(1,n) on BASIS-INDEPENDENT data (dim HH^n, dim B_n, dim End). The #1 sign
    risk: p=1 collapses all Koszul signs, so L_D == the shipped (1,n) bracket."""
    A = (truncated_polynomial(4, field=GF(p)) if build == "trunc4"
         else _kron(GF(p)))
    L = lie_module_action(A, top=2)
    inv = _bracket_action_invariants(A.gerstenhaber_brackets(top=2))
    mine = _action_invariants(L)
    assert inv and all(mine[n] == inv[n] for n in inv)


# --- the STRONG entry-wise check (strictly stronger than the induced-map arbiter) ---
def _der_reps_on_engine_basis(A, E):
    """Each Der/Inn rep of HH^1 as a 1-cochain over cochain_basis(E, 1) (int64/GF(p))."""
    import numpy as np
    from quiverlab.engine.scan3 import cochain_basis
    from quiverlab.invariants.hh1_lie import derivations, inner_derivations, _hh1_reps
    B = A.unit_adapted(); dom = B.domain; m = B.dim
    reps = _hh1_reps(derivations(B), inner_derivations(B), dom)
    b1 = cochain_basis(E, 1)
    out = []
    for Dflat in reps:
        v = np.zeros(len(b1), dtype=np.int64)
        for i, (w, j) in enumerate(b1):
            v[i] = int(Dflat[j * m + w[0]]) % dom.p
        out.append(v)
    return out


def _der_matrix_on_engine(E, Dc):
    import numpy as np
    from quiverlab.engine.scan3 import cochain_basis
    m = E.m
    idx1 = {g: i for i, g in enumerate(cochain_basis(E, 1))}
    Dm = np.zeros((m, m), dtype=np.int64)
    for r in E.R:
        for j in range(m):
            Dm[j, r] = int(Dc[idx1[((r,), j)]])
    return Dm


def _L_D_matrix_on_cochains(E, Dc, n):
    """Our L_D (Plan 71 sec 1) as a matrix on the engine's degree-n cochain basis."""
    import numpy as np
    from quiverlab.engine.scan3 import cochain_basis
    m = E.m
    Dm = _der_matrix_on_engine(E, Dc)
    bn = cochain_basis(E, n)
    idx = {g: i for i, g in enumerate(bn)}
    M = np.zeros((len(bn), len(bn)), dtype=np.int64)
    for o, (w, jo) in enumerate(bn):
        for s in range(m):                                  # term 1: D o g
            c = int(Dm[jo, s])
            if c:
                M[o, idx[(w, s)]] += c
        for i in range(n):                                  # term 2: - g o D
            wi = w[i]
            for k in E.R:
                c = int(Dm[k, wi])
                if c:
                    wk = w[:i] + (k,) + w[i + 1:]
                    M[o, idx[(wk, jo)]] -= c
    return M


def _cochain_basis_vectors(E, n):
    import numpy as np
    from quiverlab.engine.scan3 import cochain_basis
    N = len(cochain_basis(E, n))
    for c in range(N):
        e = np.zeros(N, dtype=np.int64); e[c] = 1
        yield e


@xeng
@pytest.mark.parametrize("build,p,n", [("trunc3", 7, 2), ("trunc3", 7, 3), ("kron", 5, 1)])
def test_L_D_equals_circle_cochain_ENTRY_WISE(build, p, n):
    """STRONG form: the L_D operator built on the engine's degree-n COCHAIN basis equals
    the shipped Gerstenhaber bracket cochain-for-cochain over GF(p) -- the p=1 sign
    collapse is exact, not just up to induced-map invariants. (Reaches into
    engine.tt_calculus, as tests may -- the no-reach-into-engine rule is app/GUI code.)"""
    from quiverlab.engine import tt_calculus as TT
    from quiverlab.engine.adapter import to_engine
    A = (truncated_polynomial(3, field=GF(p)) if build == "trunc3" else _kron(GF(p)))
    E = to_engine(A.unit_adapted())
    reps = _der_reps_on_engine_basis(A, E)
    assert reps
    for D in reps:
        M = _L_D_matrix_on_cochains(E, D, n)
        for f in _cochain_basis_vectors(E, n):
            assert (((M @ f) % p) == (TT.gerstenhaber_bracket_cochain(E, 1, n, D, f) % p)).all()


# ---------------------------------------------------------------------------
# Task 3: char-0 weight / torus decomposition (assert the basis-independent PATTERN,
# never the literal integers -- a scalar-doubled torus generator sends {-2,0,2} to
# {-4,0,4}; the reproducible content is the distinct-weight count + PATTERN + per-weight
# dims, the P70 sl2-triple NON-NORMATIVE precedent.)
# ---------------------------------------------------------------------------
from fractions import Fraction
from quiverlab.errors import QuiverlabError


def _weight_entry(L, n):
    return next(e for e in L.weights if e["n"] == n)


def _weights_at(L, n):
    """The DISTINCT scalar weights at degree n (a rank-<=1 torus: kK2, k[x]/x^n) as
    Fractions -- the normalization-free PATTERN input."""
    ws = _weight_entry(L, n)["weights"]
    out = []
    for lam, _dim in ws:
        out.append(Fraction(lam[0]) if lam else Fraction(0))
    return out


def _weight_dims(L, n):
    return [(tuple(lam), d) for lam, d in _weight_entry(L, n)["weights"]]


def _is_sl2_string(ws):                       # symmetric {-c, 0, c}, c != 0, each mult 1
    s = sorted(ws)
    c = -s[0]
    return len(s) == 3 and c != 0 and s == [-c, 0, c]


def _is_equal_gap(ws):                         # arithmetic progression, common difference != 0
    s = sorted(set(ws))
    diffs = {b - a for a, b in zip(s, s[1:])}
    return len(s) < 2 or (len(diffs) == 1 and 0 not in diffs)


@lit
def test_kronecker_hh1_is_sl2_string():
    A = _kron(QQ)                                     # HH^1 = adjoint sl2
    L = lie_module_action(A, top=1)
    assert _is_sl2_string(_weights_at(L, 1))          # PATTERN, normalization-free
    assert set(_weights_at(L, 0)) == {0}              # HH^0 trivial L(0), one weight 0
    assert all(d == 1 for _, d in _weight_dims(L, 1))  # each weight multiplicity 1


@lit
@pytest.mark.parametrize("n", [2, 3])
def test_truncpoly_equal_gap_weights(n):
    """k[x]/x^n: each HH^m has an equal-gap weight progression (the gap is
    normalization-dependent, so assert equal-gap, NOT gap==1)."""
    A = truncated_polynomial(n, field=QQ)
    L = lie_module_action(A, top=2)
    for m in (1, 2):
        assert _is_equal_gap(_weights_at(L, m))       # PATTERN, not a literal gap


@lit
def test_dual_numbers_weight_ladder():
    """k[x]/x^2: the cleanest weight ladder -- one weight per degree, an equal-gap
    ladder across degrees (each HH^{>=1} is a 1-dim weight space)."""
    A = truncated_polynomial(2, field=QQ)
    L = lie_module_action(A, top=4)
    for m in range(1, 5):
        assert len(_weights_at(L, m)) == 1            # a single weight per degree


@selfcert
def test_weights_gated_char_p_but_action_computed():
    """WEIGHTS are char-0-gated (loud); the field-general action is still computed."""
    A = truncated_polynomial(3, field=GF(3))          # W_1 in char 3
    L = lie_module_action(A, top=2)
    assert L.module_axiom_ok is True                  # field-general part still computed
    assert L.weights is None and L.char0_note          # WEIGHTS gated + explained
    with pytest.raises(QuiverlabError):
        lie_module_action(A, top=2, require_char0=True)  # the weight/torus accessor raises


@selfcert
def test_weight_spaces_sum_to_dim():
    A = truncated_polynomial(3, field=QQ)
    L = lie_module_action(A, top=3)
    for n in range(4):
        assert sum(d for _, d in _weight_dims(L, n)) == L.hh_dims[n]


@selfcert
def test_gl2_two_dim_torus():
    """k[x,y]/(x,y)^2: HH^1 = gl2, a 2-dim torus -> weights are 2-tuples; the adjoint
    HH^1 carries the gl2 root pattern (+-1 off-diagonal, 0 with multiplicity)."""
    A = _rad2loops(QQ)
    L = lie_module_action(A, top=2)
    assert _weight_entry(L, 1)["torus_rank"] == 2
    for n in range(3):
        assert sum(d for _, d in _weight_dims(L, n)) == L.hh_dims[n]
