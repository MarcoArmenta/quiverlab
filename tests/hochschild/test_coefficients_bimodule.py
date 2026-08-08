"""Exact-Domain Bimodule. check() certifies the A^e-module axioms; the named
constructors reproduce engine/bimodule.py over GF(p) (the ported reference)."""
import pytest
import quiverlab as ql
from quiverlab.hochschild.coefficients import (
    Bimodule, twisted_homology_classes, twisted_cohomology_classes)
from quiverlab.errors import QuiverlabError

pytestmark = pytest.mark.oracle_selfcert


def _A(field):
    return ql.truncated_polynomial(3, field=field)   # k[x]/(x^3), self-injective (symmetric)


def test_regular_reproduces_algebra_multiplication():
    A = _A(ql.CC)
    M = Bimodule.regular(A)
    assert M.dim_M == A.dim
    M.check()                                    # loud QuiverlabError on any axiom failure
    # Lact is literally A.T (row s of Lact[j] = e_j . m_s = A.T[j][s])
    assert M.Lact == [[list(row) for row in A.T[j]] for j in range(A.dim)]


def test_check_rejects_broken_actions():
    # DEVIATION FROM PLAN (documented): the plan's tamper `Ract[1] = Lact[1]` is a
    # NO-OP on the COMMUTATIVE k[x]/(x^3) (there Lact == Ract already), so check()
    # would not fire. We corrupt one entry of the left action instead -- a genuine
    # generic tamper that breaks left-associativity and MUST raise.
    A = _A(ql.CC)
    M = Bimodule.regular(A)
    dom = A.domain
    M.Lact[1][0][0] = dom.add(M.Lact[1][0][0], dom.one())   # corrupt e_1 . m_0
    with pytest.raises(QuiverlabError):
        M.check()


def test_dual_matches_ported_reference_over_gfp():
    # exact-Domain dual == engine/bimodule.py::dual_bimodule (hanlab port) over
    # GF(p). SINGLE-VERTEX and already unit-adapted (k[x]/(x^3)) => A's basis ==
    # the engine basis, so Lact/Ract are directly comparable (across the transpose).
    from quiverlab.engine import bimodule as eb
    from quiverlab.engine.adapter import to_engine
    A = _A(ql.GF(32003))
    eng = to_engine(A.unit_adapted())
    assert Bimodule.dual(A).same_actions_mod_p(eb.dual_bimodule(eng), 32003)


def test_regular_matches_ported_reference_over_gfp():
    from quiverlab.engine import bimodule as eb
    from quiverlab.engine.adapter import to_engine
    A = _A(ql.GF(32003))
    eng = to_engine(A.unit_adapted())
    assert Bimodule.regular(A).same_actions_mod_p(eb.regular_bimodule(eng), 32003)


def test_twisted_nakayama_uses_A_basis_multivertex_gfp():
    # DD1b: a 3-vertex NON-symmetric self-injective algebra over GF(p) where
    # unit_adapted() is a GENUINE change of basis and nu != id. twisted_by_nakayama
    # must build Lact/Ract with nu in A'S OWN basis: check() certifies the
    # A^e-module axioms and FAILS if the engine unit-adapted matrix were used.
    A = ql.NakayamaAlgebra(n=3, l=3, cyclic=True, field=ql.GF(32003))
    assert A.is_symmetric() is False               # nu != identity: nontrivial twist
    Bimodule.twisted_by_nakayama(A).check()        # loud if nu is in the wrong basis
    # CONTROL: the raw ENGINE-basis matrix must NOT validate against the A.T-built
    # actions (A is genuinely not unit-adapted, so the two bases differ).
    nu_engine = A.nakayama_automorphism()          # engine unit-adapted basis
    with pytest.raises(QuiverlabError):
        Bimodule.twisted(A, psi=nu_engine).check()


@pytest.mark.oracle_crossengine
def test_floats_refused_in_from_actions():
    A = ql.Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=ql.CC)
    with pytest.raises(QuiverlabError):
        Bimodule.from_actions(A, 2, left_maps={"a": [[0.5]]}, right_maps={"a": [[0]]})


def test_from_actions_folds_regular():
    # from_actions must FOLD the generator actions to every basis element: feeding
    # the regular bimodule's generator actions back in reconstructs regular exactly.
    A = ql.Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=ql.CC)
    reg = Bimodule.regular(A)
    labels = A.basis_labels
    gens = [i for i, lab in enumerate(labels) if lab.startswith("e_") or lab in A.quiver.arrows]
    left_maps = {labels[i]: [list(r) for r in reg.Lact[i]] for i in gens}
    right_maps = {labels[i]: [list(r) for r in reg.Ract[i]] for i in gens}
    M = Bimodule.from_actions(A, A.dim, left_maps, right_maps)
    assert M.Lact == reg.Lact and M.Ract == reg.Ract


def test_mod_socle_dim():
    A = _A(ql.CC)                                  # k[x]/(x^3): soc_{A^e}A = <x^2>, dim 1
    M = Bimodule.mod_socle(A)
    M.check()
    assert M.dim_M == A.dim - 1


def test_invariants_coinvariants_degree0():
    A = _A(ql.CC)
    reg = Bimodule.regular(A)
    # k[x]/(x^3) is commutative: Z(A)=A, A/[A,A]=A, both dim 3.
    assert reg.invariants_dim() == 3
    assert reg.coinvariants_dim() == 3


def test_twisted_classes_convention():
    # P54 cross-plan contract (chain-level): basis ordering pinned, b_{n-1} b_n = 0,
    # and the class count equals dim HH_n(A, {}_1A_nu) computed from the SAME
    # boundaries. k[x]/(x^3) with {}_1A_nu.
    A = _A(ql.CC)
    M = Bimodule.twisted_by_nakayama(A)
    top = 4
    data = twisted_homology_classes(A, M, top)
    dom = A.domain
    m = A.dim
    from quiverlab.hochschild import bar
    from quiverlab.fields.linalg import rank
    # basis ordering: M-index s outermost, bar multi-index J lexicographic
    for n in range(top + 1):
        assert data[n]["basis"] == bar._cochain_basis(m, n, M.dim_M)
    # b_{n-1} o b_n = 0
    for n in range(2, top + 1):
        bn = data[n]["boundary"]
        bnm1 = data[n - 1]["boundary"]
        prod = [[dom.zero()] * len(bn[0]) for _ in range(len(bnm1))]
        for i in range(len(bnm1)):
            for kk in range(len(bn)):
                a = bnm1[i][kk]
                if dom.is_zero(a):
                    continue
                for j in range(len(bn[0])):
                    prod[i][j] = dom.add(prod[i][j], dom.mul(a, bn[kk][j]))
        assert all(dom.is_zero(x) for row in prod for x in row)
    # class count == dims computed from the boundaries (dim C_n - rk b_n - rk b_{n+1})
    ranks = {}
    for n in range(top + 2):
        bn, _, _ = bar.boundary_matrix(A.unit_adapted(), n, 4_000_000, coefficients=M) \
            if n >= 1 else (None, 0, 0)
        ranks[n] = rank(bn, dom) if bn else 0
    for n in range(top + 1):
        cn = M.dim_M * (m - 1) ** n
        expected = cn - ranks[n] - ranks[n + 1]
        assert len(data[n]["classes"]) == expected


def test_twisted_cohomology_classes_smoke():
    A = _A(ql.CC)
    M = Bimodule.dual(A)
    data = twisted_cohomology_classes(A, M, 3)
    m = A.dim
    from quiverlab.hochschild import bar
    for n in range(4):
        assert data[n]["basis"] == bar._cochain_basis(m, n, M.dim_M)
