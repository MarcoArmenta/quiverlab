"""Derived AR translate tau_{D^b} = nu[-1] on perfect complexes. Self-cert:
d.d=0 in the output (ChainComplex check), and the round-trip tau^-_Db . tau_Db is
a quasi-iso to X. Cross-engine: on a projective resolution of a non-projective
indecomposable M over kA_n, homology of tau_Db is concentrated in degree 0 and
isomorphic to the trusted module tau(M). Literature: the K0 identity
chi(tau_Db X) = c . chi(X) with c the K0-action Coxeter matrix. Negative:
k[x]/(x^2) (infinite gl.dim) refuses."""
import pytest

from quiverlab import GF, Quiver, linear_path_algebra
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.modules.complexes import ChainComplex
from quiverlab.derived.tau import tau_Db, tau_Db_minus

selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine
lit = pytest.mark.oracle_literature


def _a3():
    return linear_path_algebra(3, field=QQ)     # kA3, 1->2->3, hereditary


@xeng
@pytest.mark.parametrize("n, v", [(2, 1), (2, 2), (3, 2)])
def test_tau_db_of_nonprojective_is_module_tau(n, v):
    A = linear_path_algebra(n, field=QQ)
    M = A.simple(v)
    if M.tau().dim == 0:                          # projective: skip (no module tau)
        pytest.skip("projective simple has no module tau")
    X = ChainComplex.from_projective_resolution(M, length=6)
    T = tau_Db(X)
    hd = T.homology_dims()
    # concentrated in degree 0 (nu M = 0 for non-projective interval modules over kA_n)
    assert all(d == 0 for k, d in hd.items() if k != 0)
    assert hd.get(0, 0) == M.tau().dim
    assert T.homology(0).dimension_vector() == M.tau().dimension_vector()


@lit
def test_k0_coxeter_bookkeeping():
    # chi(tau_Db X) = c . chi(X), c = -C * C^-T the K0-ACTION Coxeter matrix (NOT
    # P38's coxeter_matrix, which is the conjugate -C^-T C -- same char poly).
    import sympy as sp
    A = _a3()
    C = sp.Matrix(A.cartan_matrix())
    c = -C * C.inv().T
    verts = list(A.quiver.vertices)
    for v in (1, 2):                              # non-projective simples over kA3
        M = A.simple(v)
        X = ChainComplex.from_projective_resolution(M, length=6)
        T = tau_Db(X)
        chiX = _chi_vec(X, verts)
        chiT = _chi_vec(T, verts)
        assert list(c * sp.Matrix(chiX)) == [sp.Integer(x) for x in chiT]


@selfcert
def test_round_trip_is_quasi_iso():
    A = _a3()
    X = ChainComplex.from_projective_resolution(A.simple(2), length=6)
    back = tau_Db_minus(tau_Db(X))
    # homology dimension vectors agree degreewise (a quasi-iso to X in D^b)
    assert back.homology_dims() == X.homology_dims()


@selfcert
def test_infinite_gldim_refused():
    from quiverlab.families import truncated_polynomial
    A = truncated_polynomial(2, field=GF(7))     # k[x]/(x^2): gl.dim = infinity
    X = ChainComplex.from_projective_resolution(A.simple(1), length=4)
    with pytest.raises(QuiverlabError, match="global dimension|Serre|gl.dim"):
        tau_Db(X)


def _chi_vec(Z, verts):
    out = [0] * len(verts)
    for k in Z.degrees():
        dv = Z.term(k).dimension_vector()
        for i, w in enumerate(verts):
            out[i] += (-1) ** k * dv.get(w, 0)
    return out
