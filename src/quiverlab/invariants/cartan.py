"""Cartan and Coxeter data from the quiver presentation (any field).

For a monomial bound quiver algebra with basis the irreducible paths,
C[i][j] = dim e_i A e_j = #(basis paths from vertex i to vertex j) -- an
integer matrix independent of the ground field. The Coxeter matrix is
Phi = -C^{-T} C (requires C invertible over Q, e.g. finite global dimension);
the Coxeter polynomial is charpoly(Phi)."""
import sympy

from quiverlab.errors import QuiverlabError


def cartan_matrix(A):
    Q = A.quiver
    if Q is None or A.basis_labels is None:
        raise QuiverlabError(
            "Cartan matrix needs the quiver presentation",
            hint="construct the algebra via Quiver.algebra(...); "
                 "structure-constant algebras carry no path basis",
        )
    verts = list(Q.vertices)
    vindex = {v: k for k, v in enumerate(verts)}
    n = len(verts)
    C = [[0] * n for _ in range(n)]
    for label in A.basis_labels:
        if label.startswith("e_"):
            # trivial path at vertex v: label 'e_<v>'; recover v by matching reprs
            v = next(w for w in verts if f"e_{w}" == label)
            C[vindex[v]][vindex[v]] += 1
        else:
            word = tuple(label.split("*"))
            C[vindex[Q.word_source(word)]][vindex[Q.word_target(word)]] += 1
    return C


def coxeter_matrix(A):
    C = sympy.Matrix(cartan_matrix(A))
    if C.det() == 0:
        raise QuiverlabError(
            "Cartan matrix is singular: the Coxeter matrix -C^{-T} C is undefined",
            hint="this happens e.g. for infinite global dimension with |det C| != 1",
        )
    Phi = -C.inv().T * C
    if any(x.q != 1 for x in Phi):  # entries should be integers for det = +-1
        # keep exact rationals when det != +-1; return exact sympy entries otherwise
        return [[sympy.nsimplify(Phi[i, j]) for j in range(Phi.cols)] for i in range(Phi.rows)]
    return [[int(Phi[i, j]) for j in range(Phi.cols)] for i in range(Phi.rows)]


def coxeter_polynomial(A):
    """Characteristic polynomial of the Coxeter matrix Phi = -C^{-T} C, as an exact
    sympy Poly in t (no numerical root-finding; the exact spectral_radius / mahler_measure
    invariants live in invariants/spectral.py).

    DOMAIN CAVEAT (det C not in {0, +-1}): when the Cartan is non-unimodular, Phi may have
    rational entries, so the Coxeter polynomial can be over QQ (e.g. t^2 + 3t/2 + 1 for
    C = [[2,1],[0,1]]). It is still EXACT, but such a rational Phi is NOT the classical
    integral Coxeter transformation -- a caveat the coxeter_matrix() sibling surfaces via
    its rational-entry branch. The domain follows the actual COEFFICIENTS, not det C: a
    non-unimodular Cartan may still be integral (k[x]/(x^2) -> t+1 over ZZ; diag(1,2) ->
    (t+1)^2 over ZZ). A singular Cartan (det C == 0) has no Coxeter transformation and
    raises loudly."""
    C = sympy.Matrix(cartan_matrix(A))
    if C.det() == 0:
        raise QuiverlabError(
            "Cartan matrix is singular: no Coxeter polynomial",
            hint="see coxeter_matrix",
        )
    Phi = -C.inv().T * C
    t = sympy.Symbol("t")
    return sympy.Poly(Phi.charpoly(t).as_expr(), t)   # domain inferred from coefficients
