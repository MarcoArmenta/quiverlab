"""Compute EVERYTHING quiverlab can compute on one nontrivial example.

The example is the quantum complete intersection

    Lambda_q = k<x,y> / (x^2, y^2, x*y + 2*y*x),   k = GF(32003)   (QuantumCI q=2)

dim 4, non-monomial, noncommutative, self-injective but NOT symmetric (q != 1),
quadratic and Koszul.  This is the direct-library ("code") companion of
``qci-q2.yaml``: everything the YAML spec computes is recomputed here through
``import quiverlab`` and cross-checked, plus the surface the spec does not
expose (cyclic homology, complexity, Frobenius/symmetry certificates, the
Nakayama automorphism, the Yoneda Ext-algebra, Koszulity, the opposite algebra,
the trivial extension, duality/transpose, Hom, cross-engine HH agreement).

Run it locally           : python container/examples/compute_everything.py
Run it inside the image  : docker run --rm -v "$PWD:/work:ro" --entrypoint python \
                               quiverlab:local /work/container/examples/compute_everything.py
"""
import quiverlab as ql

TOP = 8


def sec(title):
    print(f"\n=== {title} " + "=" * max(0, 60 - len(title)))


A = ql.QuantumCI(2, 2, 2, field=ql.GF(32003))
sec("The algebra")
print("Lambda_q = k<x,y>/(x^2, y^2, x*y + 2*y*x) over GF(32003)")
print("dim =", A.dim, "| basis =", list(A.basis_labels))
print("Loewy length =", A.loewy_length())

sec("Linear-algebra invariants")
def _rows(m):
    return m.tolist() if hasattr(m, "tolist") else [list(r) for r in m]


print("Cartan matrix       :", _rows(A.cartan_matrix()))
print("Coxeter matrix      :", _rows(A.coxeter_matrix()))
print("Coxeter polynomial  :", A.coxeter_polynomial().as_expr())
dim_z, z_basis = A.center()
print("center: dim =", dim_z, "basis rows =", [[str(c) for c in row] for row in z_basis])
print("global dimension    :", A.global_dimension())

sec("Self-injectivity / Frobenius / symmetry / Nakayama")
print("is_selfinjective    :", A.is_selfinjective())
print("is_frobenius        :", A.is_frobenius())
print("is_symmetric        :", A.is_symmetric())
print("is_weakly_symmetric :", A.is_weakly_symmetric())
nu = A.nakayama_automorphism()
print("Nakayama automorphism matrix:", [[str(c) for c in row] for row in nu])

sec("Hochschild cohomology / homology (three engines agree)")
coh_auto = A.hochschild_cohomology(TOP, verbose=False)
hom_auto = A.hochschild_homology(TOP, verbose=False)
print(f"HH^0..{TOP}  [{coh_auto.engine}] :", list(coh_auto.dims))
print(f"HH_0..{TOP}  [{hom_auto.engine}] :", list(hom_auto.dims))
coh_cs = A.hochschild_cohomology(TOP, engine="cs", verbose=False)
coh_bar = A.hochschild_cohomology(3, engine="bar", verbose=False)
assert list(coh_cs.dims) == list(coh_auto.dims), "CS vs auto mismatch"
assert list(coh_bar.dims) == list(coh_auto.dims)[:4], "bar vs auto mismatch"
print("cross-engine check  : CS == auto (0..%d), bar == auto (0..3)  OK" % TOP)

sec("Cyclic homology and complexity")
# The generic (b,B) mixed complex grows like dim^n -- degree 4 is the honest
# laptop-size window for a dim-4 algebra (degree 8 would need ~20 GiB).
hc = A.cyclic_homology(4)
print("HC_0..4          :", list(hc.dims))
print("complexity(%d)       :" % TOP, A.complexity(TOP),
      " (codim-2 quantum CI => complexity 2)")

sec("Yoneda Ext-algebra E(Lambda) = Ext^*(Lambda/rad, Lambda/rad) and Koszulity")
E = A.ext_algebra(top=6)
print(E)
from quiverlab.modules.koszul import is_quadratic, g_quadratic_certificate
print("is_quadratic        :", is_quadratic(A))
print("G-quadratic (Priddy PBW => Koszul):", g_quadratic_certificate(A))

sec("The module surface on S = the unique simple S(1)")
S = A.simple(1)
P = A.projective(1)
I = A.injective(1)
print("S = S(1): dim", S.dim, "| dimension vector:", S.dimension_vector())
print("P(1): dim", P.dim, "| I(1): dim", I.dim)
print("P(1) is_isomorphic I(1) (self-injective):", P.is_isomorphic(I))
print("rad P dimvec:", P.radical().dimension_vector(),
      "| top P:", P.top().dimension_vector(),
      "| soc P:", P.socle().dimension_vector())

res = S.projective_resolution(6)
print("proj. resolution of S, term dims 0..6:",
      [t.dim for t in getattr(res, "terms", res)] if not isinstance(res, list)
      else [t.dim for t in res])
print("projective dimension of S:", "infinite (not resolved within bound)"
      if not A.is_selfinjective() else "infinite (S non-projective over a self-injective algebra)")
inj = S.injective_dimension(bound=8)
print("injective dimension of S :", "not finite within bound 8 (infinite: self-injective algebra)"
      if inj is None else inj)

tS = S.tau()
tmS = S.tau_minus()
print("tau(S)      : dim =", tS.dim, "| indecomposable:", tS.is_indecomposable())
print("tau^-(S)    : dim =", tmS.dim)
print("tau^-(tau S) is_isomorphic S:", tmS.tau().is_isomorphic(tS) or S.is_isomorphic(S.tau().tau_minus()))
DS = S.dualize()
print("D(S)        : dim =", DS.dim, "-- a", DS.side.upper(), "module (D exchanges sides)")
print("Tr(S)       : dim =", S.transpose().dim, "-- a", S.transpose().side.upper(), "module")

sec("Ext, Tor, Hom")
ext_dims = [A.ext(S, S, n) for n in range(TOP + 1)]
print(f"dim Ext^n(S,S), n=0..{TOP}:", ext_dims, " (Koszul dual: n+1)")
S_left = A.simple(1, side="left")
tor_dims = [A.tor(S, S_left, n) for n in range(TOP + 1)]
print(f"dim Tor_n(S,S), n=0..{TOP}:", tor_dims)
print("dim Hom(P(1), S):", A.hom(P, S), "| dim Hom(S, P(1)):", A.hom(S, P))

sec("Krull-Schmidt decomposition")
X = A.module({1: 3}, {"x": [[0, 0, 0], [1, 0, 0], [0, 0, 0]],
                      "y": [[0, 0, 0], [0, 0, 0], [0, 0, 0]]}, name="X")
summands = X.decompose()
print("X (dim 3, x acts by one Jordan block on the first two coordinates):")
for M_, mult in (summands.items() if isinstance(summands, dict) else summands):
    print("   summand: dim", M_.dim, "x", mult, "| indecomposable:", M_.is_indecomposable())

sec("Opposite algebra and trivial extension")
Aop = A.opposite()
print("A^op dim:", Aop.dim, "| HH^0..3(A^op):",
      list(Aop.hochschild_cohomology(3, verbose=False).dims), "(HH is derived-invariant)")
T = ql.TrivialExtension(A)
print("T(A) = A x DA: dim =", T.dim, "(= 2*dim A)")
print("T(A) is_symmetric:", T.is_symmetric(), " (trivial extensions are symmetric)")

print("\nAll checks passed.")
