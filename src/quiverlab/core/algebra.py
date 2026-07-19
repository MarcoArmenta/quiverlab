"""The structure-constant Algebra: quiverlab's internal currency (spec §5).
T[i][j] is the coordinate vector of b_i * b_j. 'Unit-adapted' means b_0 = 1_A
(hanlab's convention), which the bar complex requires."""
from quiverlab.errors import QuiverlabError
from quiverlab.fields.linalg import rank, solve


class Algebra:
    def __init__(self, domain, T, unit, basis_labels=None, is_unit_adapted=None, _quiver=None,
                 _relations=None):
        self.domain = domain
        self.T = T
        self.unit = unit
        self.dim = len(T)
        self.basis_labels = basis_labels
        self.quiver = _quiver
        self.relations = _relations
        if is_unit_adapted is None:
            one = domain.one()
            is_unit_adapted = (
                not domain.is_zero(unit[0])
                and domain.eq(unit[0], one)
                and all(domain.is_zero(c) for c in unit[1:])
            )
        self.is_unit_adapted = is_unit_adapted

    # -- arithmetic -----------------------------------------------------------
    def multiply(self, u, v):
        dom = self.domain
        out = [dom.zero()] * self.dim
        for i, ui in enumerate(u):
            if dom.is_zero(ui):
                continue
            for j, vj in enumerate(v):
                if dom.is_zero(vj):
                    continue
                c = dom.mul(ui, vj)
                for t, w in enumerate(self.T[i][j]):
                    if not dom.is_zero(w):
                        out[t] = dom.add(out[t], dom.mul(c, w))
        return out

    # -- construction ---------------------------------------------------------
    @classmethod
    def from_structure_constants(cls, T, unit, field=None, check=True, basis_labels=None):
        if field is None:
            from quiverlab.fields import CC
            field = CC
        m = len(T)
        raw = [x for row in T for vec in row for x in vec] + list(unit)
        parsed = [field.parse_entry(x) for x in raw]
        dom = field.make_domain(parsed)
        Tc = [[[dom.coerce(field.parse_entry(x)) for x in T[i][j]] for j in range(m)]
              for i in range(m)]
        unit_c = [dom.coerce(field.parse_entry(x)) for x in unit]
        A = cls(dom, Tc, unit_c, basis_labels=basis_labels)
        if check:
            A._validate()
        return A

    def _basis_vec(self, i):
        dom = self.domain
        v = [dom.zero()] * self.dim
        v[i] = dom.one()
        return v

    def _validate(self):
        dom = self.domain
        for i in range(self.dim):
            bi = self._basis_vec(i)
            left = self.multiply(self.unit, bi)
            right = self.multiply(bi, self.unit)
            if left != bi or right != bi:
                raise QuiverlabError(
                    f"the given unit vector is not a two-sided unit (fails on basis {i})",
                    hint="check the structure constants and the unit coordinates",
                )
        for i in range(self.dim):
            for j in range(self.dim):
                ij = self.T[i][j]
                for k in range(self.dim):
                    lhs = self.multiply(ij, self._basis_vec(k))
                    rhs = self.multiply(self._basis_vec(i), self.T[j][k])
                    if lhs != rhs:
                        raise QuiverlabError(
                            f"structure constants are not associative: (b{i}·b{j})·b{k} != b{i}·(b{j}·b{k})",
                            hint="re-derive the multiplication table; quiverlab never guesses",
                        )

    # -- base change ----------------------------------------------------------
    def change_of_basis(self, P):
        """New algebra in the basis whose j-th vector has old coordinates column j of P."""
        dom = self.domain
        m = self.dim
        if rank(P, dom) != self.dim:
            raise QuiverlabError("change of basis matrix is singular",
                                 hint="columns must form a basis")
        cols = [[P[i][j] for i in range(m)] for j in range(m)]
        newT = []
        for i in range(m):
            row = []
            for j in range(m):
                prod_old = self.multiply(cols[i], cols[j])
                x = solve(P, prod_old, dom)
                if x is None:
                    raise QuiverlabError("change of basis matrix is singular",
                                         hint="columns must form a basis")
                row.append(x)
            newT.append(row)
        new_unit = solve(P, list(self.unit), dom)
        if new_unit is None:
            raise QuiverlabError("change of basis matrix is singular",
                                 hint="columns must form a basis")
        return Algebra(dom, newT, new_unit, basis_labels=None,
                       _quiver=self.quiver, _relations=self.relations)

    def unit_adapted(self):
        """Return an isomorphic copy whose basis vector 0 is 1_A (spec §5, component 4)."""
        if self.is_unit_adapted:
            return self
        dom = self.domain
        m = self.dim
        j = next(i for i, c in enumerate(self.unit) if not dom.is_zero(c))
        P = [[dom.one() if r == c else dom.zero() for c in range(m)] for r in range(m)]
        for r in range(m):
            P[r][j] = self.unit[r]
        if j != 0:
            for r in range(m):
                P[r][0], P[r][j] = P[r][j], P[r][0]
        out = self.change_of_basis(P)
        labels = None
        if self.basis_labels is not None:
            labels = list(self.basis_labels)
            old0 = labels[j]
            labels[j] = old0 if j == 0 else labels[0]
            labels[0] = "1"
            if j == 0:
                labels[0] = "1"
        out.basis_labels = labels
        out.is_unit_adapted = True
        return out

    def _use_fast_engine(self, engine):
        from quiverlab.fields.primefield import PrimeField
        return engine == "fast" or (
            engine == "auto" and isinstance(self.domain, PrimeField)
        )

    def hochschild_cohomology(self, top, max_cells=4_000_000, engine="auto"):
        """Dimensions of HH^0..HH^top, exact. engine: 'auto' (fast over GF(p),
        bar otherwise), 'bar' (pure, any field), 'fast' (GF(p) only, loud otherwise)."""
        from quiverlab.hochschild.bar import hochschild_cohomology_dims
        from quiverlab.hochschild.table import HHTable

        if engine not in ("auto", "bar", "fast"):
            raise QuiverlabError(f"unknown engine {engine!r}",
                                 hint="choose 'auto', 'bar', or 'fast'")
        if self._use_fast_engine(engine):
            from quiverlab.engine.adapter import engine_cohomology_dims
            dims = engine_cohomology_dims(self, top, max_cells=max_cells)
            return HHTable(dims, "HH^", repr(self).splitlines()[0],
                           engine="hanlab engine (F_p fast rank)")
        return hochschild_cohomology_dims(self, top, max_cells=max_cells)

    def hochschild_homology(self, top, max_cells=4_000_000, engine="auto"):
        """Dimensions of HH_0..HH_top, exact. Same engine semantics as cohomology."""
        from quiverlab.hochschild.bar import hochschild_homology_dims
        from quiverlab.hochschild.table import HHTable

        if engine not in ("auto", "bar", "fast"):
            raise QuiverlabError(f"unknown engine {engine!r}",
                                 hint="choose 'auto', 'bar', or 'fast'")
        if self._use_fast_engine(engine):
            from quiverlab.engine.adapter import engine_homology_dims
            dims = engine_homology_dims(self, top, max_cells=max_cells)
            return HHTable(dims, "HH_", repr(self).splitlines()[0],
                           engine="hanlab engine (F_p fast rank)")
        return hochschild_homology_dims(self, top, max_cells=max_cells)

    # -- modules --------------------------------------------------------------
    def simple(self, v):
        """The simple right module S_v (spec §3.6)."""
        from quiverlab.modules.builders import simple
        return simple(self, v)

    def projective(self, v):
        """The indecomposable projective right module P_v = e_v A (spec §3.6)."""
        from quiverlab.modules.builders import projective
        return projective(self, v)

    def injective(self, v):
        """The indecomposable injective right module I_v = D(A e_v) (spec §3.6)."""
        from quiverlab.modules.builders import injective
        return injective(self, v)

    def hom(self, M, N):
        """dim Hom_A(M, N) for right A-modules M, N (spec §3.6)."""
        from quiverlab.modules.hom import hom_dim
        return hom_dim(M, N)

    def ext(self, M, N, n):
        """dim Ext^n_A(M, N) for right A-modules M, N (spec §3.6)."""
        from quiverlab.modules.ext import ext
        return ext(self, M, N, n)

    def global_dimension(self):
        """Global dimension: exact value or a labeled certified lower bound (spec §3.5)."""
        from quiverlab.modules.ext import global_dimension
        return global_dimension(self)

    def is_selfinjective(self):
        """True iff every indecomposable projective is injective (self-injective =
        Frobenius for a f.d. algebra); exact over any field (spec §3.5)."""
        from quiverlab.modules.ext import is_selfinjective
        return is_selfinjective(self)

    # -- invariants -----------------------------------------------------------
    def cartan_matrix(self):
        """Integer Cartan matrix from the quiver presentation (any field)."""
        from quiverlab.invariants.cartan import cartan_matrix
        return cartan_matrix(self)

    def coxeter_matrix(self):
        """Coxeter matrix -C^{-T} C (exact; loud if the Cartan matrix is singular)."""
        from quiverlab.invariants.cartan import coxeter_matrix
        return coxeter_matrix(self)

    def coxeter_polynomial(self):
        """Characteristic polynomial of the Coxeter matrix, as an exact sympy Poly."""
        from quiverlab.invariants.cartan import coxeter_polynomial
        return coxeter_polynomial(self)

    def _require_prime_field(self, what):
        from quiverlab.errors import FieldError
        from quiverlab.fields.primefield import PrimeField
        if not isinstance(self.domain, PrimeField):
            raise FieldError(
                f"{what} is available over GF(p) today (fast engine); "
                f"this algebra is over {self.domain.name}",
                hint="construct the algebra over a prime field, or wait for the "
                     "later phase that generalizes this invariant",
            )

    def cyclic_homology(self, top):
        """Dimensions of HC_0..HC_top (Connes mixed complex; GF(p) via the engine)."""
        self._require_prime_field("cyclic homology")
        from quiverlab.engine.adapter import to_engine
        from quiverlab.engine.cyclic import cyclic_homology_dims
        from quiverlab.hochschild.table import HHTable
        p = self.domain.p
        out = cyclic_homology_dims(to_engine(self.unit_adapted()), top, primes=(p,))
        dims = [int(d) for d in out[p]]
        return HHTable(dims, "HC_", repr(self).splitlines()[0],
                       engine="hanlab engine (F_p fast rank)")

    def nakayama_automorphism(self):
        """Nakayama automorphism nu as an integer matrix (columns = images) in the
        unit-adapted basis, over GF(p); loud if not Frobenius / not GF(p)."""
        self._require_prime_field("the Nakayama automorphism")
        from quiverlab.engine.adapter import to_engine
        from quiverlab.engine.coxeter import nakayama_automorphism
        S, _ = nakayama_automorphism(to_engine(self.unit_adapted()), self.domain.p)
        return [[int(S[i, j]) for j in range(S.shape[1])] for i in range(S.shape[0])]

    def is_frobenius(self):
        """Is the algebra Frobenius? (GF(p) via the engine.)"""
        self._require_prime_field("the Frobenius test")
        from quiverlab.engine.adapter import to_engine
        from quiverlab.engine.coxeter import is_frobenius
        return bool(is_frobenius(to_engine(self.unit_adapted()), self.domain.p))

    def is_symmetric(self):
        """Is the algebra symmetric? (Frobenius with identity Nakayama automorphism; GF(p).)"""
        self._require_prime_field("the symmetry test")
        if not self.is_frobenius():
            return False
        from quiverlab.engine.adapter import to_engine
        from quiverlab.engine.coxeter import is_identity, nakayama_automorphism
        E = to_engine(self.unit_adapted())
        S, _ = nakayama_automorphism(E, self.domain.p)
        return bool(is_identity(S, self.domain.p))

    def __repr__(self):
        base = f"Algebra of dimension {self.dim} over {self.domain.name}"
        if self.basis_labels:
            base += "\nbasis: " + ", ".join(self.basis_labels)
        return base
