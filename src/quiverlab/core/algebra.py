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

    def __repr__(self):
        base = f"Algebra of dimension {self.dim} over {self.domain.name}"
        if self.basis_labels:
            base += "\nbasis: " + ", ".join(self.basis_labels)
        return base
