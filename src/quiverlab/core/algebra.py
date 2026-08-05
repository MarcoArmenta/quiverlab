"""The structure-constant Algebra: quiverlab's internal currency (spec §5).
`T[i][j]` is the coordinate vector of b_i * b_j. 'Unit-adapted' means b_0 = 1_A
(hanlab's convention), which the bar complex requires."""
from quiverlab.errors import DepthLimitError, QuiverlabError
from quiverlab.fields.linalg import rank, solve


class Algebra:
    def __init__(self, domain, T, unit, basis_labels=None, is_unit_adapted=None, _quiver=None,
                 _relations=None, _family_citations=()):
        self.domain = domain
        self.T = T
        self.unit = unit
        self.dim = len(T)
        self.basis_labels = basis_labels
        self.quiver = _quiver
        self.relations = _relations
        self._family_citations = tuple(_family_citations)
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

    # -- citations ------------------------------------------------------------
    def _engine_citations(self):
        # HH^*/HH_* dimensions are produced by the (normalized/fast-rank) bar complex,
        # whatever the presentation -- so the engine key is always "bar" (Hochschild1945).
        # Resolution-specific keys (bardzell / chouhy_solotar) are carried by the family
        # stamp where the family author declares them relevant (e.g. QuantumCI), and by
        # the Plan-04/05 resolution/ops result objects that actually run those engines.
        return ("bar",)

    def citations(self):
        """Registry keys relevant to this algebra: its family stamp plus the HH engine
        (spec §3.9). Every key resolves via quiverlab.citations."""
        seen, out = set(), []
        for k in tuple(getattr(self, "_family_citations", ())) + self._engine_citations():
            if k not in seen:
                seen.add(k)
                out.append(k)
        return tuple(out)

    def _auto_cs_routes(self):
        """True iff this is a NON-MONOMIAL admissible presentation — the ONLY case the
        opt-in auto_cs mode hands to the Chouhy-Solotar engine. Monomial presentations
        (Bardzell already covers them via the fast engine) and anything whose reduction
        system cannot be certified admissible fall through to the frozen auto dispatch,
        so enabling auto_cs never changes a result the shipped 'auto' already produces."""
        if self.quiver is None or self.relations is None:
            return False
        try:
            from quiverlab.resolutions_cs.build import reduction_system_of
            rs = reduction_system_of(self)
        except QuiverlabError:
            return False
        if not rs.is_confluent or not rs.irreducibles:
            return False
        return not all(len(r.tail) == 0 for r in rs.rules)

    def _route_to_cs(self, engine, auto_cs):
        """Explicit engine='cs' always routes to CS; engine='auto' routes ONLY when the
        caller opts in with auto_cs=True AND the algebra is non-monomial admissible.
        Default engine='auto' (auto_cs=False) reaches CS in exactly one further way:
        the depth FALLBACK (see _cs_depth_fallback) at the point where bar/fast would
        raise DepthLimitError — the Marco-2026-07-26 amendment to Pillar-4 (recorded
        in the dispatch trace, in-window results byte-unchanged)."""
        if engine == "cs":
            return True
        return engine == "auto" and auto_cs and self._auto_cs_routes()


    def _cs_depth_fallback(self, side, rec, top, max_cells, cause):
        """The Marco-2026-07-26 dispatch amendment: engine='auto' no longer DIES at
        the bar/fast depth wall when the algebra carries a quiver presentation -- it
        reroutes to the Chouhy-Solotar engine (recorded in the dispatch trace, never
        silent; every in-window result is byte-unchanged because the fallback fires
        only where auto previously raised DepthLimitError). Presentation-less
        algebras still refuse honestly (CS needs the presentation)."""
        from quiverlab.trace.events import Dispatch
        if rec is not None:
            rec.append(Dispatch(
                route="chouhy-solotar",
                reason="the %s exceeded max_cells at this depth; auto reroutes to the "
                       "Chouhy-Solotar resolution over the admissible presentation "
                       "(exact, certified per instance)" % cause,
                n_relations=len(self.relations or ())))
        from quiverlab.resolutions_cs.homology import (cs_cohomology_dims,
                                                       cs_homology_dims)
        fn = cs_cohomology_dims if side == "coh" else cs_homology_dims
        return fn(self, top, max_cells=max_cells, trace=rec)

    def hochschild_cohomology(self, top, max_cells=4_000_000, engine="auto",
                              auto_cs=False, verbose=None, trace=None):
        """Dimensions of HH^0..HH^top, exact. engine: 'auto' (fast over GF(p),
        bar otherwise), 'bar' (pure, any field), 'fast' (GF(p) only, loud otherwise),
        'cs' (Chouhy-Solotar, any admissible presentation over any field). Set
        auto_cs=True to let engine='auto' route non-monomial admissible algebras to CS
        up front. Default 'auto' additionally FALLS BACK to CS at the exact depth
        where the bar/fast route would raise DepthLimitError, for quiver-presented
        algebras (recorded in the dispatch trace; in-window results byte-unchanged;
        presentation-less algebras still refuse). verbose: per-call override
        of quiverlab.verbose (None defers to it); when on, a recorder captures the
        worked steps. trace: an explicit event sink (list or Trace) the engines fill
        with resolution/rank/dispatch events; passing it is programmatic and does not
        by itself request a file."""
        import quiverlab
        from quiverlab.hochschild.bar import hochschild_cohomology_dims
        from quiverlab.hochschild.table import HHTable
        from quiverlab.trace.events import Dispatch
        from quiverlab.trace.recorder import Trace, resolve_verbose

        if engine not in ("auto", "bar", "fast", "cs"):
            raise QuiverlabError(f"unknown engine {engine!r}",
                                 hint="choose 'auto', 'bar', 'fast', or 'cs'")
        want = resolve_verbose(verbose, quiverlab.verbose)
        rec = trace if trace is not None else (Trace() if want else None)
        if self._route_to_cs(engine, auto_cs):
            from quiverlab.resolutions_cs.homology import cs_cohomology_dims
            if rec is not None:
                rec.append(Dispatch(
                    route="chouhy-solotar",
                    reason="general Chouhy-Solotar resolution over the admissible presentation",
                    n_relations=len(self.relations or ())))
            table = cs_cohomology_dims(self, top, max_cells=max_cells, trace=rec)  # CS fills rec
        elif self._use_fast_engine(engine):
            from quiverlab.engine.adapter import engine_cohomology_dims
            if rec is not None:
                rec.append(Dispatch(
                    route="hanlab fast GF(p) rank",
                    reason="domain is a prime field; the exact mod-p rank engine applies",
                    n_relations=len(self.relations or ())))
            try:
                dims = engine_cohomology_dims(self, top, max_cells=max_cells)  # plain list[int]
                table = HHTable(dims, "HH^", repr(self).splitlines()[0],
                                engine="hanlab engine (F_p fast rank)")     # WRAP the list
            except DepthLimitError:
                if engine != "auto" or self.quiver is None:
                    raise
                table = self._cs_depth_fallback("coh", rec, top, max_cells,
                                                "fast GF(p) bar basis")
        else:
            if rec is not None:
                rec.append(Dispatch(
                    route="normalized bar complex",
                    reason="domain is not a prime field; the exact bar oracle is used",
                    n_relations=len(self.relations or ())))
            try:
                table = hochschild_cohomology_dims(self, top, max_cells=max_cells,
                                                   trace=rec)
            except DepthLimitError:
                if engine != "auto" or self.quiver is None:
                    raise
                table = self._cs_depth_fallback("coh", rec, top, max_cells,
                                                "bar oracle")
        table.references = self.citations()   # FROZEN contract (family+engine keys); Task 11 must NOT change it
        if want and trace is None and rec is not None:
            from quiverlab.trace.provenance import references_for, resolve_references
            from quiverlab.trace.writer import write_trace
            # References SECTION = engine keys implied by the trace's Dispatch, resolved
            # through bibliography(); table.references stays self.citations() (untouched).
            write_trace(list(rec), table, algebra=self, kind="HH^", top=top,
                        references=resolve_references(references_for(rec)))
        return table

    def hochschild_homology(self, top, max_cells=4_000_000, engine="auto",
                            auto_cs=False, verbose=None, trace=None):
        """Dimensions of HH_0..HH_top, exact. Same engine semantics as cohomology
        (including 'cs', auto_cs, verbose, and the trace event sink)."""
        import quiverlab
        from quiverlab.hochschild.bar import hochschild_homology_dims
        from quiverlab.hochschild.table import HHTable
        from quiverlab.trace.events import Dispatch
        from quiverlab.trace.recorder import Trace, resolve_verbose

        if engine not in ("auto", "bar", "fast", "cs"):
            raise QuiverlabError(f"unknown engine {engine!r}",
                                 hint="choose 'auto', 'bar', 'fast', or 'cs'")
        want = resolve_verbose(verbose, quiverlab.verbose)
        rec = trace if trace is not None else (Trace() if want else None)
        if self._route_to_cs(engine, auto_cs):
            from quiverlab.resolutions_cs.homology import cs_homology_dims
            if rec is not None:
                rec.append(Dispatch(
                    route="chouhy-solotar",
                    reason="general Chouhy-Solotar resolution over the admissible presentation",
                    n_relations=len(self.relations or ())))
            table = cs_homology_dims(self, top, max_cells=max_cells, trace=rec)  # CS fills rec
        elif self._use_fast_engine(engine):
            from quiverlab.engine.adapter import engine_homology_dims
            if rec is not None:
                rec.append(Dispatch(
                    route="hanlab fast GF(p) rank",
                    reason="domain is a prime field; the exact mod-p rank engine applies",
                    n_relations=len(self.relations or ())))
            try:
                dims = engine_homology_dims(self, top, max_cells=max_cells)  # plain list[int]
                table = HHTable(dims, "HH_", repr(self).splitlines()[0],
                                engine="hanlab engine (F_p fast rank)")   # WRAP the list
            except DepthLimitError:
                if engine != "auto" or self.quiver is None:
                    raise
                table = self._cs_depth_fallback("hom", rec, top, max_cells,
                                                "fast GF(p) bar basis")
        else:
            if rec is not None:
                rec.append(Dispatch(
                    route="normalized bar complex",
                    reason="domain is not a prime field; the exact bar oracle is used",
                    n_relations=len(self.relations or ())))
            try:
                table = hochschild_homology_dims(self, top, max_cells=max_cells,
                                                 trace=rec)
            except DepthLimitError:
                if engine != "auto" or self.quiver is None:
                    raise
                table = self._cs_depth_fallback("hom", rec, top, max_cells,
                                                "bar oracle")
        table.references = self.citations()   # FROZEN contract (family+engine keys); Task 11 must NOT change it
        if want and trace is None and rec is not None:
            from quiverlab.trace.provenance import references_for, resolve_references
            from quiverlab.trace.writer import write_trace
            # References SECTION = engine keys implied by the trace's Dispatch, resolved
            # through bibliography(); table.references stays self.citations() (untouched).
            write_trace(list(rec), table, algebra=self, kind="HH_", top=top,
                        references=resolve_references(references_for(rec)))
        return table

    # -- modules --------------------------------------------------------------
    def _sided_builder(self, name, v, side):
        """Build S/P/I on the requested side (Plan 24). A left A-module IS a right
        A^op-module: route left construction through A^op and re-tag "left"."""
        if side == "right":
            from quiverlab.modules import builders
            return getattr(builders, name)(self, v)
        if side == "left":
            from quiverlab.modules.opposite import opposite_algebra
            return getattr(opposite_algebra(self), name)(v).with_side("left")
        from quiverlab.errors import QuiverlabError
        raise QuiverlabError(f'side must be "right" or "left", got {side!r}')

    def simple(self, v, side="right"):
        """The simple S_v (spec §3.6). ``side="left"`` for the left A-module (Plan 24)."""
        return self._sided_builder("simple", v, side)

    def projective(self, v, side="right"):
        """The indecomposable projective P_v: right ``e_v A`` (default) or, with
        ``side="left"``, the left projective ``A e_v`` (Plan 24)."""
        return self._sided_builder("projective", v, side)

    def injective(self, v, side="right"):
        """The indecomposable injective I_v = D(A e_v) (right, default) or the left
        injective with ``side="left"`` (Plan 24)."""
        return self._sided_builder("injective", v, side)

    def module(self, dimension_vector, arrow_action, side="right", name="M"):
        """Build a module from a dimension vector + one exact matrix per arrow
        (Plan 05 `Module.from_arrow_action`). ``side="right"`` (default) reads the
        matrices as a right A-module (arrow a: s->t acting M_s -> M_t); ``side="left"``
        builds a left A-module = right A^op-module, so the matrices are the opposite-
        quiver representation (Plan 24)."""
        from quiverlab.modules.module import Module
        if side == "right":
            return Module.from_arrow_action(self, dimension_vector, arrow_action, name=name)
        if side == "left":
            from quiverlab.modules.opposite import opposite_algebra
            m = Module.from_arrow_action(opposite_algebra(self), dimension_vector,
                                         arrow_action, name=name)
            return m.with_side("left")
        from quiverlab.errors import QuiverlabError
        raise QuiverlabError(f'side must be "right" or "left", got {side!r}')

    def opposite(self):
        """The opposite algebra A^op (reversed quiver, transposed structure
        constants), as a first-class Algebra. Involutive: A.opposite().opposite()
        is A (Plan 23)."""
        from quiverlab.modules.opposite import opposite_algebra
        return opposite_algebra(self)

    def hom(self, M, N):
        """dim Hom_A(M, N) for right A-modules M, N (spec §3.6)."""
        from quiverlab.modules.hom import hom_dim
        return hom_dim(M, N)

    def ext(self, M, N, n):
        """dim Ext^n_A(M, N) for right A-modules M, N (spec §3.6)."""
        from quiverlab.modules.ext import ext
        return ext(self, M, N, n)

    def ext_algebra(self, top=6):
        """The Yoneda / Ext-algebra E(A) = Ext^*_A(A/J, A/J) as a graded
        quiver-with-relations presentation over R = k^{Q_0}, through degree `top`
        (or complete through gl.dim when finite); a YonedaPresentation (Plan 27)."""
        from quiverlab.modules.ext_algebra import ext_algebra
        return ext_algebra(self, top)

    def crosscheck(self, what="hochschild", *args, **kwargs):
        """Independently recompute an invariant via the optional QPA backend and
        compare (spec §5 c.12). Requires `pip install quiverlab[qpa]`; raises
        QpaUnavailableError otherwise. Examples:
            A.crosscheck("hochschild", 3)          # HH^0..HH^3 vs QPA enveloping route
            A.crosscheck("module_ext", M, 4)       # Ext^0..Ext^4(M,M) vs QPA (self-Ext)
        Returns a CrosscheckReport; call .assert_agree() to fail loudly on mismatch."""
        from quiverlab.qpa.crosscheck import crosscheck as _cc
        return _cc(self, what, *args, **kwargs)

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

    def euler_form(self, d, e):
        """Euler bilinear form <d, e> = d C^{-1} e^T on integer dimension vectors
        (vertex order); for finite gl.dim, sum (-1)^i dim Ext^i (Plan 38 / C2)."""
        from quiverlab.invariants.forms import euler_form
        return euler_form(self, d, e)

    def tits_form(self, d):
        """Tits quadratic form q(d) = <d, d> of the Euler form (Plan 38 / C2)."""
        from quiverlab.invariants.forms import tits_form
        return tits_form(self, d)

    def form_type(self):
        """'finite' / 'tame' / 'wild' by exact definiteness of the Tits form; a
        representation-type theorem for hereditary algebras, the signature
        otherwise (Plan 38 / C2)."""
        from quiverlab.invariants.forms import form_type
        return form_type(self)

    def dynkin_type(self):
        """Orientation-blind Dynkin / Euclidean type of the underlying quiver:
        ("A"|"D"|"E"|"~A"|"~D"|"~E", n) or None (Plan 38 / C2). Loud if this
        algebra has no quiver presentation."""
        if self.quiver is None:
            from quiverlab.errors import QuiverlabError
            raise QuiverlabError(
                "dynkin_type needs the quiver presentation",
                hint="build the algebra via Quiver.algebra(...); structure-constant "
                     "algebras carry no quiver")
        from quiverlab.invariants.dynkin_type import dynkin_type
        return dynkin_type(self.quiver)

    def positive_roots(self):
        """Positive roots of the Tits form (= dimension vectors of the
        indecomposables, Gabriel) for a hereditary Dynkin algebra; loud on
        affine/wild/non-hereditary input (Plan 38 / C2)."""
        from quiverlab.invariants.roots import positive_roots
        return positive_roots(self)

    def loewy_length(self):
        """Loewy length = nilpotency index of rad A (exact, any field) (spec §3.5)."""
        from quiverlab.invariants.scalar import loewy_length
        return loewy_length(self)

    def center(self):
        """(dim, basis) of the center Z(A), exact over any field (spec §3.5)."""
        from quiverlab.invariants.scalar import center
        return center(self)

    def complexity(self, n):
        """Apparent complexity from the minimal A^e resolution's growth (GF(p) only).

        Exact for path-basis algebras of any vertex count (Plan 13: the engine builds
        the corner-typed minimal projective resolution). Remaining caveat: a
        memory-truncated build adds a silent prefix (read the number as a lower bound
        in that case). See ``invariants.scalar.complexity``.
        """
        from quiverlab.invariants.scalar import complexity
        return complexity(self, n)

    def draw(self, file=None):
        """Draw the quiver (matplotlib): vertices by depth, loops as self-arcs,
        parallel arrows fanned out, the relation list below (spec §3.7). Returns
        the Figure; pass file="A.png"/"A.svg" to also save it."""
        if self.quiver is None:
            from quiverlab.errors import QuiverlabError
            raise QuiverlabError(
                "this Algebra has no quiver to draw",
                hint="build it via Quiver.algebra(...) rather than from_structure_constants")
        from quiverlab.viz.draw import draw_quiver
        return draw_quiver(self.quiver, self.relations or [], file=file)

    def tikz(self):
        """TikZ source for the quiver, sharing draw()'s exact layout coordinates
        (spec §3.7). Paste into a LaTeX document."""
        if self.quiver is None:
            from quiverlab.errors import QuiverlabError
            raise QuiverlabError(
                "this Algebra has no quiver to render as TikZ",
                hint="build it via Quiver.algebra(...) rather than from_structure_constants")
        from quiverlab.viz.tikz import tikz_quiver
        return tikz_quiver(self.quiver, self.relations or [])

    def cyclic_homology(self, top, max_cells=4_000_000, with_reps=False):
        """Dimensions of HC_0..HC_top (Connes (b, B) mixed complex).

        GF(p): the fast engine (int64 rank). Any other exact Domain: the
        generic mixed complex on the normalized bar basis. Both are exponential
        in ``top`` (dim C_n = m*(m-1)^n); ``max_cells`` guards the blow-up on
        BOTH paths, refusing loudly (DepthLimitError) before any matrix is
        allocated — raise it to compute a bigger case. Works for any unital
        algebra.

        Plan 35 wave 3b — ``with_reps=True`` returns ``(table, payload)`` where
        ``payload`` carries the explicit HC representatives (``basis_classes`` /
        ``chain_basis`` / ``differentials`` / ``column_structure`` keyed by
        ``str(degree)``) captured from the SAME total complex; the default path is
        byte-unchanged."""
        from quiverlab.fields.primefield import PrimeField
        if isinstance(self.domain, PrimeField):
            from quiverlab.engine.adapter import to_engine
            from quiverlab.engine.cyclic import cyclic_homology_dims
            from quiverlab.hochschild.table import HHTable
            p = self.domain.p
            AU = self.unit_adapted()
            E = to_engine(AU)
            res = cyclic_homology_dims(E, top, primes=(p,), with_reps=with_reps,
                                       max_cells=max_cells)
            out, raw = res if with_reps else (res, None)
            dims = [int(d) for d in out[p]]
            table = HHTable(dims, "HC_", repr(self).splitlines()[0],
                            engine="hanlab engine (F_p fast rank)")
            if not with_reps:
                return table
            from quiverlab.hochschild.cyclic_reps import gfp_payload
            return table, gfp_payload(E, AU, top, raw, dims)
        from quiverlab.hochschild.cyclic import cyclic_homology_dims
        if not with_reps:
            return cyclic_homology_dims(self, top, max_cells=max_cells)
        table, raw = cyclic_homology_dims(self, top, max_cells=max_cells, with_reps=True)
        from quiverlab.hochschild.cyclic_reps import generic_payload
        return table, generic_payload(self, top, raw, list(table.dims))

    def _product_dispatch(self, kind, top, engine, max_cells):
        """Shared Plan-35 routing: GF(p) -> bar/tt tables; quiver-presented ->
        CS native (any Domain; also the DepthLimitError fallback target);
        presentation-less off GF(p) -> loud refusal. One engine end-to-end."""
        from quiverlab.errors import DepthLimitError
        from quiverlab.fields.primefield import PrimeField
        from quiverlab.hochschild.products import gfp_product_tables
        if engine not in ("auto", "bar", "cs"):
            raise QuiverlabError(
                f"unknown engine {engine!r} for {kind} tables",
                hint="choose 'auto', 'bar', or 'cs'")
        is_gfp = isinstance(self.domain, PrimeField)
        presented = self.quiver is not None and self.relations is not None
        if engine == "cs" or (engine == "auto" and not is_gfp):
            if kind == "bracket":
                raise QuiverlabError(
                    "the Gerstenhaber bracket is served over GF(p) only "
                    "(bar window; no CS-native brace machinery in v1)",
                    hint="construct the algebra over GF(p)")
            if not presented:
                raise QuiverlabError(
                    f"{kind} tables off GF(p) need a quiver presentation "
                    "(the CS route); this algebra has structure constants only",
                    hint="build the algebra via Quiver.algebra, or use GF(p)")
            from quiverlab.resolutions_cs.products import cs_product_tables
            return cs_product_tables(self, kind, top, max_cells)
        if not is_gfp:            # engine == "bar" explicitly, off GF(p)
            raise QuiverlabError(
                f"engine='bar' {kind} tables need GF(p) (the tt facade)",
                hint="use engine='auto' (routes CS for presented algebras)")
        try:
            return gfp_product_tables(self, kind, top, max_cells)
        except DepthLimitError:
            if engine != "auto" or not presented or kind == "bracket":
                raise
            from quiverlab.resolutions_cs.products import cs_product_tables
            return cs_product_tables(self, kind, top, max_cells)

    def cup_products(self, top, engine="auto", max_cells=4_000_000):
        """Structure-constant tables of the cup product HH^p (x) HH^q ->
        HH^{p+q} for every p+q <= top, on the recorded basis. Exact. engine:
        'auto' (GF(p) -> bar/tt, else CS for presented algebras, with the CS
        depth fallback), 'bar' (GF(p) tt facade, loud otherwise), 'cs'
        (Chouhy-Solotar native diagonal, presented algebras, any Domain)."""
        return self._product_dispatch("cup", top, engine, max_cells)

    def cap_products(self, top, engine="auto", max_cells=4_000_000):
        """Structure-constant tables of the cap action HH^p (x) HH_n ->
        HH_{n-p} for p <= n <= top. Same engine semantics as cup_products."""
        return self._product_dispatch("cap", top, engine, max_cells)

    def gerstenhaber_brackets(self, top, engine="auto", max_cells=4_000_000):
        """Structure-constant tables of the Gerstenhaber bracket HH^p (x)
        HH^q -> HH^{p+q-1} for pairs p, q >= 1 with p+q-1 <= top. GF(p) only
        and window-bounded (the result records the served window); the
        degree-0 insertion action is out of scope."""
        return self._product_dispatch("bracket", top, engine, max_cells)

    def connes_differentials(self, top, max_cells=4_000_000):
        """Induced Connes differentials B : HH_n -> HH_{n+1} (matrices +
        ranks) for 0 <= n < top. GF(p) via the engine (b,B); any other exact
        Domain via the generic mixed complex — no engine choice to make."""
        from quiverlab.hochschild.products import connes_b_tables
        return connes_b_tables(self, top, max_cells=max_cells)

    def nakayama_automorphism(self):
        """Nakayama automorphism nu as a matrix (columns = images) in the
        algebra's basis. GF(p): integer matrix via the engine (unit-adapted
        basis). Other exact Domains: Domain-element matrix on the path-type
        basis (Plan 19). Loud if not Frobenius."""
        from quiverlab.fields.primefield import PrimeField
        if isinstance(self.domain, PrimeField):
            from quiverlab.engine.adapter import to_engine
            from quiverlab.engine.coxeter import nakayama_automorphism
            S, _ = nakayama_automorphism(to_engine(self.unit_adapted()), self.domain.p)
            return [[int(S[i, j]) for j in range(S.shape[1])] for i in range(S.shape[0])]
        from quiverlab.invariants.frobenius import nakayama_automorphism_generic
        return nakayama_automorphism_generic(self)

    def is_frobenius(self):
        """Is the algebra Frobenius? GF(p): engine form search. Other exact
        Domains: the exact socle criterion on a path-type basis (Plan 19)."""
        from quiverlab.fields.primefield import PrimeField
        if isinstance(self.domain, PrimeField):
            from quiverlab.engine.adapter import to_engine
            from quiverlab.engine.coxeter import is_frobenius
            return bool(is_frobenius(to_engine(self.unit_adapted()), self.domain.p))
        from quiverlab.invariants.frobenius import is_frobenius_generic
        return is_frobenius_generic(self)

    def is_symmetric(self):
        """Is the algebra symmetric — does it carry a nondegenerate trace form
        lambda(ab) = lambda(ba) (equivalently A ~= DA as bimodules)? Exact over
        every exact Domain via the trace-form certificate on the path-type basis
        (Plan 29; Skowronski–Yamagata). Loud refusal on a presentation-less
        algebra, like the other path-basis invariants.

        Plan 29 replaced the former GF(p) shortcut ``is_frobenius and the
        engine's Nakayama automorphism == identity matrix``, which was
        sufficient-not-necessary and returned a silent wrong False on
        multi-vertex symmetric Nakayama (Brauer star) algebras."""
        from quiverlab.invariants.frobenius import is_symmetric_generic
        return is_symmetric_generic(self)

    def is_weakly_symmetric(self):
        """Is the algebra weakly symmetric — Frobenius with the identity Nakayama
        permutation (soc P_v = top P_v for every indecomposable projective)?
        Exact over every exact Domain (Plan 29). For self-injective Nakayama
        kZ_n/J^L this is n | (L - 1) (Skowronski–Yamagata, Frobenius Algebras I).
        Every symmetric algebra is weakly symmetric; the converse can fail."""
        from quiverlab.invariants.frobenius import is_weakly_symmetric_generic
        return is_weakly_symmetric_generic(self)

    def tor(self, M, N, n):
        """dim Tor_n^A(M, N) for a RIGHT A-module M and a LEFT A-module N (Plan 29).

        The homological sibling of ``ext``: H_n(P_* (x)_A N) with P_* the minimal
        projective resolution of the right module M. Certified against ``ext`` by the
        duality dim Tor_n(M, N) = dim Ext^n(M, DN) (D side-aware, Plan 24)."""
        from quiverlab.modules.tor import tor
        return tor(self, M, N, n)

    def __repr__(self):
        base = f"Algebra of dimension {self.dim} over {self.domain.name}"
        lines = [base]
        if self.basis_labels:
            lines.append("basis: " + ", ".join(self.basis_labels))
        q = self.quiver
        if q is not None:  # spec 3.7: plain-text shows vertices, arrows, relations
            lines.append("vertices: " + ", ".join(str(v) for v in q.vertices))
            lines.append("arrows: " + "; ".join(
                f"{n}: {s} -> {t}" for n, (s, t) in q.arrows.items()))
            rels = self.relations or []
            if rels:
                lines.append("relations: " + "; ".join(repr(r) for r in rels))
        return "\n".join(lines)
