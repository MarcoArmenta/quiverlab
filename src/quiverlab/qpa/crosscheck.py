"""A.crosscheck(...): independent QPA recomputation of Hochschild dims / module Ext
for validation workflows (spec §5 c.12, §8 ring 3). Returns a CrosscheckReport;
never silently disagrees -- .assert_agree() raises on mismatch."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QpaUnavailableError, QuiverlabError
from quiverlab.qpa import scripts, session


@dataclass
class CrosscheckReport:
    what: str                 # "hochschild" | "module_ext"
    ours: list                # quiverlab dims
    qpa: list                 # QPA dims
    agree: bool

    def assert_agree(self):
        if not self.agree:
            raise AssertionError(
                f"QPA cross-check DISAGREES on {self.what}: quiverlab {self.ours} "
                f"vs QPA {self.qpa}")
        return self


def _read_int_list(gap_value) -> list:
    """Convert a GAP list of integers into a Python list[int] (exact; no floats)."""
    return [int(x) for x in gap_value]


def crosscheck_hochschild(algebra, top: int) -> CrosscheckReport:
    session.require_gap()
    ours = algebra.hochschild_cohomology(top).dims
    # the trailing read-back must be its OWN line: session.run evals per line
    gap = session.run(scripts.hochschild_dims_script(algebra, top) + "\nhh;")
    qpa = _read_int_list(gap)
    return CrosscheckReport("hochschild", list(ours), qpa, list(ours) == qpa)


def crosscheck_module_ext(algebra, M, top: int) -> CrosscheckReport:
    """Self-Ext Ext^*(M, M) vs QPA (via ExtAlgebraGenerators). Distinct-module
    Ext(M, N) is a flagged post-v1 extension (needs ExtOverAlgebra + syzygies)."""
    session.require_gap()
    ours = [algebra.ext(M, M, n) for n in range(top + 1)]   # ext() returns dim (int)
    dimvec = M.dimension_vector()                           # dict {vertex: dim}
    dims = [dimvec[v] for v in algebra.quiver.vertices]     # QPA order = quiver order
    gap = session.run(
        scripts.module_self_ext_dims_script(algebra, dims, top) + "\next;")
    qpa = _read_int_list(gap)
    return CrosscheckReport("module_ext", list(ours), qpa, list(ours) == qpa)


def crosscheck_symmetric(algebra) -> CrosscheckReport:
    """``is_symmetric`` / ``is_weakly_symmetric`` vs QPA's ``IsSymmetricAlgebra`` /
    ``IsWeaklySymmetricAlgebra`` (Plan 29). Both verdicts are compared at once; ``ours``
    and ``qpa`` are ``[symmetric, weakly_symmetric]`` boolean pairs. Reproduces the fixed
    live bug: multi-vertex symmetric Nakayama (Brauer star) algebras kZ_n/J^L with
    n | (L-1) that the former engine shortcut wrongly reported non-symmetric.

    Scope: the algebra must carry a quiver presentation over QQ or a prime GF(p) (QPA's
    ``kQ/rels`` route). A presented ``TrivialExtension(A)`` (Plan 31 double-quiver build)
    qualifies and can be fed here directly; the ``trivial_extension`` crosscheck below
    additionally compares it against QPA's native ``TrivialExtensionOfQuiverAlgebra``."""
    session.require_gap()
    ours = [bool(algebra.is_symmetric()), bool(algebra.is_weakly_symmetric())]
    base = scripts.symmetric_predicates_script(algebra)
    qpa = [bool(session.run(base + "\nIsSymmetricAlgebra(A);")),
           bool(session.run(base + "\nIsWeaklySymmetricAlgebra(A);"))]
    return CrosscheckReport("symmetric", ours, qpa, ours == qpa)


def crosscheck_trivial_extension(algebra) -> CrosscheckReport:
    """Our certified double-quiver ``TrivialExtension(algebra)`` vs QPA's native
    ``TrivialExtensionOfQuiverAlgebra`` (Plan 31). Compares the 5-tuple
    ``[dim, #arrows, IsSymmetric, IsWeaklySymmetric, IsSelfinjective]``.

    ``algebra`` (the BASE A) must carry a quiver presentation over QQ or a prime
    GF(p); our TrivialExtension then takes the presented route, so ``B.quiver`` is
    non-None and its arrows are countable. QPA labels its dual arrows differently
    (``te_a1_i_j``), so the arrow COUNT is compared, not names."""
    session.require_gap()
    from quiverlab.families.trivial_extension import TrivialExtension
    B = TrivialExtension(algebra)
    if B.quiver is None:
        raise QuiverlabError(
            "crosscheck_trivial_extension needs a presented base algebra",
            hint="give a Quiver(...).algebra(...) over QQ or a prime GF(p) so "
                 "TrivialExtension takes its double-quiver route")
    ours = [B.dim, len(B.quiver.arrows), bool(B.is_symmetric()),
            bool(B.is_weakly_symmetric()), bool(B.is_selfinjective())]
    base = scripts.trivial_extension_script(algebra)
    qpa = [int(session.run(base + "\nDimension(TE);")),
           int(session.run(base + "\nNumberOfArrows(QuiverOfPathAlgebra(TE));")),
           bool(session.run(base + "\nIsSymmetricAlgebra(TE);")),
           bool(session.run(base + "\nIsWeaklySymmetricAlgebra(TE);")),
           bool(session.run(base + "\nIsSelfinjectiveAlgebra(TE);"))]
    return CrosscheckReport("trivial_extension", ours, qpa, ours == qpa)


@dataclass
class ModuleCrosscheckReport:
    """A module-level QPA crosscheck (tau/tau^-, resolutions, injective dimension).
    `iso` is None when no isomorphism-class check applies."""
    what: str
    ours: object
    qpa: object
    agree: bool
    iso: object = None

    def assert_agree(self):
        if not self.agree:
            raise AssertionError(
                f"QPA cross-check DISAGREES on {self.what}: quiverlab {self.ours} "
                f"vs QPA {self.qpa}" + ("" if self.iso is None else f" (iso={self.iso})"))
        return self


def _dv_list(algebra, M):
    dv = M.dimension_vector()
    return [dv[v] for v in algebra.quiver.vertices]


def _dv_list_of_dict(algebra, d):
    return [d.get(v, 0) for v in algebra.quiver.vertices]


def _graded(algebra, M):
    from quiverlab.modules.qpa_module import graded_form
    return graded_form(M)


def crosscheck_tau(algebra, M, minus: bool = False) -> ModuleCrosscheckReport:
    """tau M (or tau^- M) vs QPA DTr/TrD: dimension vectors AND isomorphism class
    (our translate is emitted as a QPA module and compared via IsomorphicModules)."""
    session.require_gap()
    what = "tau_minus" if minus else "tau"
    op = "TrD" if minus else "DTr"
    tr = M.tau_minus() if minus else M.tau()
    ours_dv = _dv_list(algebra, tr)
    dvM, arrM = _graded(algebra, M)
    base = scripts.quiver_and_algebra_script(algebra)
    base += "\n" + scripts.module_decl(algebra, dvM, arrM, "M")
    qpa_dv = _read_int_list(session.run(base + f"\nt := {op}(M);;\nDimensionVector(t);"))
    # isomorphism class: build our translate as a QPA module, compare to op(M)
    dvT, arrT = _graded(algebra, tr)
    iso_script = base + "\n" + scripts.module_decl(algebra, dvT, arrT, "TR")
    iso_val = session.run(iso_script + f"\nt := {op}(M);;\nIsomorphicModules(TR, t);")
    iso = bool(iso_val)
    return ModuleCrosscheckReport(what, ours_dv, qpa_dv, ours_dv == qpa_dv and iso, iso)


def crosscheck_proj_resolution(algebra, M, top: int) -> ModuleCrosscheckReport:
    """Projective resolution term dimension vectors vs QPA ProjectiveResolution."""
    session.require_gap()
    res = M.projective_resolution(top)
    dvs = res.dimension_vectors()
    ours = [_dv_list_of_dict(algebra, dvs[n]) if n < len(dvs) else [0] * len(_dv_list(algebra, M))
            for n in range(top + 1)]
    dvM, arrM = _graded(algebra, M)
    base = scripts.quiver_and_algebra_script(algebra)
    base += "\n" + scripts.module_decl(algebra, dvM, arrM, "M") + "\npr := ProjectiveResolution(M);;"
    qpa = [_read_int_list(session.run(base + f"\nDimensionVector(ObjectOfComplex(pr,{n}));"))
           for n in range(top + 1)]
    return ModuleCrosscheckReport("proj_resolution", ours, qpa, ours == qpa)


def crosscheck_inj_resolution(algebra, M, top: int) -> ModuleCrosscheckReport:
    """Injective resolution term dimension vectors vs QPA: the defining identity
    inj.res_A(M) <-> proj.res_{A^op}(DM) (term dim-vectors agree since D preserves
    dimension vectors)."""
    session.require_gap()
    res = M.injective_resolution(top)
    dvs = res.dimension_vectors()
    ours = [_dv_list_of_dict(algebra, dvs[n]) if n < len(dvs) else [0] * len(_dv_list(algebra, M))
            for n in range(top + 1)]
    dvM, arrM = _graded(algebra, M)
    base = scripts.quiver_and_algebra_script(algebra)
    base += "\n" + scripts.module_decl(algebra, dvM, arrM, "M")
    base += "\nDM := DualOfModule(M);;\npr := ProjectiveResolution(DM);;"
    qpa = [_read_int_list(session.run(base + f"\nDimensionVector(ObjectOfComplex(pr,{n}));"))
           for n in range(top + 1)]
    return ModuleCrosscheckReport("inj_resolution", ours, qpa, ours == qpa)


def crosscheck_inj_dimension(algebra, M, bound: int) -> ModuleCrosscheckReport:
    """Injective dimension vs QPA InjDimensionOfModule (false <-> infinite <-> None)."""
    session.require_gap()
    ours = M.injective_dimension(bound=bound)
    dvM, arrM = _graded(algebra, M)
    base = scripts.quiver_and_algebra_script(algebra)
    base += "\n" + scripts.module_decl(algebra, dvM, arrM, "M")
    val = session.run(base + f"\nInjDimensionOfModule(M, {bound});")
    try:
        qpa = int(val)
    except (TypeError, ValueError):
        qpa = None                           # GAP `false` => injective dimension > bound
    return ModuleCrosscheckReport("inj_dimension", ours, qpa, ours == qpa)


# ---------------------------------------------------------------------------
# Plan 40: homological-dimension family crosschecks. QPA exposes
# GlobalDimensionOfAlgebra / DominantDimensionOfAlgebra / GorensteinDimensionOfAlgebra
# (all bound-parametrised, returning an int or GAP `infinity`); it exposes NO
# Igusa-Todorov surface (a live NamesGVars() probe finds none -- test_homdims_qpa),
# so phi/psi have no QPA oracle (the Task-B literature battery covers them). The
# verbs are inlined here, like crosscheck_inj_dimension, over quiver_and_algebra_script.
# ---------------------------------------------------------------------------
def _qpa_dim_or_none(val):
    """A QPA homological dimension as ``int`` or ``None`` -- GAP ``infinity`` / ``false``
    (dimension beyond the bound, i.e. our infinite / unresolved marker) both map to
    ``None``, mirroring crosscheck_inj_dimension."""
    try:
        return int(val)
    except (TypeError, ValueError):
        return None


def crosscheck_global_dimension(algebra, bound: int = 20) -> CrosscheckReport:
    """gl.dim vs QPA ``GlobalDimensionOfAlgebra(A, n)`` (int, or ``infinity`` -> None).
    Our unresolved/infinite verdict (not exact) maps to ``None`` on both sides."""
    session.require_gap()
    g = algebra.global_dimension()
    ours = g.value if g.exact else None
    base = scripts.quiver_and_algebra_script(algebra)
    qpa = _qpa_dim_or_none(session.run(base + f"\nGlobalDimensionOfAlgebra(A, {bound});"))
    return CrosscheckReport("global_dimension", ours, qpa, ours == qpa)


def crosscheck_dominant_dimension(algebra, bound: int = 20) -> CrosscheckReport:
    """Dominant dimension vs QPA ``DominantDimensionOfAlgebra(A, n)`` (int, or
    ``infinity`` -> None). Our ``infinite`` verdict (self-injective) maps to ``None``."""
    session.require_gap()
    dd = algebra.dominant_dimension()
    ours = None if dd.infinite else dd.value
    base = scripts.quiver_and_algebra_script(algebra)
    qpa = _qpa_dim_or_none(session.run(base + f"\nDominantDimensionOfAlgebra(A, {bound});"))
    return CrosscheckReport("dominant_dimension", ours, qpa, ours == qpa)


def crosscheck_gorenstein(algebra, bound: int = 20) -> CrosscheckReport:
    """Gorenstein dimension vs QPA ``GorensteinDimensionOfAlgebra(A, n)`` (int, or
    ``infinity`` -> None). Ours = ``max(right inj.dim, left inj.dim)`` when Gorenstein
    (``is_gorenstein`` True), else ``None`` (the bounded engine did not prove finiteness
    -- QPA likewise returns ``infinity``)."""
    session.require_gap()
    gd = algebra.gorenstein_dimension()
    ours = (max(gd.right_id, gd.left_id) if gd.is_gorenstein else None)
    base = scripts.quiver_and_algebra_script(algebra)
    qpa = _qpa_dim_or_none(session.run(base + f"\nGorensteinDimensionOfAlgebra(A, {bound});"))
    return CrosscheckReport("gorenstein", ours, qpa, ours == qpa)


# ---------------------------------------------------------------------------
# Plan 30: Krull-Schmidt decomposition crosschecks (finite fields only -- QPA's
# DecomposeModule requires GF(p))
# ---------------------------------------------------------------------------
def _flat_dimvec_multiset(algebra, pairs):
    """A ``(summand_dimvec, multiplicity)`` list expanded into the SORTED multiset of
    per-summand dimension-vector tuples (each summand repeated by its multiplicity).
    This is the order- and grouping-independent invariant compared against QPA -- robust
    even if two non-isomorphic summands happen to share a dimension vector."""
    flat = []
    for dv, mult in pairs:
        flat.extend([tuple(dv)] * mult)
    return sorted(flat)


def crosscheck_decompose(algebra, M) -> ModuleCrosscheckReport:
    """Krull-Schmidt summand dimension-vectors + multiplicities vs QPA's
    ``DecomposeModuleWithMultiplicities`` (Plan 30). Both sides are reduced to the sorted
    multiset of per-indecomposable-summand dimension vectors (expanded by multiplicity),
    so ordering and grouping conventions do not matter. GF(p) only (QPA requirement)."""
    session.require_gap()
    ours_pairs = [(_dv_list(algebra, s), mult) for s, mult in M.decompose()]
    ours = _flat_dimvec_multiset(algebra, ours_pairs)
    dvM, arrM = _graded(algebra, M)
    base = scripts.decompose_multiplicities_script(algebra, dvM, arrM)
    qpa_dimvecs = [_read_int_list(row)
                   for row in session.run(base + "\nList(d[1], DimensionVector);")]
    qpa_mults = _read_int_list(session.run(base + "\nd[2];"))
    qpa = _flat_dimvec_multiset(algebra, list(zip(qpa_dimvecs, qpa_mults)))
    return ModuleCrosscheckReport("decompose", ours, qpa, ours == qpa)


def crosscheck_indecomposable(algebra, M) -> ModuleCrosscheckReport:
    """``M.is_indecomposable()`` vs QPA's ``IsIndecomposableModule`` (Plan 30). GF(p)
    only. (Over ``char <= dim M`` our engine may refuse loudly rather than answer -- run
    the crosscheck over a characteristic ``> dim M`` so both sides decide.)"""
    session.require_gap()
    ours = M.is_indecomposable()
    dvM, arrM = _graded(algebra, M)
    base = scripts.is_indecomposable_script(algebra, dvM, arrM)
    qpa = bool(session.run(base + "\nIsIndecomposableModule(M);"))
    return ModuleCrosscheckReport("indecomposable", ours, qpa, ours == qpa)


# ---------------------------------------------------------------------------
# Plan 37: C1 categorical glue -- Hom dims (+ kernel/image/cokernel dim-vectors
# where QPA exposes them: KernelInclusion / ImageInclusion / CoKernelProjection)
# ---------------------------------------------------------------------------
def crosscheck_hom_glue(algebra, M, N) -> ModuleCrosscheckReport:
    """``dim Hom_A(M, N)`` (our ``hom_basis``) vs QPA ``Length(HomOverAlgebra(M, N))``,
    and -- whenever ``dim Hom = 1`` (the nonzero hom is then unique up to a scalar, so
    its kernel/image/cokernel dimension vectors are well-defined invariants both bases
    agree on) -- the kernel/image/cokernel dimension vectors of that canonical hom vs
    QPA's ``KernelInclusion`` / ``ImageInclusion`` / ``CoKernelProjection`` (QPA Ch. 7).
    ``ours`` / ``qpa`` are dicts keyed ``hom_dim`` (always) and, when comparable,
    ``kernel`` / ``image`` / ``cokernel`` (dimension vectors in quiver-vertex order)."""
    session.require_gap()
    ours_hom = algebra.hom(M, N)
    dvM, arrM = _graded(algebra, M)
    dvN, arrN = _graded(algebra, N)
    base = scripts.hom_glue_script(algebra, dvM, arrM, dvN, arrN)
    qpa_hom = int(session.run(base + "\nLength(homs);"))
    ours = {"hom_dim": ours_hom}
    qpa = {"hom_dim": qpa_hom}
    if ours_hom == 1 and qpa_hom == 1:
        f = algebra.hom_basis(M, N)[0]
        K, _ = f.kernel()
        I, _, _ = f.image()
        C, _ = f.cokernel()
        ours["kernel"] = _dv_list(algebra, K)
        ours["image"] = _dv_list(algebra, I)
        ours["cokernel"] = _dv_list(algebra, C)
        qpa["kernel"] = _read_int_list(session.run(
            base + "\nf := homs[1];;\nDimensionVector(Source(KernelInclusion(f)));"))
        qpa["image"] = _read_int_list(session.run(
            base + "\nf := homs[1];;\nDimensionVector(Source(ImageInclusion(f)));"))
        qpa["cokernel"] = _read_int_list(session.run(
            base + "\nf := homs[1];;\nDimensionVector(Range(CoKernelProjection(f)));"))
    return ModuleCrosscheckReport("hom_glue", ours, qpa, ours == qpa)


# ---------------------------------------------------------------------------
# Plan 27: the Yoneda / Ext-algebra E(A) = Ext^*(A/J, A/J) crosschecks
# ---------------------------------------------------------------------------
def crosscheck_ext_algebra_dims(algebra, top: int) -> CrosscheckReport:
    """Total graded Ext-algebra dimensions ``dim E^n = sum_{i,j} dim Ext^n(S_i, S_j)``
    (our engine, summed over corners) vs QPA's ``ExtAlgebraGenerators(M, top)[1]`` for
    ``M = (+) SimpleModules(A)`` (Plan 27). Both indexed n = 0..top; degree 0 is
    ``|Q_0|`` on both sides."""
    session.require_gap()
    ours = algebra.ext_algebra(top=top).graded_dims_through(top)
    script = scripts.ext_algebra_generators_script(algebra, top)
    qpa = _read_int_list(session.run(script + "\ninfo[1];"))
    return CrosscheckReport("ext_algebra_dims", list(ours), qpa, list(ours) == qpa)


def crosscheck_ext_generator_degrees(algebra, top: int) -> CrosscheckReport:
    """Per-degree counts of NEW minimal E-algebra generators vs QPA's
    ``ExtAlgebraGenerators(M, top)[2]`` (Plan 27).

    DEGREE-0 MAPPING (documented, honest): our engine treats degree 0 as the semisimple
    base R = k^{Q_0} (the augmentation ideal starts in degree 1) and so does NOT list it
    in ``generators_by_degree``; QPA instead counts the ``|Q_0|`` vertex idempotents as
    the degree-0 generators. We align by reporting our degree-0 count as
    ``len(A.quiver.vertices)``; for n >= 1 the count is
    ``len(generators_by_degree[n])``. This makes the Koszulity discriminator explicit:
    rad^2 = 0 A_3 and the cubic A_4 BOTH have ``dim Ext^2 = 1``, but A_3's degree-2
    class is decomposable (generator count 0) while A_4's is a genuine new generator
    (count 1)."""
    session.require_gap()
    P = algebra.ext_algebra(top=top)
    nv = len(algebra.quiver.vertices)
    ours = [nv] + [len(P.generators_by_degree.get(n, [])) for n in range(1, top + 1)]
    script = scripts.ext_algebra_generators_script(algebra, top)
    qpa = _read_int_list(session.run(script + "\ninfo[2];"))
    return CrosscheckReport("ext_generator_degrees", ours, qpa, ours == qpa)


def crosscheck_ext_quiver(algebra) -> CrosscheckReport:
    """The Ext-quiver / Ext^1 corner matrix ``dim Ext^1(S_i, S_j)`` (our corner (i, j))
    vs QPA's ``Length(ExtOverAlgebra(S_i, S_j)[2])`` (Plan 27). Pins the corner /
    direction convention: our corner (i, j) must equal QPA's ``ExtOverAlgebra(S_i, S_j)``.
    Both matrices are indexed in quiver-vertex order."""
    session.require_gap()
    eng = algebra.ext_algebra(top=2)._eng
    verts = list(algebra.quiver.vertices)
    ours = [[eng.ext_dim(vi, 1, vj) for vj in verts] for vi in verts]
    base = scripts.ext_quiver_script(algebra)
    qpa = [[int(session.run(base + "\n" + scripts.ext_quiver_entry(i + 1, j + 1)))
            for j in range(len(verts))] for i in range(len(verts))]
    return CrosscheckReport("ext_quiver", ours, qpa, ours == qpa)


def crosscheck_quadratic(algebra) -> CrosscheckReport:
    """``modules.koszul.is_quadratic(A)`` (defining ideal generated in degree 2) vs
    QPA's ``IsQuadraticIdeal(rels)`` (Plan 27)."""
    session.require_gap()
    from quiverlab.modules.koszul import is_quadratic
    ours = is_quadratic(algebra)
    qpa = bool(session.run(scripts.quadratic_ideal_script(algebra) + "\nisquad;"))
    return CrosscheckReport("quadratic", ours, qpa, ours == qpa)


def crosscheck_koszul_derived(algebra, top: int) -> CrosscheckReport:
    """The DERIVED Koszul verdict (Plan 27). QPA ships no ``IsKoszul``/``KoszulDual``,
    so the QPA-side verdict is derived from primitives it DOES compute: A is
    (QPA-)derived-Koszul iff ``IsQuadraticIdeal(rels)`` AND no new E-algebra generator
    appears in any degree 2..top (``ExtAlgebraGenerators(M, top)[2][n] = 0`` for
    n >= 2). This is compared against our ``YonedaPresentation.koszul``:

    * ours is not None (a certified True/False): agree iff ``ours == qpa_derived``.
    * ours is None ("no obstruction through degree top", not certified): we only
      require QPA NOT CONTRADICT us -- i.e. QPA must not exhibit a degree->=2 generator
      (which would be a certain not-Koszul that our engine, computing the same
      generators, would also have caught as False). Note the derived verdict is the
      honest one QPA can back; the certifier itself (G-quadratic / Priddy PBW) is our
      theory oracle, which QPA cannot compare."""
    session.require_gap()
    ours = algebra.ext_algebra(top=top).koszul
    script = scripts.ext_algebra_generators_script(algebra, top)
    gens = _read_int_list(session.run(script + "\ninfo[2];"))
    qpa_high_gen = any(gens[n] for n in range(2, len(gens)))
    qpa_quadratic = bool(session.run(scripts.quadratic_ideal_script(algebra) + "\nisquad;"))
    qpa_derived = qpa_quadratic and not qpa_high_gen
    agree = (not qpa_high_gen) if ours is None else (ours == qpa_derived)
    return CrosscheckReport("koszul_derived", ours, qpa_derived, agree)


def crosscheck(algebra, what: str, *args, **kwargs) -> CrosscheckReport:
    """Dispatch. what="hochschild"|"module_ext" (Plan 08); "symmetric" (Plan 29);
    "trivial_extension" (Plan 31); "tau"|"tau_minus"|"proj_resolution"|
    "inj_resolution"|"inj_dimension" (Plan 23 module surface);
    "ext_algebra_dims"|"ext_generator_degrees"|"ext_quiver"|"quadratic"|
    "koszul_derived" (Plan 27 Yoneda/Ext-algebra)."""
    if what == "hochschild":
        return crosscheck_hochschild(algebra, *args, **kwargs)
    if what == "module_ext":
        return crosscheck_module_ext(algebra, *args, **kwargs)
    if what == "symmetric":
        return crosscheck_symmetric(algebra, *args, **kwargs)
    if what == "trivial_extension":
        return crosscheck_trivial_extension(algebra, *args, **kwargs)
    if what == "tau":
        return crosscheck_tau(algebra, *args, minus=False, **kwargs)
    if what == "tau_minus":
        return crosscheck_tau(algebra, *args, minus=True, **kwargs)
    if what == "proj_resolution":
        return crosscheck_proj_resolution(algebra, *args, **kwargs)
    if what == "inj_resolution":
        return crosscheck_inj_resolution(algebra, *args, **kwargs)
    if what == "inj_dimension":
        return crosscheck_inj_dimension(algebra, *args, **kwargs)
    if what == "decompose":
        return crosscheck_decompose(algebra, *args, **kwargs)
    if what == "indecomposable":
        return crosscheck_indecomposable(algebra, *args, **kwargs)
    if what == "hom_glue":
        return crosscheck_hom_glue(algebra, *args, **kwargs)
    if what == "ext_algebra_dims":
        return crosscheck_ext_algebra_dims(algebra, *args, **kwargs)
    if what == "ext_generator_degrees":
        return crosscheck_ext_generator_degrees(algebra, *args, **kwargs)
    if what == "ext_quiver":
        return crosscheck_ext_quiver(algebra, *args, **kwargs)
    if what == "quadratic":
        return crosscheck_quadratic(algebra, *args, **kwargs)
    if what == "koszul_derived":
        return crosscheck_koszul_derived(algebra, *args, **kwargs)
    # An unrecognized `what` is a usage error, NOT "QPA unavailable".
    raise QuiverlabError(f"unknown cross-check {what!r}",
                         hint='use "hochschild", "module_ext", "symmetric", '
                              '"trivial_extension", "tau", "tau_minus", '
                              '"proj_resolution", "inj_resolution", '
                              '"inj_dimension", "decompose", "indecomposable", '
                              '"ext_algebra_dims", "ext_generator_degrees", '
                              '"ext_quiver", "quadratic", or "koszul_derived"')
