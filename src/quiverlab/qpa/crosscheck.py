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


def crosscheck(algebra, what: str, *args, **kwargs) -> CrosscheckReport:
    """Dispatch. what="hochschild"|"module_ext" (Plan 08); "tau"|"tau_minus"|
    "proj_resolution"|"inj_resolution"|"inj_dimension" (Plan 23 module surface)."""
    if what == "hochschild":
        return crosscheck_hochschild(algebra, *args, **kwargs)
    if what == "module_ext":
        return crosscheck_module_ext(algebra, *args, **kwargs)
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
    # An unrecognized `what` is a usage error, NOT "QPA unavailable".
    raise QuiverlabError(f"unknown cross-check {what!r}",
                         hint='use "hochschild", "module_ext", "tau", "tau_minus", '
                              '"proj_resolution", "inj_resolution", or "inj_dimension"')
