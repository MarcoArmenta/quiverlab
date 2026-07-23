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


def crosscheck(algebra, what: str, *args, **kwargs) -> CrosscheckReport:
    """Dispatch. what="hochschild" -> crosscheck_hochschild(algebra, top);
    what="module_ext" -> crosscheck_module_ext(algebra, M, top) (self-Ext)."""
    if what == "hochschild":
        return crosscheck_hochschild(algebra, *args, **kwargs)
    if what == "module_ext":
        return crosscheck_module_ext(algebra, *args, **kwargs)
    # An unrecognized `what` is a usage error, NOT "QPA unavailable".
    raise QuiverlabError(f"unknown cross-check {what!r}",
                         hint='use "hochschild" or "module_ext"')
