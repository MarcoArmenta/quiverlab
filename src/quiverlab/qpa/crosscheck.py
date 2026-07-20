"""A.crosscheck(...): independent QPA recomputation of Hochschild dims / module Ext
for validation workflows (spec §5 c.12, §8 ring 3). Returns a CrosscheckReport;
never silently disagrees -- .assert_agree() raises on mismatch."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QpaUnavailableError
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
    gap = session.run(scripts.hochschild_dims_script(algebra, top) + " hh;")
    qpa = _read_int_list(gap)
    return CrosscheckReport("hochschild", list(ours), qpa, list(ours) == qpa)


def crosscheck_module_ext(algebra, M, top: int) -> CrosscheckReport:
    """Self-Ext Ext^*(M, M) vs QPA (via ExtAlgebraGenerators). Distinct-module
    Ext(M, N) is a flagged post-v1 extension (needs ExtOverAlgebra + syzygies)."""
    session.require_gap()
    ours = [algebra.ext(M, M, n).dimension() for n in range(top + 1)]  # Plan 05 surface
    gap = session.run(
        scripts.module_self_ext_dims_script(algebra, M.dimension_vector(), top) + " ext;")
    qpa = _read_int_list(gap)
    return CrosscheckReport("module_ext", list(ours), qpa, list(ours) == qpa)


def crosscheck(algebra, what: str, *args, **kwargs) -> CrosscheckReport:
    """Dispatch. what="hochschild" -> crosscheck_hochschild(algebra, top);
    what="module_ext" -> crosscheck_module_ext(algebra, M, top) (self-Ext)."""
    if what == "hochschild":
        return crosscheck_hochschild(algebra, *args, **kwargs)
    if what == "module_ext":
        return crosscheck_module_ext(algebra, *args, **kwargs)
    raise QpaUnavailableError(f"unknown cross-check {what!r}",
                              hint='use "hochschild" or "module_ext"')
