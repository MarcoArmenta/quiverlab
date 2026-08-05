"""Marked bordered surfaces (Plan 48). A MarkedSurface is the combinatorial datum
(genus, boundary marked-point counts, punctures) of a compact oriented surface with
marked points. The ideal-arc count n = 6g-6+3(b+p)+Sum k_i (Fomin-Shapiro-Thurston 2008)
and the Euler invariants are derived exactly (derivation in the Plan-48 doc, Task 1);
construction refuses the FST-inadmissible surfaces loudly. Float-free (all counts are
ints)."""
from __future__ import annotations

from dataclasses import dataclass

from quiverlab.errors import QuiverlabError

# FST 2008 explicit exclusions with n > 0 (the n <= 0 degeneracies -- monogon/digon/
# triangle -- are refused by the arc-count check itself):
#   (genus, sorted-boundary-tuple, punctures)
_FST_EXCLUDED = {
    (0, (), 1), (0, (), 2), (0, (), 3),      # sphere with 1/2/3 punctures
    (0, (1,), 1),                            # once-punctured monogon (n = 1 but excluded)
}


@dataclass(frozen=True)
class MarkedSurface:
    """A compact connected oriented surface (S, M): genus ``genus``, boundary components
    carrying ``boundary_marked = (k1, ..., kb)`` marked points each (k_i >= 1), and
    ``punctures`` interior marked points. Refuses the FST-inadmissible surfaces on
    construction. See Fomin-Shapiro-Thurston 2008 (``fomin_shapiro_thurston``)."""

    genus: int
    boundary_marked: tuple
    punctures: int = 0

    def __post_init__(self):
        object.__setattr__(self, "boundary_marked", tuple(self.boundary_marked))
        if not isinstance(self.genus, int) or isinstance(self.genus, bool):
            raise QuiverlabError("MarkedSurface: genus must be an int",
                                 hint=f"got genus={self.genus!r}")
        if not isinstance(self.punctures, int) or isinstance(self.punctures, bool):
            raise QuiverlabError("MarkedSurface: punctures must be an int",
                                 hint=f"got punctures={self.punctures!r}")
        if self.genus < 0:
            raise QuiverlabError("MarkedSurface: genus must be >= 0",
                                 hint=f"got genus={self.genus}")
        if self.punctures < 0:
            raise QuiverlabError("MarkedSurface: punctures must be >= 0",
                                 hint=f"got punctures={self.punctures}")
        for k in self.boundary_marked:
            # bool is an int subclass -- reject True/False smuggled in as a "1"/"0" count.
            if not isinstance(k, int) or isinstance(k, bool) or k < 1:
                raise QuiverlabError(
                    "MarkedSurface: each boundary component needs >= 1 marked point",
                    hint=f"got boundary_marked={self.boundary_marked}")
        b, p, c = self.num_boundary_components, self.punctures, sum(self.boundary_marked)
        if b == 0 and p == 0:
            raise QuiverlabError("MarkedSurface: no marked points -- cannot triangulate",
                                 hint="a closed surface needs >= 1 puncture")
        key = (self.genus, tuple(sorted(self.boundary_marked)), self.punctures)
        n = 6 * self.genus - 6 + 3 * (b + p) + c
        if n < 1 or key in _FST_EXCLUDED:
            raise QuiverlabError(
                "MarkedSurface: not FST-admissible (no nontrivial ideal triangulation)",
                hint="excluded (FST 2008): the monogon/digon/triangle, spheres with "
                     "<=3 punctures, and the once-punctured monogon")

    @property
    def num_boundary_components(self) -> int:
        """b = the number of boundary components."""
        return len(self.boundary_marked)

    @property
    def total_marked(self) -> int:
        """p + Sum k_i = all 0-cells (punctures + boundary marked points)."""
        return self.punctures + sum(self.boundary_marked)

    def euler_characteristic(self) -> int:
        """chi(S) = 2 - 2g - b for a compact surface with b boundary components."""
        return 2 - 2 * self.genus - self.num_boundary_components

    def arc_count(self) -> int:
        """n = 6g - 6 + 3(b + p) + Sum k_i -- the number of interior arcs in any ideal
        triangulation (FST 2008; boundary segments are NOT arcs)."""
        b, p, c = self.num_boundary_components, self.punctures, sum(self.boundary_marked)
        return 6 * self.genus - 6 + 3 * (b + p) + c

    def triangle_count(self) -> int:
        """t = (2n + c)/3 -- an exact integer by the side-counting identity 3t = 2n + c."""
        n, c = self.arc_count(), sum(self.boundary_marked)
        assert (2 * n + c) % 3 == 0                       # exact by the derivation
        return (2 * n + c) // 3

    def in_v1_scope(self) -> bool:
        """True iff the surface is in Plan-48 v1 scope: unpunctured with non-empty
        boundary (the ABCP/LFS gentle regime)."""
        return self.punctures == 0 and self.num_boundary_components >= 1
