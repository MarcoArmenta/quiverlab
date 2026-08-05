"""Ideal triangulations of marked surfaces (Plan 48). A triangulation is a tuple of
ordered triangles, each an ANTICLOCKWISE 3-tuple of SIDES (interior-arc ints + boundary-
segment strings). Construction validates the FST arc-adjacency: interior arcs in exactly
2 triangles, boundary segments in exactly 1, arc count == surface.arc_count(). Canonical
constructors for the disc fan, the annulus band (built by cutting along one arc into a
strip and re-gluing -- gives an ACYCLIC quiver), the hexagon-with-internal-triangle, and
the once-punctured torus. Float-free (all combinatorics are integers/strings)."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from quiverlab.errors import QuiverlabError
from quiverlab.surfaces.marked import MarkedSurface


def _is_arc(side) -> bool:
    """True iff ``side`` is an interior arc id (a positive int, not a bool, not a
    boundary-segment string)."""
    return isinstance(side, int) and not isinstance(side, bool)


@dataclass(frozen=True)
class Triangulation:
    """An ideal triangulation: ``surface`` + ``triangles`` (a tuple of anticlockwise
    3-tuples of sides). Construction self-certifies the FST arc-adjacency and arc count."""

    surface: MarkedSurface
    triangles: tuple

    def __post_init__(self):
        object.__setattr__(self, "triangles",
                           tuple(tuple(t) for t in self.triangles))
        for tri in self.triangles:
            if len(tri) != 3:
                raise QuiverlabError("Triangulation: every triangle needs 3 sides",
                                     hint=f"got {tri!r}")
            if len(set(tri)) != 3:                        # a side repeated => self-folded
                raise QuiverlabError(
                    "Triangulation: a side appears twice in one triangle (self-folded)",
                    hint="self-folded triangles need a puncture and are out of v1 scope "
                         "(successor P48.1)")
        occ = Counter(s for tri in self.triangles for s in tri)
        for side, k in occ.items():
            want = 2 if _is_arc(side) else 1
            if k != want:
                raise QuiverlabError(
                    f"Triangulation: side {side!r} borders {k} triangle(s), expected {want}",
                    hint="interior arcs must border exactly 2 triangles, boundary "
                         "segments exactly 1")
        n_arcs = len(self.arcs())
        if n_arcs != self.surface.arc_count():
            raise QuiverlabError(
                f"Triangulation: {n_arcs} arcs but surface has {self.surface.arc_count()}",
                hint="the arc-count self-certificate failed (bad triangulation)")
        if len(self.triangles) != self.surface.triangle_count():
            raise QuiverlabError(
                f"Triangulation: {len(self.triangles)} triangles but surface has "
                f"{self.surface.triangle_count()}", hint="arc-count self-cert failed")

    def arcs(self) -> tuple:
        """The interior-arc ids, sorted ascending."""
        return tuple(sorted({s for tri in self.triangles for s in tri if _is_arc(s)}))

    def boundary_segments(self) -> tuple:
        """The boundary-segment ids, sorted by string."""
        return tuple(sorted({s for tri in self.triangles for s in tri if not _is_arc(s)},
                            key=str))

    def triangles_containing(self, arc) -> tuple:
        """The (<= 2) triangles having ``arc`` as a side."""
        return tuple(tri for tri in self.triangles if arc in tri)


def fan_triangulation(marked) -> Triangulation:
    """The fan of the ``marked``-gon from vertex 0 (vertices 0..marked-1 anticlockwise).
    Diagonal 0-j (j = 2..marked-2) is arc id ``j-1`` (arcs 1..marked-3). Triangle
    ``(0, i, i+1)`` has anticlockwise sides ``(side 0-i, side i-(i+1), side (i+1)-0)``.
    The disc arbiter (qp.quiver_of) pins the angle orientation so this gives linear
    ``A_{marked-3}``: 1 -> 2 -> ... -> (marked-3)."""
    if not _is_arc(marked) or marked < 4:
        raise QuiverlabError("fan_triangulation: need marked >= 4 (a polygon)",
                             hint=f"got marked={marked!r}; a triangle (3) has no diagonal")
    S = MarkedSurface(genus=0, boundary_marked=(marked,), punctures=0)

    def diag(j):                                          # arc id of diagonal 0-j
        return j - 1

    def bseg(i):                                          # boundary segment i-(i+1)
        return f"b0_{i}"

    tris = []
    for i in range(1, marked - 1):                        # triangle (0, i, i+1)
        left = bseg(0) if i == 1 else diag(i)             # side 0-i (bd 0-1 for i==1)
        base = bseg(i)                                    # side i-(i+1)
        right = bseg(marked - 1) if i == marked - 2 else diag(i + 1)  # side (i+1)-0
        tris.append((left, base, right))                  # anticlockwise 0->i->(i+1)->0
    return Triangulation(S, tuple(tris))


def annulus_triangulation(n, m) -> Triangulation:
    """A standard band triangulation of the annulus C(n, m) (n outer marked points, m
    inner). Built by cutting the annulus along one bridging arc into a STRIP between the
    outer chain (n+1 corners) and the inner chain (m+1 corners), triangulating the strip
    by an integer balanced merge (advance the outer side while ``jt*m <= jb*n``), then
    re-gluing the two copies of the bridging arc (arc id 1). The resulting quiver is the
    ACYCLICALLY-oriented affine A~_{n+m-1} (arc count n+m, no oriented cycle -- so the P44
    Jacobian is finite; see qp.jacobian_of). Anticlockwise side orders: an outer-boundary
    triangle is ``(cur_arc, new_arc, outer_seg)``, an inner-boundary triangle is
    ``(cur_arc, inner_seg, new_arc)``."""
    if not (_is_arc(n) and _is_arc(m)) or n < 1 or m < 1:
        raise QuiverlabError("annulus_triangulation: need n >= 1 and m >= 1",
                             hint=f"got n={n!r}, m={m!r}")
    S = MarkedSurface(genus=0, boundary_marked=(n, m), punctures=0)
    # strip merge over states (jt, jb): jt counts outer edges (0..n), jb inner (0..m).
    steps = []                                            # (typ, cur, new, bnd)
    jt = jb = 0
    for _ in range(n + m):
        if jt >= n:
            adv_top = False
        elif jb >= m:
            adv_top = True
        else:
            adv_top = jt * m <= jb * n                    # balanced; tie -> outer
        cur = (jt, jb)
        if adv_top:
            new = (jt + 1, jb)
            steps.append(("top", cur, new, f"b0_{jt}"))   # outer boundary comp 0
            jt += 1
        else:
            new = (jt, jb + 1)
            steps.append(("bot", cur, new, f"b1_{jb}"))   # inner boundary comp 1
            jb += 1

    def key(st):                                          # (n,m) is the same arc as (0,0)
        return (0, 0) if st == (n, m) else st

    arc_id = {}
    for _typ, cur, new, _bnd in steps:
        for st in (cur, new):
            k = key(st)
            if k not in arc_id:
                arc_id[k] = len(arc_id) + 1
    tris = []
    for typ, cur, new, bnd in steps:
        c, nw = arc_id[key(cur)], arc_id[key(new)]
        if typ == "top":
            tris.append((c, nw, bnd))                     # (cur_arc, new_arc, outer_seg)
        else:
            tris.append((c, bnd, nw))                     # (cur_arc, inner_seg, new_arc)
    return Triangulation(S, tuple(tris))


def hexagon_with_internal_triangle() -> Triangulation:
    """The hexagon (6 boundary marked points 0..5) triangulated by the central triangle
    on the diagonals {1: 1-3, 2: 3-5, 3: 5-1} plus three boundary "ear" triangles. The
    one internal triangle (all sides arcs) yields a single 3-cycle potential, so
    jacobian_of is the triangle algebra of dim 6 (P44's pin). Anticlockwise orders:
    central = (1, 2, 3); ears (corners 1-2-3, 3-4-5, 5-0-1)."""
    S = MarkedSurface(genus=0, boundary_marked=(6,), punctures=0)
    tris = (
        (1, 2, 3),                       # central internal triangle (corners 1,3,5)
        ("b0_1", "b0_2", 1),             # ear at corner 2 (triangle 1-2-3)
        ("b0_3", "b0_4", 2),             # ear at corner 4 (triangle 3-4-5)
        ("b0_5", "b0_0", 3),             # ear at corner 0 (triangle 5-0-1)
    )
    return Triangulation(S, tris)


def once_punctured_torus() -> Triangulation:
    """The once-punctured torus (g=1, p=1): 3 loop-arcs at the single puncture, 2
    triangles each using all three arcs. FST-admissible; OUT of v1 scope -- quiver_of /
    jacobian_of refuse it (it is the pinned loud-refusal oracle, Task 3)."""
    S = MarkedSurface(genus=1, boundary_marked=(), punctures=1)
    tris = ((1, 2, 3), (1, 3, 2))
    return Triangulation(S, tris)
