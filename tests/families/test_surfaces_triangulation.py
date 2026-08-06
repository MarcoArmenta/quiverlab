"""Ideal triangulations (Plan 48, deep). Self-cert: arc-adjacency counts + the
arc-count identity on every constructor; loud on a self-folded configuration and on a
triangulation whose arc count contradicts the surface."""
import pytest

from quiverlab.errors import QuiverlabError
from quiverlab.surfaces.marked import MarkedSurface
from quiverlab.surfaces.triangulation import (Triangulation, annulus_triangulation,
                                             fan_triangulation, once_punctured_torus)

selfcert = pytest.mark.oracle_selfcert


@selfcert
@pytest.mark.parametrize("marked", [4, 5, 6, 7, 8])
def test_fan_adjacency_and_arc_count(marked):
    T = fan_triangulation(marked)
    assert len(T.arcs()) == T.surface.arc_count() == marked - 3
    assert len(T.triangles) == T.surface.triangle_count() == marked - 2
    for a in T.arcs():                                   # every interior arc in 2 triangles
        assert len(T.triangles_containing(a)) == 2
    for seg in T.boundary_segments():                    # every boundary segment in 1
        assert sum(seg in tri for tri in T.triangles) == 1


@selfcert
@pytest.mark.parametrize("n, m", [(2, 1), (2, 2), (3, 1)])
def test_annulus_adjacency_and_arc_count(n, m):
    T = annulus_triangulation(n, m)
    assert len(T.arcs()) == T.surface.arc_count() == n + m
    for a in T.arcs():
        assert len(T.triangles_containing(a)) == 2


@selfcert
def test_wrong_arc_count_refused():
    S = MarkedSurface(0, (5,), 0)                        # pentagon: needs exactly 2 arcs
    with pytest.raises(QuiverlabError):
        Triangulation(S, triangles=((1, "b0_0", "b0_1"),))   # one triangle, one arc: wrong


@selfcert
def test_self_folded_refused():
    S = MarkedSurface(0, (4,), 0)
    with pytest.raises(QuiverlabError):
        Triangulation(S, triangles=((1, 1, "b0_0"),))        # arc 1 twice in one triangle


def test_boundary_topology_guard_rejects_disc_shaped_list_on_annulus():
    # Devil's-advocate fold (2026-08-05): count checks alone accepted a
    # disc-shaped triangle list on a 2-boundary surface; the per-component
    # boundary-segment guard must refuse it loudly.
    import pytest
    from quiverlab.errors import QuiverlabError
    from quiverlab.surfaces import MarkedSurface, Triangulation

    S = MarkedSurface(0, (2, 1), 0)                # annulus, arcs = 3, tris = 3
    # a triangle list using ONLY b0_* segments (disc-like topology)
    tris = ((1, 2, "b0_0"), (2, 3, "b0_1"), (3, 1, "b0_2"))
    with pytest.raises(QuiverlabError, match="boundary"):
        Triangulation(S, tris)
