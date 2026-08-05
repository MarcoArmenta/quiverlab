"""MarkedSurface validation + the FST admissibility exclusion list (Plan 48, deep).
Loud refusals: the monogon/digon/triangle (n<=0), the small punctured spheres, the
once-punctured monogon; bad marked-point counts."""
import pytest

from quiverlab.errors import QuiverlabError
from quiverlab.surfaces.marked import MarkedSurface

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@lit
@pytest.mark.parametrize("g, bm, p", [
    (0, (1,), 0),        # unpunctured monogon (n = -2)
    (0, (1,), 1),        # once-punctured monogon (n = 1 but FST-excluded)
    (0, (2,), 0),        # unpunctured digon (n = -1)
    (0, (3,), 0),        # unpunctured triangle (n = 0)
    (0, (), 1),          # sphere, 1 puncture
    (0, (), 2),          # sphere, 2 punctures
    (0, (), 3),          # sphere, 3 punctures (n = 3 but FST-excluded)
])
def test_fst_inadmissible_surfaces_refused(g, bm, p):
    with pytest.raises(QuiverlabError):
        MarkedSurface(genus=g, boundary_marked=bm, punctures=p)


@selfcert
@pytest.mark.parametrize("g, bm, p", [(-1, (4,), 0), (0, (0,), 0), (0, (2, 0), 0),
                                      (0, (4,), -1)])
def test_bad_counts_refused(g, bm, p):
    with pytest.raises(QuiverlabError):
        MarkedSurface(genus=g, boundary_marked=bm, punctures=p)


@lit
def test_admissible_surfaces_construct():
    # square (n=1), pentagon (n=2), annulus (n+m), 4-punctured sphere, once-punctured torus
    assert MarkedSurface(0, (4,), 0).arc_count() == 1
    assert MarkedSurface(0, (2, 2), 0).arc_count() == 4
    assert MarkedSurface(0, (), 4).arc_count() == 6        # 4-punctured sphere is admissible
    assert MarkedSurface(1, (), 1).arc_count() == 3
