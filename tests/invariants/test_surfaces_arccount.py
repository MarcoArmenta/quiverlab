"""FST ideal-arc count + Euler invariants (Plan 48). Pure combinatorics -- fast bucket.
Self-cert: n = 6g-6+3(b+p)+Sum k_i and t=(2n+c)/3 are integers matching the derivation;
literature: disc(n+3)->A_n arc count, annulus(n,m)->n+m, once-punctured torus->3."""
import pytest

from quiverlab.surfaces.marked import MarkedSurface

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature


@lit
@pytest.mark.parametrize("marked, n", [(4, 1), (5, 2), (6, 3), (7, 4), (10, 7)])
def test_disc_arc_count_is_A_n(marked, n):
    # disc with `marked` boundary points -> A_{marked-3}: n = marked - 3 arcs.
    S = MarkedSurface(genus=0, boundary_marked=(marked,), punctures=0)
    assert S.arc_count() == n == marked - 3
    assert S.triangle_count() == n + 1          # (marked-2) triangles
    assert S.euler_characteristic() == 1        # 2 - 0 - 1


@lit
@pytest.mark.parametrize("n, m", [(1, 1), (2, 1), (2, 2), (3, 1)])
def test_annulus_arc_count_is_n_plus_m(n, m):
    S = MarkedSurface(genus=0, boundary_marked=(n, m), punctures=0)
    assert S.arc_count() == n + m
    assert S.triangle_count() == n + m          # (2(n+m)+(n+m))/3
    assert S.euler_characteristic() == 0        # 2 - 0 - 2


@lit
def test_once_punctured_torus_has_three_arcs():
    S = MarkedSurface(genus=1, boundary_marked=(), punctures=1)
    assert S.arc_count() == 3 and S.triangle_count() == 2
    assert not S.in_v1_scope()                  # closed + punctured -> out of v1 scope


@selfcert
@pytest.mark.parametrize("g, bm, p", [(0, (5,), 0), (0, (2, 2), 0), (1, (1,), 0),
                                      (0, (4, 4), 0), (2, (1,), 0)])
def test_side_counting_identity_holds(g, bm, p):
    S = MarkedSurface(genus=g, boundary_marked=bm, punctures=p)
    n, t, c = S.arc_count(), S.triangle_count(), sum(bm)
    assert 3 * t == 2 * n + c                    # identity (II)
    assert p - n + t == S.euler_characteristic() # identity (I)
