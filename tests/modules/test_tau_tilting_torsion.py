"""Torsion lattice + bricks/semibricks (Plan 45 / C4). Literature: the AIR four-way count
identity #stau-tilt = #f.f. torsion classes = #semibricks = Catalan(n+1) on kA2/kA3 (the
silting leg is added in Task G). Self-cert: (A,0) is the unique Hasse source, (0,A) the
unique sink; every brick label has end_dim 1; semibricks are pairwise Hom-orthogonal."""
import pytest

from quiverlab import linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.hom import end_dim, hom_dim
from quiverlab.tautilting.mutation import exchange_graph
from quiverlab.tautilting.torsion import (bricks, hasse_orientation, semibricks,
                                          torsion_class_data)

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert
xeng = pytest.mark.oracle_crossengine


@selfcert
def test_hasse_has_unique_source_and_sink():
    A = linear_path_algebra(3, field=QQ)
    eg = exchange_graph(A)
    o = hasse_orientation(eg)
    indeg = {i: 0 for i in range(len(eg.vertices))}
    outdeg = dict(indeg)
    for (i, j), d in o.items():
        a, b = (i, j) if d == "down" else (j, i)     # a -> b downward
        outdeg[a] += 1
        indeg[b] += 1
    sources = [i for i in indeg if indeg[i] == 0]
    sinks = [i for i in outdeg if outdeg[i] == 0]
    assert len(sources) == 1 and len(sinks) == 1
    assert eg.vertices[sources[0]]["is_initial"]      # (A,0) on top


@selfcert
def test_brick_labels_are_bricks():
    A = linear_path_algebra(3, field=QQ)
    for B in bricks(A):
        assert end_dim(B) == 1                         # End(B) = k over QQ (alg closed base)


@selfcert
def test_semibricks_are_hom_orthogonal():
    A = linear_path_algebra(2, field=QQ)
    for sb in semibricks(A):
        mods = list(sb)
        for i in range(len(mods)):
            for j in range(len(mods)):
                if i != j:
                    assert hom_dim(mods[i], mods[j]) == 0    # pairwise Hom-orthogonal


@lit
@pytest.mark.parametrize("n, count", [(2, 5), (3, 14)])
def test_air_four_way_identity_torsion_and_semibricks(n, count):
    # #support tau-tilting pairs == #f.f. torsion classes == #semibricks (silting leg: G).
    A = linear_path_algebra(n, field=QQ)
    eg = exchange_graph(A)
    n_pairs = len(eg.vertices)
    n_torsion = len({tuple(torsion_class_data(v["pair"])["gen_dimvecs"])
                     for v in eg.vertices})
    n_semibricks = len(semibricks(A))
    assert n_pairs == n_torsion == n_semibricks == count


@xeng
def test_torsion_classes_are_distinct_per_pair():
    # the pair -> Gen(M) map is injective (AIR bijection): fingerprints all distinct.
    A = linear_path_algebra(3, field=QQ)
    eg = exchange_graph(A)
    fps = [tuple(torsion_class_data(v["pair"])["gen_dimvecs"]) for v in eg.vertices]
    assert len(set(fps)) == len(fps)
