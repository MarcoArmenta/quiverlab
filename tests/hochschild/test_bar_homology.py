from quiverlab import CC, GF, Quiver


def _dual(field):
    Q = Quiver(vertices=[1], arrows={"x": (1, 1)})
    return Q.algebra(relations=["x^2"], field=field)


def test_dual_numbers_char0_homology():
    # HH_0 = 2, HH_n = 1 for n >= 1 in char 0
    assert _dual(CC).hochschild_homology(4).dims == [2, 1, 1, 1, 1]


def test_dual_numbers_char2_homology():
    assert _dual(GF(2)).hochschild_homology(4).dims == [2, 2, 2, 2, 2]


def test_kA2_homology_vanishes_positively():
    Q = Quiver(vertices=[1, 2], arrows={"a": (1, 2)})
    A = Q.algebra()
    # HH_0 = k^{#vertices} for an acyclic monomial algebra; higher vanish (hereditary, acyclic)
    assert A.hochschild_homology(3).dims == [2, 0, 0, 0]


def test_symmetric_algebra_duality_smoke():
    # k[x]/(x^2) is symmetric: HH^n and HH_n dimensions agree in every characteristic
    for field in (CC, GF(2), GF(3)):
        A = _dual(field)
        assert A.hochschild_cohomology(3).dims == A.hochschild_homology(3).dims
