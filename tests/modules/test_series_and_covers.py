import pytest

from quiverlab import GF, Quiver

pytestmark = pytest.mark.oracle_selfcert


def _a3():
    return Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)}).algebra(
        relations=["a*b"], field=GF(5))


def test_projective_cover_hom_is_epi_with_superfluous_kernel():
    A = _a3()
    S1 = A.simple(1)
    f = S1.projective_cover_hom()
    assert f.is_epi() and f.tgt.dim == S1.dim
    K, _ = f.kernel()
    # cover kernel lies in rad(P): K + rad P = rad P
    assert K.dimension_vector() == f.src.radical().dimension_vector()


def test_injective_envelope_hom_is_mono():
    A = _a3()
    S2 = A.simple(2)
    f = S2.injective_envelope_hom()
    assert f.is_mono() and f.src.dim == S2.dim
    assert f.tgt.socle().dimension_vector() == S2.dimension_vector()


def test_radical_series_strictly_decreases_to_zero():
    A = _a3()
    P1 = A.projective(1)
    series = P1.radical_series()
    dims = [X.dim for X in series]
    assert dims[0] == P1.dim and dims[-1] == 0
    assert all(a > b for a, b in zip(dims, dims[1:]))
    assert len(series) - 1 <= A.loewy_length()


def test_socle_series_reaches_the_module():
    A = _a3()
    P1 = A.projective(1)
    socs = P1.socle_series()
    assert socs[0].dim == 0 and socs[-1].dim == P1.dim


def test_loewy_layers_sum_to_composition_factors():
    A = _a3()
    P1 = A.projective(1)
    layers = P1.loewy_layers()
    total = {}
    for layer in layers:
        for k, v in layer.items():
            total[k] = total.get(k, 0) + v
    assert total == P1.composition_factors()
    # all simples here are 1-dimensional, so factor count == module dim
    assert sum(total.values()) == P1.dim
