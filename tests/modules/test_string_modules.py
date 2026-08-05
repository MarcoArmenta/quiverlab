

def test_band_module_refuses_proper_power_and_non_band():
    # Devil's-advocate MEDIUM (2026-08-05): b^2 previously materialised as a
    # silently DECOMPOSABLE module; both bad inputs must refuse loudly now.
    import pytest
    from quiverlab import GF, Quiver
    from quiverlab.errors import QuiverlabError
    from quiverlab.strings.modules import band_module
    from quiverlab.strings.walks import find_bands

    K = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(field=GF(5))
    (band,) = [w for w in find_bands(K) if len(w) == 2][:1]
    band_module(K, band, "2", mult=1)                 # primitive: fine
    with pytest.raises(QuiverlabError, match="proper power|PRIMITIVE"):
        band_module(K, band + band, "2", mult=1)
    with pytest.raises(QuiverlabError, match="not a band"):
        band_module(K, (band[0],), "2", mult=1)       # single letter: no closure


def test_band_module_mult_two_is_indecomposable():
    # Devil's-advocate coverage (2026-08-05): the Jordan-block claim at mult=2.
    from quiverlab import GF, Quiver
    from quiverlab.strings.modules import band_module
    from quiverlab.strings.walks import find_bands

    K = Quiver([1, 2], {"a": (1, 2), "b": (1, 2)}).algebra(field=GF(32003))
    (band,) = [w for w in find_bands(K) if len(w) == 2][:1]
    M = band_module(K, band, "2", mult=2)
    assert M.dim == 4
    assert M.is_indecomposable()
