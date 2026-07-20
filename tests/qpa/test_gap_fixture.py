"""The kA_2/GF(2) enveloping-algebra HH fixture: quiverlab vs QPA both give
[1,0,0] (spec §8 ring 3). Skips locally (no GAP); mandatory in the qpa.yml CI job.

The oracle HH^*(kA_2) = [1, 0, 0] is Happel-certain: kA_2 = path algebra of
1 --a--> 2 (no relations) is hereditary with a tree quiver and no oriented cycles,
so HH^0 = k (center = scalars) and HH^{>=1} = 0; over GF(2) the dimensions are
identical (characteristic-independent here). The QPA enveloping-algebra route
(HH^n(A) = Ext^n_{A^e}(A, A)) must reproduce [1, 0, 0] -- VERIFY AT EXECUTION.
"""
import pytest

from quiverlab import Quiver, GF
from quiverlab.qpa import session, scripts


def _ka2():
    return Quiver([1, 2], {"a": (1, 2)}).algebra(relations=[], field=GF(2))


@pytest.mark.skipif(session.should_skip_qpa(), reason="[qpa] backend not installed")
def test_ka2_gf2_hochschild_matches_qpa():
    A = _ka2()
    assert A.hochschild_cohomology(2).dims == [1, 0, 0]     # Happel; char-independent
    report = A.crosscheck("hochschild", 2)
    report.assert_agree()
    assert report.qpa == [1, 0, 0]                          # VERIFY AT EXECUTION


@pytest.mark.fast
def test_ka2_gf2_script_is_wellformed_even_without_gap():
    """The emitted GAP source is stable and mentions the enveloping route -- this
    part runs everywhere (no GAP needed), pinning the script the CI job will run."""
    src = scripts.hochschild_dims_script(_ka2(), 2)
    assert 'Quiver(2, [[1, 2, "a"]])' in src
    assert "PathAlgebra(GF(2), Q)" in src
    assert "EnvelopingAlgebra(A)" in src
