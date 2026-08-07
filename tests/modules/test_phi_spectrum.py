"""phi-spectrum + gaps (Plan 53 / R23b; Barrios-Mata-Rama 1810.12112). Self-cert:
0 and phidim are attained; when 0 < phidim < infinity, 1 and phidim-1 are attained
(BMR theorem). Rep-finite only; a partial spectrum claims no gaps. Over QQ."""
import pytest

from quiverlab import TruncatedPathAlgebra, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.homdims import phi_spectrum

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


@selfcert
def test_bmr_endpoints_attained():
    A = TruncatedPathAlgebra("A4", 2, field=QQ)          # phidim = 3
    S = phi_spectrum(A)
    assert S.complete
    assert 0 in S.values and S.phidim in S.values
    if 0 < S.phidim < 10**6:
        assert 1 in S.values and (S.phidim - 1) in S.values   # BMR


@lit
def test_hereditary_spectrum_is_zero_one():
    A = linear_path_algebra(3, field=QQ)                 # hereditary: phi = pd in {0,1}
    S = phi_spectrum(A)
    assert S.values == [0, 1] and S.gaps == [] and S.phidim == 1
