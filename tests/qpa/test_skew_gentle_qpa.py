"""QPA as the oracle for skew-gentle algebras (Plan 68). QPA has NO skew-gentle
constructor / recognizer -- so the crosschecks are on the SPLIT algebra as a plain
kQ/I: dim, Cartan, selfinjectivity (Chen Cor 1.2c), and module decompose. qpa-marked."""
import pytest

from quiverlab import Quiver
from quiverlab.fields import QQ
from quiverlab.qpa import session
from quiverlab.skewgentle.split import SkewGentleAlgebra
from quiverlab.skewgentle.triple import SkewGentleTriple

pytestmark = pytest.mark.skipif(session.should_skip_qpa(),
                                reason="[qpa] backend not installed")


def _t():
    return SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}), relations=[], special={2})


def test_qpa_has_no_skew_gentle_surface():
    lg = session.libgap_handle()
    for name in ("SkewGentleAlgebra", "IsSkewGentleAlgebra", "SkewGentleTriple"):
        assert not bool(lg.eval(f'IsBoundGlobal("{name}")')), \
            f"QPA now ships {name} -- add a real crosscheck (honest-scope changed)"


def test_split_algebra_dim_vs_qpa():
    A = SkewGentleAlgebra(_t(), field=QQ)
    A.crosscheck("dim").assert_agree()                   # split is a plain kQ/I


def test_mesh_split_dim_vs_qpa():
    # the non-monomial mesh split (dim 9) crosschecks as an ordinary presented kQ/I.
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)})
    A = SkewGentleAlgebra(quiver=Q, relations=["a*b"], special={2}, field=QQ)
    A.crosscheck("dim").assert_agree()


def test_selfinjective_crosscheck():
    # Chen Cor 1.2(c): the 2-cycle gentle (Sp = empty) is selfinjective; QPA agrees.
    Q2 = Quiver([1, 2], {"a": (1, 2), "b": (2, 1)})
    A = SkewGentleAlgebra(quiver=Q2, relations=["a*b", "b*a"], special=(), field=QQ)
    A.crosscheck("is_selfinjective").assert_agree()      # IsSelfinjectiveAlgebra
