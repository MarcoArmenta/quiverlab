"""Standard/costandard modules Delta(i)/Nabla(i) (Plan 47). Literature (Dlab-Ringel): kA_n
natural order Delta(i)=S_i, Nabla(i)=[1..i]; opposite order Delta(i)=P(i). Self-cert: top
Delta(i)=S_i, [Delta(i):S(i)]=1, socle Nabla(i)=S_i. Char-clean (a GF(2) cell)."""
import pytest

from quiverlab import GF, linear_path_algebra
from quiverlab.fields import QQ
from quiverlab.modules.hom import is_isomorphic
from quiverlab.modules.quasihereditary import costandard_modules, standard_modules

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def _a3(field=QQ):
    return linear_path_algebra(3, field=field)


@lit
def test_natural_order_standards_are_simple():
    A = _a3()
    D = standard_modules(A)                       # order=None -> natural
    for v in (1, 2, 3):
        assert D[v].dim == 1                       # Delta(v) = S_v
        assert is_isomorphic(D[v], A.simple(v))


@lit
def test_natural_order_costandards_are_intervals():
    A = _a3()
    N = costandard_modules(A)
    for v in (1, 2, 3):
        assert N[v].dimension_vector() == {w: (1 if w <= v else 0) for w in (1, 2, 3)}
        assert is_isomorphic(N[v], A.injective(v))  # Nabla(v) = I(v) = [1..v]


@lit
def test_opposite_order_standards_are_projectives():
    A = _a3()
    order = [3, 2, 1]                              # n < ... < 1
    D = standard_modules(A, order)
    for v in (1, 2, 3):
        assert is_isomorphic(D[v], A.projective(v))  # Delta(v) = P(v) = [v..3]


@selfcert
def test_delta_top_is_simple_and_mult_one():
    A = _a3()
    for order in (None, [3, 2, 1], [2, 1, 3]):
        D = standard_modules(A, order)
        for v in (1, 2, 3):
            assert D[v].top().dimension_vector() == {w: (1 if w == v else 0)
                                                     for w in (1, 2, 3)}   # top Delta(i)=S_i
            cf = D[v].composition_factors()
            assert cf.get(str(v), 0) == 1          # [Delta(i):S(i)] = 1


@selfcert
def test_costandard_socle_is_simple():
    A = _a3()
    N = costandard_modules(A)
    for v in (1, 2, 3):
        soc = N[v].socle()
        assert soc.dimension_vector() == {w: (1 if w == v else 0) for w in (1, 2, 3)}


@selfcert
def test_standard_char_clean_gf2():
    A = _a3(GF(2))
    D = standard_modules(A)
    assert all(D[v].dim == 1 for v in (1, 2, 3))   # no char caveat on Delta


@selfcert
def test_bad_order_refused():
    from quiverlab.errors import QuiverlabError
    A = _a3()
    with pytest.raises(QuiverlabError):
        standard_modules(A, [1, 2])                # not a permutation of the vertices
