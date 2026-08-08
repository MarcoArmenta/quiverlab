"""The idempotent-split skew-gentle algebra (Plan 68 / Chen 2212.06467 sec 3).
Self-cert: dim(split) == dim(associated gentle) (HZZ Lemma 1.5); n_split == |Q_0|+|Sp|
(HZZ Rmk 1.4 / K_0). Literature: Sp = empty BYTE-REDUCES to the plain gentle algebra;
Chen Cor 1.2(c) an INDECOMPOSABLE skew-gentle algebra is selfinjective iff Sp = empty and
the gentle A(Q,I) is selfinjective (the tests use connected algebras only). Loud: a broken
triple; a dimension-certificate failure."""
import pytest

from quiverlab import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import GF, QQ
from quiverlab.skewgentle.split import SkewGentleAlgebra, split_quiver
from quiverlab.skewgentle.triple import SkewGentleTriple, associated_gentle

selfcert = pytest.mark.oracle_selfcert
lit = pytest.mark.oracle_literature
xeng = pytest.mark.oracle_crossengine


def _t_arrow(special={2}):            # Q: 1 --a--> 2 , Sp = {2}
    return SkewGentleTriple.make(Quiver([1, 2], {"a": (1, 2)}),
                                 relations=[], special=special)


@selfcert
def test_dim_law_and_rank():
    t = _t_arrow()
    A = SkewGentleAlgebra(t, field=QQ)
    assert A.dim == associated_gentle(t, field=QQ).dim == 5      # HZZ Lemma 1.5
    Qhat, _rels, _meta = split_quiver(t)
    assert len(list(Qhat.vertices)) == 2 + 1                     # |Q_0| + |Sp| = 3
    assert A.dim == 5 and len(list(A.quiver.vertices)) == 3


@selfcert
@pytest.mark.parametrize("field", [QQ, GF(32003), GF(2)])       # GF(2): char-free!
def test_split_is_characteristic_free(field):
    A = SkewGentleAlgebra(_t_arrow(), field=field)
    assert A.dim == 5                                            # even in char 2


# --- the MESH oracle: a special vertex in the MIDDLE of a length-2 relation (H2a) ------
def _t_mesh():
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)})
    return SkewGentleTriple.make(Q, relations=["a*b"], special={2})   # a*b IN I


@selfcert
def test_mesh_dim_law_is_necessary_not_sufficient():
    t = _t_mesh()
    A = SkewGentleAlgebra(t, field=QQ)
    assert A.dim == associated_gentle(t, field=QQ).dim == 9      # HZZ Lemma 1.5 (mesh)
    assert len(list(A.quiver.vertices)) == 4


@lit
def test_mesh_gentleness_forces_relation():
    # a*b NOT in I at a special middle vertex is NOT a gentle pair -> not a valid triple.
    from quiverlab.skewgentle.triple import is_skew_gentle_triple
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)})
    assert is_skew_gentle_triple(Q, relations=["a*b"], special={2}) is True
    assert is_skew_gentle_triple(Q, relations=[], special={2}) is False    # VERIFIED


@xeng
def test_mesh_ar_and_tau_tilting_counts():
    # the SUFFICIENCY oracle for the mesh construction (dim alone is not enough).
    from quiverlab.modules.ar import knit_ar_quiver
    from quiverlab.tautilting.mutation import exchange_graph
    A = SkewGentleAlgebra(_t_mesh(), field=QQ)
    ar = knit_ar_quiver(A)
    assert ar.is_complete and len(ar.vertices) == 11            # VERIFIED live
    eg = exchange_graph(A)
    assert eg.is_complete and len(eg.vertices) == 46            # VERIFIED live


@lit
def test_sp_empty_byte_reduces_to_plain_gentle():
    Q = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3)})
    A = SkewGentleAlgebra(quiver=Q, relations=["a*b"], special=(), field=QQ)
    plain = Q.algebra(relations=["a*b"], field=QQ)
    assert A.dim == plain.dim
    assert A.cartan_matrix() == plain.cartan_matrix()
    assert sorted(A.quiver.arrows) == sorted(plain.quiver.arrows)


@lit
def test_chen_selfinjective_iff():
    # Chen Cor 1.2(c) (INDECOMPOSABLE hypothesis): both algebras below are connected.
    Q2 = Quiver([1, 2], {"a": (1, 2), "b": (2, 1)})
    self_inj = SkewGentleAlgebra(quiver=Q2, relations=["a*b", "b*a"], special=(),
                                 field=QQ)
    assert self_inj.is_selfinjective() is True
    not_self = SkewGentleAlgebra(quiver=Q2, relations=["a*b", "b*a"], special={1},
                                 field=QQ)
    assert not_self.is_selfinjective() is False


@lit
def test_kxk_decomposable_selfinjective_with_special():
    # M2 reference: k x k = the Sp={1} triple on a one-vertex quiver. It is semisimple
    # (selfinjective) WITH Sp != empty -- NOT a counterexample to Cor 1.2(c) (which
    # requires INDECOMPOSABLE). VERIFIED live: dim 2, is_selfinjective() True.
    t = SkewGentleTriple.make(Quiver([1], {}), relations=[], special={1})
    A = SkewGentleAlgebra(t, field=QQ)
    assert A.dim == 2 and A.is_selfinjective() is True


@selfcert
def test_broken_triple_refused():
    with pytest.raises(QuiverlabError):
        SkewGentleAlgebra(quiver=Quiver([1], {}), relations=[], special={99}, field=QQ)
