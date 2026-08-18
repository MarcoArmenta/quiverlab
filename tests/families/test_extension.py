"""Arrow-removal subalgebra ``B <= A`` -- the "extension by arrows and relations"
input model (CLMS ``2101.02597`` Def. 5.1/5.2). Plan 73 / R7, Task 0.

The two literature examples are CLMS Ex. 5.3 (a BOUNDED extension, ``dim`` 47/20/27)
and Ex. 5.5 (NOT bounded, ``dim`` 10/6/4). Composition is left-to-right, so every
CLMS right-to-left path is translated (see the plan's Mathematical foundation).
"""
import pytest

from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import QuiverlabError
from quiverlab.fields import QQ
from quiverlab.families.extension import Extension, arrow_removal_subalgebra
from quiverlab.families.tensor import TensorProduct

lit = pytest.mark.oracle_literature
selfcert = pytest.mark.oracle_selfcert


def build_ex53(field=QQ):
    """CLMS Ex. 5.3 (bounded). ``F = {a}``; ``dim`` 47/20/27; ``gl.dim B = 2``."""
    return Quiver(
        [1, 2, 3, 4, 5],
        {"al": (5, 1), "be": (1, 5), "d": (1, 4), "c": (4, 3), "b": (3, 2), "a": (2, 1)},
    ).algebra(relations=["al*be", "d*c*b*a - be*al"], field=field)


def build_ex55(field=QQ):
    """CLMS Ex. 5.5 (NOT bounded). ``F = {d}``; ``dim`` 10/6/4."""
    return Quiver(
        [1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (1, 3), "d": (3, 1)}
    ).algebra(relations=["a*b", "d*a", "b*d", "d*c*d"], field=field)


@lit
def test_ex53_dims():
    ext = arrow_removal_subalgebra(build_ex53(), ("a",))
    assert isinstance(ext, Extension)
    assert (ext.dim_A, ext.dim_B, ext.dim_quotient) == (47, 20, 27)
    assert tuple(ext.new_arrows) == ("a",)
    assert ext.B.dim == 20


@lit
def test_ex55_dims_and_relative_paths():
    ext = arrow_removal_subalgebra(build_ex55(), ("d",))
    assert (ext.dim_A, ext.dim_B, ext.dim_quotient) == (10, 6, 4)
    # A/B basis = the A-basis paths using >= 1 new arrow (CLMS Def. 5.9-5.10)
    assert set(ext.relative_paths()) == {"d", "c*d", "d*c", "c*d*c"}
    assert len(ext.relative_paths()) == 4


@selfcert
def test_dim_B_certified_equals_ffree_count():
    """The extension axiom ``J cap B = 0`` == ``dim B`` matches the F-free basis count
    (present_from_pi certifies dim(kQ_B/ker pi))."""
    ext = arrow_removal_subalgebra(build_ex53(), ("a",))
    assert ext.dim_B == 20
    assert ext.dim_A - ext.dim_B == ext.dim_quotient == 27


@selfcert
def test_al_prefix_not_matched_as_a():
    """The 'al' vs 'a' whole-token trap: removing 'a' must NOT treat 'al' as new."""
    ext = arrow_removal_subalgebra(build_ex53(), ("a",))
    rel = ext.relative_paths()
    assert "a" in rel                     # the new arrow itself is a relative path
    assert "be*al" not in rel             # F-free (no whole 'a' token) -> in B
    assert "a*be*al" in rel               # uses the whole 'a' token -> relative


def test_structure_constant_algebra_refused():
    """A quiver-less (structure-constant) algebra is refused loudly."""
    sc = TensorProduct(build_ex55(), build_ex55())   # TensorProduct -> quiver is None
    assert sc.quiver is None
    with pytest.raises(QuiverlabError):
        arrow_removal_subalgebra(sc, ("d",))


def test_unknown_arrow_refused():
    with pytest.raises(QuiverlabError):
        arrow_removal_subalgebra(build_ex55(), ("nope",))


# --------------------------------------------------------------------------- #
# the adjudicated P72/P73 seam gate (DD-B2): P72's inert-arrow remove_arrows is
# the special case of P73's general arrow_removal_subalgebra on an INERT F.
# --------------------------------------------------------------------------- #
@selfcert
def test_seam_remove_arrows_equals_subalgebra_on_inert_F():
    from quiverlab.hochschild.arrow_removal import inert_arrows, remove_arrows

    # c is parallel to a and appears in NO relation -> inert (CLMS Def. 3.1)
    A = Quiver([1, 2, 3], {"a": (1, 2), "b": (2, 3), "c": (1, 2)}
               ).algebra(relations=["a*b"], field=QQ)
    assert inert_arrows(A) == ["c"]
    B_p72 = remove_arrows(A, ("c",))                       # P72 inert reduction
    B_p73 = arrow_removal_subalgebra(A, ("c",)).B          # P73 general constructor
    # on an inert F the two constructions agree exactly (I_B = I_A verbatim)
    assert B_p72.dim == B_p73.dim
    assert B_p72.basis_labels == B_p73.basis_labels
    assert [repr(r) for r in B_p72.relations] == [repr(r) for r in B_p73.relations]
    assert (B_p72.hochschild_homology(3, verbose=False).dims
            == B_p73.hochschild_homology(3, verbose=False).dims)
