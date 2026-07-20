"""print(A) shows vertices, arrows, relations (spec §3.7) without changing the
HHTable title contract (repr line 0 is used verbatim as the table title)."""
from quiverlab import Quiver, CC


def test_repr_line0_is_unchanged():
    A = Quiver([1, 2], {"a": (1, 2)}).algebra(field=CC)
    assert repr(A).splitlines()[0] == f"Algebra of dimension {A.dim} over {A.domain.name}"


def test_repr_appends_quiver_summary():
    Q = Quiver([1, 2, 3, 4], {"a": (1, 2), "b": (2, 4), "c": (1, 3), "d": (3, 4)})
    A = Q.algebra(relations=["a*b - c*d"], field=CC)
    text = repr(A)
    assert "vertices: 1, 2, 3, 4" in text
    assert "a: 1 -> 2" in text and "d: 3 -> 4" in text
    assert "relations:" in text and "a*b - c*d" in text


def test_repr_without_quiver_has_no_summary():
    # An Algebra built without a quiver (e.g. the structure-constant escape hatch)
    # has .quiver = None, so the vertices/arrows/relations block is NOT appended.
    from quiverlab.core.algebra import Algebra
    from quiverlab.fields import QQ
    one = QQ.one()
    A = Algebra(QQ, [[[one]]], [one])   # 1-dim algebra k; _quiver defaults to None
    assert A.quiver is None
    assert "vertices:" not in repr(A)
