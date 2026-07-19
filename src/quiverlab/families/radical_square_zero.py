"""RadicalSquareZero(Q) = kQ/rad^2 (spec §3.4): the truncation at r = 2, whose basis
is the vertices plus the arrows. Always finite-dimensional, even when kQ is not."""
from quiverlab.families.truncated import TruncatedPathAlgebra


def RadicalSquareZero(quiver, field=None):
    A = TruncatedPathAlgebra(quiver, 2, field=field)
    A._family_citations = ("path_algebra", "bardzell")
    return A
