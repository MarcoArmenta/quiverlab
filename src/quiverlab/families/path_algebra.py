"""PathAlgebra: the hereditary kQ of a Dynkin/Euclidean quiver (spec §3.4).
relations=[]; finite-dimensional exactly when the chosen orientation is acyclic --
a cyclic orientation is refused loudly by the monomial route."""
from quiverlab.combinat.quiver import Quiver
from quiverlab.families.dynkin import dynkin_quiver


def PathAlgebra(type_or_quiver, orientation="linear", field=None):
    Q = type_or_quiver if isinstance(type_or_quiver, Quiver) else \
        dynkin_quiver(type_or_quiver, orientation)
    A = Q.algebra(relations=[], field=field)     # hereditary; loud if cyclic orientation
    A._family_citations = ("path_algebra", "happel_question", "assem_book")
    return A
