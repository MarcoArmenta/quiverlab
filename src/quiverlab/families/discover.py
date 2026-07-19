"""The v1 family catalog + families() discoverability (spec §3.4)."""
from dataclasses import dataclass


@dataclass(frozen=True)
class FamilyInfo:
    name: str
    signature: str
    route: str          # "monomial" | "general" | "structure-constant" | "iterator"
    citations: tuple
    summary: str


CATALOG = (
    FamilyInfo("NakayamaAlgebra", "NakayamaAlgebra([c1,..,cn]) | (n=, l=, cyclic=)",
               "monomial", ("nakayama", "assem_book"),
               "Serial algebra by Kupisch series (linear or cyclic)."),
    FamilyInfo("PathAlgebra", "PathAlgebra(type, orientation='linear')",
               "monomial", ("path_algebra", "happel_question", "assem_book"),
               "Hereditary path algebra of a Dynkin/Euclidean quiver."),
    FamilyInfo("TruncatedPathAlgebra", "TruncatedPathAlgebra(type_or_Q, r)",
               "monomial", ("path_algebra", "bardzell"),
               "kQ / rad^r."),
    FamilyInfo("RadicalSquareZero", "RadicalSquareZero(Q)",
               "monomial", ("path_algebra", "bardzell"),
               "kQ / rad^2 from any quiver."),
    FamilyInfo("IncidenceAlgebra", "IncidenceAlgebra(covers)",
               "general", ("incidence", "assem_book", "hodge"),
               "Incidence algebra of a finite poset (Hasse quiver + commutativity)."),
    FamilyInfo("QuantumCI", "QuantumCI(q=...)",
               "general", ("quantum_ci", "qci_hh_oracle", "bardzell", "chouhy_solotar"),
               "Quantum complete intersection k<x,y>/(x^2,y^2,xy+q yx)."),
    FamilyInfo("ExteriorAlgebra", "ExteriorAlgebra(n)",
               "general", ("quantum_ci", "chouhy_solotar"),
               "Exterior algebra Lambda(k^n), dim 2^n."),
    FamilyInfo("PreprojectiveAlgebra", "PreprojectiveAlgebra(type)",
               "general", ("preprojective", "chouhy_solotar", "assem_book"),
               "Preprojective algebra of a Dynkin quiver."),
    FamilyInfo("TrivialExtension", "TrivialExtension(A)",
               "structure-constant", ("assem_book",),
               "Trivial extension A |x D(A) (symmetric)."),
    FamilyInfo("TensorProduct", "TensorProduct(A, B)",
               "structure-constant", ("tensor_product", "hodge"),
               "Tensor product A (x)_k B."),
    FamilyInfo("zoo", "zoo(dim_max=12)",
               "iterator", ("han_conjecture", "chouhy_solotar"),
               "Iterator over the curated exact zoo of open (Han-conjecture) algebras."),
)


class FamilyListing:
    def __init__(self, catalog=CATALOG):
        self._catalog = tuple(catalog)

    def __iter__(self):
        return iter(self._catalog)

    def names(self):
        return tuple(f.name for f in self._catalog)

    def by_name(self, name):
        for f in self._catalog:
            if f.name == name:
                return f
        raise KeyError(name)

    def to_dict(self):
        return {"families": [f.__dict__ | {"citations": list(f.citations)}
                             for f in self._catalog]}

    def __str__(self):
        w = max(len(f.name) for f in self._catalog)
        rows = [f"{f.name:<{w}}  [{f.route}]  {f.signature}" for f in self._catalog]
        return "quiverlab families (spec 3.4):\n" + "\n".join("  " + r for r in rows)


def families():
    return FamilyListing()
