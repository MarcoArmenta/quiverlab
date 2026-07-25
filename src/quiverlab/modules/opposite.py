"""The opposite algebra A^op as a first-class Algebra (Plan 23, Tier 1b).

A^op is the SAME k-space with the reversed product ``x .^op y := y . x``. We
realize it with a reversed quiver, transposed structure constants, reversed basis
labels, and reversed relations, so that all the module machinery (projectives,
resolutions, Hom, rad/top/soc) runs over A^op unchanged.

Path composition is left-to-right (Assem-Simson-Skowronski): reversing an arrow
``a: s -> t`` to ``a: t -> s`` flips source/target, AND the product order flips
(``T^op[i][j] = T[j][i]``). The two flips are consistent: for paths p, q,
``reverse(q.p) = reverse(p) * reverse(q)`` in Q^op, matching ``p .^op q = q . p``.

The index set is PRESERVED (index i in A <-> index i in A^op, same vector,
reversed label). This makes the duality D and the transpose Tr coordinate-clean
(``modules/duality.py``) and avoids re-running Groebner completion for A^op.
"""
from quiverlab.combinat.quiver import Quiver
from quiverlab.combinat.relations import Relation
from quiverlab.core.algebra import Algebra
from quiverlab.errors import QuiverlabError


def reverse_label(label):
    """Reverse a quiverlab basis label: 'e_v' is fixed, 'a*b*c' -> 'c*b*a'."""
    if label.startswith("e_"):
        return label
    return "*".join(reversed(label.split("*")))


def _reverse_relations(relations, quiver_op):
    """Reverse each relation word; the parallel (source, target) pair swaps too.
    Only ``.terms`` and ``repr`` are consumed downstream (Module.check_module)."""
    out = []
    for rel in (relations or []):
        terms = tuple((c, tuple(reversed(w))) for c, w in rel.terms)
        out.append(Relation(terms, source=rel.target, target=rel.source))
    return out


def opposite_algebra(A):
    """Return A^op as an Algebra (reversed quiver, transposed structure constants).

    Requires the quiver/path provenance (structure-constant-only algebras carry no
    path basis to reverse). Cached and cross-linked so ``opposite_algebra`` is an
    exact involution: ``opposite_algebra(opposite_algebra(A)) is A``."""
    cached = getattr(A, "_opposite_algebra", None)
    if cached is not None:
        return cached
    if A.quiver is None or A.basis_labels is None:
        raise QuiverlabError(
            "opposite algebra needs the quiver presentation",
            hint="construct the algebra via Quiver.algebra(...); structure-constant "
                 "algebras carry no path basis to reverse")
    Q = A.quiver
    Qop = Quiver(list(Q.vertices),
                 {name: (t, s) for name, (s, t) in Q.arrows.items()})
    m = A.dim
    # T^op[i][j] = coords of b_i .^op b_j = b_j . b_i = T[j][i]
    Top = [[list(A.T[j][i]) for j in range(m)] for i in range(m)]
    labels = [reverse_label(lab) for lab in A.basis_labels]
    rels_op = _reverse_relations(A.relations, Qop)
    Aop = Algebra(A.domain, Top, list(A.unit), basis_labels=labels,
                  is_unit_adapted=A.is_unit_adapted, _quiver=Qop, _relations=rels_op,
                  _family_citations=getattr(A, "_family_citations", ()))
    A._opposite_algebra = Aop
    Aop._opposite_algebra = A
    return Aop
