"""Path-type basis extraction: the split A = E (+) r read off basis_labels.

A *path-type basis* (every Quiver(...).algebra(...) instance) has labels
'e_v' for a complete set of orthogonal idempotents spanning a Wedderburn
complement E, and radical labels spanning r = rad A, each with a unique
source and target idempotent. This is the structural hypothesis behind the
generic any-Domain engine-backed invariants (Plan 19); it is VERIFIED here
— structurally, via multiplication, never by parsing path labels — and the
refusal for non-path-type bases is the honest one (no 'later phase')."""
from quiverlab.errors import FieldError

_HINT = ("present the algebra via Quiver(...).algebra(...) (path-type basis), "
         "or construct it over a prime field GF(p), where the fast engine "
         "serves any basis")


def path_type_basis(A, what="this invariant"):
    """(idem, rad, src, tgt): idempotent and radical basis-index lists plus,
    per radical index i, the unique v, w in idem with e_v f_i = f_i = f_i e_w.
    Raises FieldError if the basis is not path-type."""
    dom = A.domain
    labels = A.basis_labels
    if A.quiver is None or not labels:
        raise FieldError(
            f"{what} over {dom.name} needs a quiver presentation "
            f"(path-type basis); this algebra has none", hint=_HINT)

    def is_basis_vec(vec, j):
        return all(dom.eq(c, dom.one() if t == j else dom.zero())
                   for t, c in enumerate(vec))

    idem = [i for i, lab in enumerate(labels) if lab.startswith("e_")]
    rad = [i for i, lab in enumerate(labels) if not lab.startswith("e_")]
    bad = FieldError(
        f"{what} over {dom.name}: basis is not path-type "
        f"(idempotent/vertex structure not visible in the labels)", hint=_HINT)
    # complete orthogonal idempotents summing to 1
    for v in idem:
        for w in idem:
            prod = A.T[v][w]
            if v == w:
                if not is_basis_vec(prod, v):
                    raise bad
            elif any(not dom.is_zero(c) for c in prod):
                raise bad
    idemset = set(idem)
    for t, c in enumerate(A.unit):
        want = dom.one() if t in idemset else dom.zero()
        if not dom.eq(c, want):
            raise bad
    src, tgt = {}, {}
    for i in rad:
        sv = [v for v in idem if is_basis_vec(A.T[v][i], i)]
        tw = [w for w in idem if is_basis_vec(A.T[i][w], i)]
        if len(sv) != 1 or len(tw) != 1:
            raise bad
        src[i], tgt[i] = sv[0], tw[0]
    return idem, rad, src, tgt
