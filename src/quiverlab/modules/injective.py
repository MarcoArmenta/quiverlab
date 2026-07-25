"""Injective resolutions and injective dimension (Plan 23, Tier 1b item 2).

The dual of Plan 05's projective machinery, on the same A^op + D engine:

  E(M)          = D(projective cover of DM over A^op)          (injective envelope)
  E^n           = D(P_n),  P_* = min. proj. resolution of DM over A^op
  inj.dim_A(M)  = pd_{A^op}(DM)

D preserves dimension vectors, so E^n has the dimension vector of P_n and betti(n)
= #injective summands of E^n. `injective_dimension` returns an int, or None when
the module is not resolved within `bound` (infinite -- e.g. self-injective
non-projective)."""
from quiverlab.modules.duality import dualize


class InjectiveResolution:
    """0 -> M -> E^0 -> E^1 -> ... with E^n = D(P_n) injective A-modules."""

    def __init__(self, M, terms, vertices):
        self.module = M
        self.terms = terms              # list of injective A-modules (E^0, E^1, ...); None = 0
        self.vertices = vertices        # summand vertices per term (E^n = (+)_u I_u)
        self.length = len(terms)

    def term(self, n):
        return self.vertices[n] if n < len(self.vertices) else []

    def betti(self, n):
        return len(self.term(n))

    def dimension_vectors(self):
        return [t.dimension_vector() if t is not None else {} for t in self.terms]

    def is_finite(self):
        return any(t is None or t.dim == 0 for t in self.terms)

    def injective_dimension(self):
        for n, t in enumerate(self.terms):
            if t is None or t.dim == 0:
                return n - 1
        return None                     # not resolved within the requested length

    def __repr__(self):
        parts = " -> ".join(
            "(+)".join(f"I_{u}" for u in v) if v else "0" for v in self.vertices)
        return f"inj. res. of {self.module.name}: {parts}"


def injective_resolution(M, length, max_term_dim=200000):
    """Minimal injective coresolution of M to the requested length: dualize the
    minimal projective resolution of DM over A^op term by term."""
    from quiverlab.modules.resolution import minimal_resolution
    DM = dualize(M)
    terms, _ = minimal_resolution(DM, length, max_term_dim=max_term_dim)
    inj_terms, verts = [], []
    for t in terms:
        if t.module is None or t.module.dim == 0:
            inj_terms.append(None)
            verts.append([])
        else:
            inj_terms.append(dualize(t.module))     # E^n = D(P_n), a right A-module
            verts.append(list(t.vertices))
    return InjectiveResolution(M, inj_terms, verts)


def injective_dimension(M, bound=32, max_term_dim=200000):
    """inj.dim_A(M) = pd_{A^op}(DM). Exact int, or None when unresolved within
    `bound` (infinite injective dimension)."""
    DM = dualize(M)
    return DM.projective_resolution(bound, max_term_dim=max_term_dim).pd()
