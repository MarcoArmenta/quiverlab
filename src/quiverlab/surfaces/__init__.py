"""Marked surfaces -> ideal triangulations -> gentle Jacobian algebras (Plan 48).

A no-code INPUT method: pick (or build) a marked bordered surface, triangulate it, and
read off the associated gentle algebra Jac(Q(T), W(T)) (Fomin-Shapiro-Thurston 2008;
Labardini-Fragoso 2009; Assem-Bruestle-Charbonneau-Jodoin-Plamondon 2010). The v1 scope
is the fully-certifiable case -- UNPUNCTURED surfaces with non-empty boundary -- where the
Jacobian is gentle (self-certified via P46's is_gentle + P44's finiteness certificate).
Punctures, closed surfaces, and self-folded triangles refuse loudly (successor P48.1).
Float-free (all surface data is exact integers / tuples)."""
from quiverlab.surfaces.marked import MarkedSurface

__all__ = ["MarkedSurface"]
