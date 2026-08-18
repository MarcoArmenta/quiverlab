"""The Batalin-Vilkovisky operator Delta on HH of Frobenius / self-injective
algebras (Plan 54, R2).

Delta: HH^n -> HH^{n-1} is Connes' operator B carried across the sigma-twisted
Frobenius duality HH^n(A) ~= D(HH_n(A, {}_1A_nu)). The public entry point is
:meth:`quiverlab.core.Algebra.bv_operator`; the pieces live here:

- :mod:`hypothesis` -- the decidable nu hypothesis gate (symmetric anchor /
  semisimple-nu LZZ route / loud refusals);
- :mod:`transport` -- the perfect pairing (dagger) and the transport (double-dagger)
  giving Delta (symmetric and twisted routes) + the ``BVOperator`` assembly;
- :mod:`bracket` -- the bracket recovered from Delta via the BV relation and the
  in-window arbiter against the independent Gerstenhaber bracket;
- :mod:`twist` -- the P52 twisted-homology adapter (cross-plan contract);
- :mod:`twisted_connes` -- the twisted Connes operator B_sigma on the twisted
  bar complex.
"""
