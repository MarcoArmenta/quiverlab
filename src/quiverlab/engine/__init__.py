"""quiverlab.engine: the F_p compute engine, ported from hanlab.

Ported from hanlab (HansConjecture, MIT (c) 2026 Marco Armenta,
github.com/marcoarmenta/hansconjecture), bank state of 2026-07-12.

INTERNAL API. The engine computes over prime fields F_p with exact int64
arithmetic (structure constants held unreduced, reduced mod p at rank time;
p = 32003 is the char-0 proxy). Public quiverlab entry points dispatch here
automatically for algebras over GF(p); everything else uses the pure Plan-01
paths. Set QUIVERLAB_NO_NUMBA=1 to force the pure-Python kernels.

Modules keep their hanlab development names (hh_engine, scan3, coxeter, ...).
"""
