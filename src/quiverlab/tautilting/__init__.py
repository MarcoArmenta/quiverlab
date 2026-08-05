"""The C4 tau-tilting engine (Plan 45, Adachi-Iyama-Reiten Compos. Math. 150 (2014)).

Support tau-tilting pairs, g-vectors, the mutation exchange graph, the torsion-class
lattice with brick/semibrick labels, 2-term silting, King theta-stability + the
wall-and-chamber fan, and maximal green sequences -- every enumeration budget-capped with
the honest complete-iff-tau-tilting-finite contract, every geometric payload EXACT
rational (no floats in src/). Public surface: ``import quiverlab.tautilting``.
"""
from quiverlab.tautilting.rigid import g_matrix, g_vector, is_tau_rigid

__all__ = ["is_tau_rigid", "g_vector", "g_matrix"]
