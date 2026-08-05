"""The C4 tau-tilting engine (Plan 45, Adachi-Iyama-Reiten Compos. Math. 150 (2014)).

Support tau-tilting pairs, g-vectors, the mutation exchange graph, the torsion-class
lattice with brick/semibrick labels, 2-term silting, King theta-stability + the
wall-and-chamber fan, and maximal green sequences -- every enumeration budget-capped with
the honest complete-iff-tau-tilting-finite contract, every geometric payload EXACT
rational (no floats in src/). Public surface: ``import quiverlab.tautilting``.
"""
from quiverlab.tautilting.block import tau_tilting_block
from quiverlab.tautilting.green import maximal_green_sequences
from quiverlab.tautilting.mutation import ExchangeGraph, exchange_graph, mutate
from quiverlab.tautilting.pairs import (SupportTauTiltingPair, initial_pair, make_pair,
                                        terminal_pair)
from quiverlab.tautilting.rigid import g_matrix, g_vector, is_tau_rigid
from quiverlab.tautilting.silting import silting_count, two_term_silting
from quiverlab.tautilting.stability import (is_theta_semistable, is_theta_stable,
                                            wall_and_chamber_fan)
from quiverlab.tautilting.torsion import (bricks, hasse_orientation, semibricks,
                                          torsion_class_data)

__all__ = [
    "is_tau_rigid", "g_vector", "g_matrix",
    "SupportTauTiltingPair", "make_pair", "initial_pair", "terminal_pair",
    "mutate", "exchange_graph", "ExchangeGraph",
    "torsion_class_data", "hasse_orientation", "bricks", "semibricks",
    "is_theta_semistable", "is_theta_stable", "wall_and_chamber_fan",
    "maximal_green_sequences",
    "two_term_silting", "silting_count",
    "tau_tilting_block",
]
