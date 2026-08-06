"""The gentle / string subsystem (Plan 46 / C5).

Butler-Ringel string & band module classification, string-module tau by
hooks/cohooks, the Avella-Alaminos-Geiss derived invariant for gentle algebras,
and the algebra-only ``strings`` compute block. Everything is exact combinatorics
+ exact linear algebra over a Domain; float-free by design.

Public surface (reach via ``import quiverlab`` or ``quiverlab.strings``):
  - walks:      Letter alphabet, ``string_signs``, ``is_valid_walk``,
                ``enumerate_strings``, ``find_bands``, ``StringCensus``.
  - modules:    ``string_module``, ``band_module``.
  - ar_strings: ``string_tau``, ``string_tau_minus``.
  - ag:         ``permitted_threads``, ``forbidden_threads``, ``ag_invariant``,
                ``AGInvariant``.
  - block:      ``strings_block``.
"""
from quiverlab.strings.ag import (AGInvariant, ag_invariant, forbidden_threads,
                                  permitted_threads)
from quiverlab.strings.ar_strings import string_tau, string_tau_minus
from quiverlab.strings.block import strings_block
from quiverlab.strings.modules import band_module, string_module
from quiverlab.strings.walks import (StringCensus, enumerate_strings, find_bands,
                                     invert, is_valid_walk, letter_source,
                                     letter_target, string_signs)

__all__ = [
    "StringCensus", "enumerate_strings", "find_bands", "invert", "is_valid_walk",
    "letter_source", "letter_target", "string_signs",
    "band_module", "string_module",
    "string_tau", "string_tau_minus",
    "AGInvariant", "ag_invariant", "forbidden_threads", "permitted_threads",
    "strings_block",
]
