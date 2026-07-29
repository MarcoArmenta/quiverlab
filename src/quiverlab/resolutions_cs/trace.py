"""Plain trace dataclasses. Plan 07 renders these; Plan 04 only populates them and
asserts their claims equal computed values."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any


@dataclass
class AmbiguityEvent:
    degree: int
    chain_words: list


@dataclass
class ResolutionTerm:
    degree: int
    n_generators: int          # |S_n|
    collapsed_dim: int         # dim C_n / dim C^n
    corners: Any = None
    """The generators' (source, target) vertex pairs, one per generator WITH
    repetition -- so the report can NAME the term as a direct sum of projective
    bimodules ``C_n = (+)_{s in S_n} A e_{o(s)} (x) e_{t(s)} A`` instead of quoting
    a bare generator count (Marco 2026-07-29).

    ``None`` when the engine's term is not vertex-decomposed (the bar resolution
    over a structure-constants algebra) -- the renderers then omit the line rather
    than invent a decomposition. Defaulted, so every existing construction and
    every recorded stream is unchanged."""


@dataclass
class DifferentialEvent:
    degree: int
    chain: Any                 # source chain word
    terms: list                # [(coeff, a_word, target_word, c_word), ...]


@dataclass
class LiftStep:
    degree: int
    kind: str                  # "delta" | "correction-solve" | "dd-check" | "order-check"
    detail: Any = None
