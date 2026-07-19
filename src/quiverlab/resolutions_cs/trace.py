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
