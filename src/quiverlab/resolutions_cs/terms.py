"""Ambiguity-chain type for the CS S-sequence. Pure combinatorics; words are tuples of
arrow NAMES read left-to-right (matching the Plan-03 reduction system)."""
from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class Chain:
    """sigma in S_n = 𝒜_{n-1}. `word`: the underlying path (arrow names). `blocks`: the
    unique CS left decomposition u_0|...|u_{n-1} (n blocks). `o`,`t`: source/target. `degree`: n."""
    word: tuple
    blocks: tuple
    o: object
    t: object
    degree: int

    @property
    def n_blocks(self) -> int:
        return len(self.blocks)
