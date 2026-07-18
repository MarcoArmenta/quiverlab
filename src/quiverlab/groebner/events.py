"""Plain-dataclass step events emitted by the Groebner lowering, as trace hooks.

These are INERT: when a caller passes trace=[...] the lowering appends these
records; nothing else consumes them yet. The formal trace subsystem (typed
events, PDF/HTML/text renderers, eliding rules, golden-file tests) is Plan 07 --
this is only the emission boundary.
"""
from dataclasses import dataclass


@dataclass
class Dispatch:
    """Which lowering route Quiver.algebra took, and why."""
    route: str          # "monomial" | "groebner"
    reason: str
    n_relations: int


@dataclass
class ReductionStep:
    """One rewrite: the word occurrence reduced, the rule's leading word, and the
    linear combination (word -> domain element) before and after the step."""
    word: tuple
    rule_lead: tuple
    before: dict
    after: dict
