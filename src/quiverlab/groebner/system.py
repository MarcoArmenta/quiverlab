"""ReductionSystem: the certified confluent reduction system for kQ/I, and the
object Chouhy-Solotar (Plan 04, arXiv:1406.2300) consumes verbatim (spec §5 c.3,
§6). Pairs (leading word, reduction) + the admissible order + the ambiguity
S-sequence, satisfying CS's reduction-finiteness requirement (finitely many rules,
finitely many irreducible words, every ambiguity resolves). LEFT TO RIGHT."""
from dataclasses import dataclass

from quiverlab.errors import RelationError
from quiverlab.groebner.order import PathOrder, path_order
from quiverlab.groebner.reduction import ReductionRule, reduce_comb, rule_from_comb
from quiverlab.groebner.overlap import all_ambiguities
from quiverlab.groebner.complete import complete
from quiverlab.groebner.certificate import (
    default_degree_bound, check_degree_bound, certified_irreducibles,
)


@dataclass(frozen=True)
class ReductionSystem:
    quiver: object
    domain: object
    order: PathOrder
    rules: tuple
    irreducibles: tuple
    degree_bound: int
    is_confluent: bool = True

    def leading_words(self):
        return tuple(r.lead for r in self.rules)

    def reduce(self, comb):
        """Normal form of a linear combination (word -> domain element)."""
        return reduce_comb(comb, self.rules, self.order, self.domain)

    def normal_form(self, word):
        """Reduce a single word to its normal form (a linear combination)."""
        return reduce_comb({word: self.domain.one()}, self.rules, self.order, self.domain)

    def ambiguities(self):
        """All overlap/inclusion ambiguities of the completed system (each resolves
        to 0); the seed of the Chouhy-Solotar S-sequence (spec §6)."""
        return tuple(all_ambiguities(self.rules, self.degree_bound))


def _relations(quiver, relations):
    if relations and isinstance(relations[0], str):
        from quiverlab.combinat.relations import parse_relations
        return parse_relations(list(relations), quiver)
    return list(relations)


def build_reduction_system(quiver, relations, field, degree_bound=None, trace=None):
    from quiverlab.errors import AdmissibilityError
    rels = _relations(quiver, relations)
    for rel in rels:
        if rel.min_length < 2:
            raise AdmissibilityError(
                f"relation {rel!r} has a path of length {rel.min_length}: the ideal is not "
                "inside the square of the arrow ideal",
                hint="admissible relations use paths of length >= 2",
            )
    # Build the domain from the coefficient entries (0 and 1 always included so the
    # domain matches the monomial route's on monomial inputs).
    raw = [field.parse_entry(0), field.parse_entry(1)]
    for rel in rels:
        for c, _w in rel.terms:
            raw.append(field.parse_entry(c))
    dom = field.make_domain(raw)
    order = path_order(quiver)
    init_rules = []
    for rel in rels:
        comb = {}
        for c, w in rel.terms:
            cc = dom.coerce(field.parse_entry(c))
            if not dom.is_zero(cc):
                comb[w] = cc
        if not comb:
            raise RelationError(
                f"relation {rel!r} vanishes in {dom.name}, so it says 0 = 0",
                hint="a relation must have a nonzero coefficient in the chosen field; "
                     "drop it or change the field",
            )
        init_rules.append(rule_from_comb(comb, order, quiver, dom))
    if degree_bound is None:
        degree_bound = default_degree_bound(rels)
    rules = complete(init_rules, order, quiver, dom, degree_bound, trace=trace)
    check_degree_bound(rules, degree_bound)                 # AdmissibilityError if too small
    irreducibles = certified_irreducibles(quiver, rules)    # NotFiniteDimensionalError if infinite
    return ReductionSystem(quiver=quiver, domain=dom, order=order, rules=tuple(rules),
                           irreducibles=irreducibles, degree_bound=degree_bound)
