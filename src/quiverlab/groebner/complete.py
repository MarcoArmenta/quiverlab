"""S-polynomials and Buchberger-Mora completion, LEFT TO RIGHT (spec §5, c.3).

For an ambiguity word W the S-polynomial reduces W two ways and subtracts:
 * OVERLAP (w1=a*b, w2=b*c, W=a*b*c): apply the w1 rule at the front (tail(w1)*c)
   minus apply the w2 rule at the back (a*tail(w2)).
 * INCLUSION (W=w_outer=a*w_inner*c): apply the outer rule (tail(w_outer)) minus
   apply the inner rule inside (a*tail(w_inner)*c).
Completion reduces each S-polynomial by the current rules; a nonzero normal form
yields a new rule. Only ambiguity words of length <= degree_bound are formed, and
there are finitely many words of length <= degree_bound over a finite arrow set,
so completion ALWAYS TERMINATES under a fixed bound."""
from quiverlab.errors import AdmissibilityError
from quiverlab.groebner.overlap import all_ambiguities
from quiverlab.groebner.reduction import lc_add, lc_sub, reduce_comb, rule_from_comb


def _apply_tail(rule, prefix, suffix, dom):
    """The element prefix * tail(rule) * suffix (each tail word wrapped in context)."""
    out = {}
    for tc, tw in rule.tail:
        lc_add(out, prefix + tw + suffix, tc, dom)
    return out


def s_polynomial(amb, dom):
    if amb.kind == "overlap":
        front = _apply_tail(amb.left, (), amb.c, dom)     # tail(w1) * c
        back = _apply_tail(amb.right, amb.a, (), dom)      # a * tail(w2)
        return lc_sub(front, back, dom)
    outer = _apply_tail(amb.left, (), (), dom)             # tail(w_outer)
    inner = _apply_tail(amb.right, amb.a, amb.c, dom)      # a * tail(w_inner) * c
    return lc_sub(outer, inner, dom)


def _is_factor(short, long):
    """True if `short` occurs as a contiguous factor of `long`."""
    n = len(short)
    return any(long[i:i + n] == short for i in range(len(long) - n + 1))


def _minimize_leads(rules):
    """Drop rules whose leading word contains ANOTHER rule's leading word as a
    proper factor (such a lead is already reducible, so the rule is redundant), and
    drop exact-duplicate leads (keep the first). The remaining leads form an
    antichain under the factor relation, each appearing once. This does not change
    the ideal, the set of irreducible words (forbidden = leads; a subsumed forbidden
    word is redundant), or confluence -- it is the standard reduced-basis cleanup,
    it keeps max-leading-length honest for the degree-bound check, and it keeps
    leading_words()/ambiguities() (the CS S-sequence seed) free of duplicates."""
    minimal, seen = [], set()
    for r in rules:
        if r.lead in seen:                          # exact-duplicate lead: keep the first
            continue
        if any(s.lead != r.lead and len(s.lead) < len(r.lead) and _is_factor(s.lead, r.lead)
               for s in rules):
            continue
        minimal.append(r)
        seen.add(r.lead)
    return minimal


def complete(init_rules, order, quiver, dom, degree_bound, trace=None, max_rules=20000):
    rules = list(init_rules)
    added = True
    while added:
        added = False
        for amb in all_ambiguities(rules, degree_bound):
            g = reduce_comb(s_polynomial(amb, dom), rules, order, dom, trace=trace)
            if g:                                          # nonzero normal form
                rules.append(rule_from_comb(g, order, quiver, dom))
                added = True
                if len(rules) > max_rules:
                    raise AdmissibilityError(
                        f"Groebner completion produced more than {max_rules} rules under "
                        f"degree_bound={degree_bound}: the reduction system is not stabilizing",
                        hint="the ideal may not be admissible; inspect the relations, or this "
                             "presentation is beyond the v1 certificate",
                    )
                break                                      # restart with the enlarged system
    return _minimize_leads(rules)
