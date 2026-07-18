"""Linear-combination arithmetic, reduction rules, and normal forms over a Domain.

An ELEMENT of kQ is a dict  word -> domain element  with no zero coefficients;
the zero element is {}. A word is a tuple of arrow names read LEFT TO RIGHT.
Reduction rewrites an occurrence of a rule's leading word W = u*lead*v into
u*tail*v; every tail word is strictly smaller than lead in the admissible order,
so reduction terminates (spec §5, component 3)."""
from dataclasses import dataclass

from quiverlab.groebner.events import ReductionStep


def lc_add(comb, word, coeff, dom):
    """comb[word] += coeff, in place; drop the entry if it becomes zero."""
    if dom.is_zero(coeff):
        return
    cur = comb.get(word)
    total = coeff if cur is None else dom.add(cur, coeff)
    if dom.is_zero(total):
        comb.pop(word, None)
    else:
        comb[word] = total


def lc_sub(a, b, dom):
    """Return a - b as a fresh element."""
    out = dict(a)
    for w, c in b.items():
        lc_add(out, w, dom.neg(c), dom)
    return out


@dataclass(frozen=True)
class ReductionRule:
    """lead -> tail: tail is a linear combination of words strictly SMALLER than
    lead, parallel to it (same source & target), with lead - sum(tail) in I.
    Stored as tail = ((coeff, word), ...) sorted descending by the order."""
    lead: tuple
    tail: tuple
    source: object
    target: object


def rule_from_comb(comb, order, quiver, dom):
    """Turn a nonzero element into a rule lead -> tail by dividing out the leading
    coefficient: lead = LM(comb), tail = -(1/LC)*(comb - LC*lead)."""
    lead = order.leading(comb)
    inv = dom.inv(comb[lead])
    tail_items = [(dom.neg(dom.mul(inv, c)), w) for w, c in comb.items() if w != lead]
    tail_items.sort(key=lambda cw: order.key(cw[1]), reverse=True)
    return ReductionRule(
        lead=lead,
        tail=tuple(tail_items),
        source=quiver.word_source(lead),
        target=quiver.word_target(lead),
    )


def first_factor(word, rules):
    """Leftmost occurrence of any rule's lead as a contiguous factor of word.
    Returns (rule, position) with the smallest position (ties: first such rule),
    or None if word is irreducible."""
    best = None
    for rule in rules:
        L = rule.lead
        n = len(L)
        for i in range(len(word) - n + 1):
            if word[i:i + n] == L:
                if best is None or i < best[1]:
                    best = (rule, i)
                break
    return best


def reduce_comb(comb, rules, order, dom, trace=None):
    """Normal form of comb under rules. Deterministic: at each step rewrite the
    LARGEST (in the order) reducible word. Terminates because every rewrite
    replaces a word by strictly smaller ones in a well-founded order; for a
    confluent system the result is independent of these choices."""
    work = {w: c for w, c in comb.items() if not dom.is_zero(c)}
    while True:
        target_word = None
        for w in sorted(work, key=order.key, reverse=True):
            if first_factor(w, rules) is not None:
                target_word = w
                break
        if target_word is None:
            return work
        rule, i = first_factor(target_word, rules)
        coeff = work.pop(target_word)
        u, v = target_word[:i], target_word[i + len(rule.lead):]
        if trace is not None:
            before = dict(work)
            before[target_word] = coeff
        for tc, tw in rule.tail:
            lc_add(work, u + tw + v, dom.mul(coeff, tc), dom)
        if trace is not None:
            trace.append(ReductionStep(word=target_word, rule_lead=rule.lead,
                                       before=before, after=dict(work)))
