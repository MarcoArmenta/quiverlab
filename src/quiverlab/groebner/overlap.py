"""Overlap and inclusion ambiguities between leading words, LEFT TO RIGHT.

For leading words w1, w2:
 * OVERLAP: a nonempty b that is a proper SUFFIX of w1 and a proper PREFIX of w2.
   Then w1 = a*b, w2 = b*c, ambiguity word W = a*b*c = w1*c = a*w2. (Right-to-left
   treatments swap prefix/suffix; we are left-to-right per the quiver convention.)
 * INCLUSION: w_inner is a proper contiguous factor of w_outer, w_outer = a*w_inner*c
   with (a, c) != ((), ()); ambiguity word W = w_outer.
Only ambiguity words of length <= degree_bound are produced (spec §5, component 3).
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class Ambiguity:
    kind: str            # "overlap" | "inclusion"
    word: tuple          # the ambiguity word containing both leading words
    left: object         # overlap: rule with lead a*b ; inclusion: the OUTER rule
    right: object        # overlap: rule with lead b*c ; inclusion: the INNER rule
    a: tuple             # left/prefix context
    c: tuple             # right/suffix context


def overlaps(r1, r2, degree_bound):
    w1, w2 = r1.lead, r2.lead
    out = []
    # b nonempty, strictly shorter than both w1 and w2 (proper suffix / proper prefix)
    for blen in range(1, min(len(w1), len(w2))):
        b = w1[len(w1) - blen:]
        if b != w2[:blen]:
            continue
        a = w1[:len(w1) - blen]
        c = w2[blen:]
        word = a + b + c                      # == w1 + c
        if len(word) <= degree_bound:
            out.append(Ambiguity(kind="overlap", word=word, left=r1, right=r2, a=a, c=c))
    return out


def inclusions(r_outer, r_inner, degree_bound):
    W, w = r_outer.lead, r_inner.lead
    if len(w) >= len(W):                      # inner must be a PROPER factor
        return []
    if len(W) > degree_bound:
        return []
    out = []
    for i in range(len(W) - len(w) + 1):
        if W[i:i + len(w)] == w:
            a, c = W[:i], W[i + len(w):]
            if (a, c) == ((), ()):
                continue
            out.append(Ambiguity(kind="inclusion", word=W, left=r_outer, right=r_inner, a=a, c=c))
    return out


def all_ambiguities(rules, degree_bound):
    out = []
    for r1 in rules:
        for r2 in rules:
            out.extend(overlaps(r1, r2, degree_bound))
            if r1 is not r2:
                out.extend(inclusions(r1, r2, degree_bound))
    return out
