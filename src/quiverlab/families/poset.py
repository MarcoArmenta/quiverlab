"""Finite posets from cover data (spec §3.4). covers = [(x, y), ...] means x is
covered by y (Hasse edge, arrow x->y). Order <= is the reflexive-transitive closure."""
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import RelationError


class Poset:
    def __init__(self, covers, elements=None):
        self.covers = [tuple(c) for c in covers]
        elems = set(elements or [])
        for x, y in self.covers:
            elems.add(x)
            elems.add(y)
        self.elements = sorted(elems, key=str)
        self._le = self._closure()

    def _closure(self):
        le = {(x, x) for x in self.elements}
        le |= set(self.covers)
        changed = True
        while changed:
            changed = False
            for (a, b) in list(le):
                for (c, d) in list(le):
                    if b == c and (a, d) not in le:
                        le.add((a, d))
                        changed = True
        for (a, b) in le:
            if a != b and (b, a) in le:
                raise RelationError(
                    f"{a} <= {b} and {b} <= {a}: not a poset (antisymmetry fails)",
                    hint="cover data must not create a directed cycle")
        return le

    def leq(self, x, y):
        return (x, y) in self._le

    def hasse_quiver(self):
        names = {}
        arrows = {}
        for k, (x, y) in enumerate(self.covers):
            name = f"c{k}"
            names[name] = (x, y)
            arrows[name] = (x, y)
        Q = Quiver(list(self.elements), arrows)
        return Q, names
