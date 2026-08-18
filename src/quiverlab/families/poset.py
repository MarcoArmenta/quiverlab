"""Finite posets from cover data (spec §3.4). covers = [(x, y), ...] means x is
covered by y (Hasse edge, arrow x->y). Order <= is the reflexive-transitive closure."""
from quiverlab.combinat.quiver import Quiver
from quiverlab.errors import RelationError


class Poset:
    def __init__(self, covers, elements=None):
        self.covers = [tuple(c) for c in covers]
        for x, y in self.covers:
            if x == y:
                raise RelationError(
                    f"{x} covers {x}: not a poset (self-cover, reflexivity is implicit)",
                    hint="cover data must not contain a loop (x, x)")
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

    def order_complex(self):
        """The order complex (nerve) ``Delta(P)``: the abstract simplicial complex whose
        ``p``-simplices are the CHAINS ``x_0 < x_1 < ... < x_p`` of ``P``. ``|Delta(P)|``
        is the classifying space ``BP``, and by Cibils 1989 (generalizing
        Gerstenhaber-Schack 1983 from face posets to arbitrary finite posets)
        ``HH^n(kP) = H^n(Delta(P); k)`` for every ``n`` and every field ``k`` (Plan 75 /
        R9). :mod:`quiverlab.hochschild.simplicial` owns the class and the engine; this
        is the delegate."""
        from quiverlab.hochschild.simplicial import OrderComplex
        return OrderComplex.of(self)

    def is_bounded(self):
        """True iff ``P`` has BOTH a global minimum and a global maximum (``0^`` and
        ``1^``). Either one alone already makes the order complex a CONE, hence
        contractible, hence ``HH^{>=1}(kP) = 0`` -- see :meth:`has_global_bound`."""
        return self._global_min() is not None and self._global_max() is not None

    def _global_min(self):
        for x in self.elements:
            if all(self.leq(x, y) for y in self.elements):
                return x
        return None

    def _global_max(self):
        for x in self.elements:
            if all(self.leq(y, x) for y in self.elements):
                return x
        return None

    def has_global_bound(self):
        """True iff ``P`` has a global minimum OR a global maximum. EITHER makes
        ``Delta(P)`` a cone (contractible), so ``HH^{>=1}(kP) = 0`` -- the theorem
        instance behind the Plan-33 Boolean-lattice ``B_3`` pin. Note this is WEAKER
        than :meth:`is_bounded`, and it is the weaker condition that the contractibility
        argument actually needs."""
        return self._global_min() is not None or self._global_max() is not None

    def is_lattice(self):
        """True iff every pair has a unique least upper bound AND a unique greatest
        lower bound (an exact finite check over the order relation). Crowns are the
        standard NON-lattice posets, and their order complexes are circles -- which is
        why they carry the nonvanishing ``H^1`` this plan pins."""
        for x in self.elements:
            for y in self.elements:
                if self._unique_bound(x, y, upper=True) is None:
                    return False
                if self._unique_bound(x, y, upper=False) is None:
                    return False
        return True

    def _unique_bound(self, x, y, *, upper):
        if upper:
            bounds = [z for z in self.elements if self.leq(x, z) and self.leq(y, z)]
            best = [z for z in bounds if all(self.leq(z, w) for w in bounds)]
        else:
            bounds = [z for z in self.elements if self.leq(z, x) and self.leq(z, y)]
            best = [z for z in bounds if all(self.leq(w, z) for w in bounds)]
        return best[0] if len(best) == 1 else None

    def hasse_quiver(self):
        names = {}
        arrows = {}
        for k, (x, y) in enumerate(self.covers):
            name = f"c{k}"
            names[name] = (x, y)
            arrows[name] = (x, y)
        Q = Quiver(list(self.elements), arrows)
        return Q, names
