"""Domain->field shim for reduction_system_of (Pillar-4). There is no
A.reduction_system(); CS re-runs Plan-03's public build_reduction_system over the
algebra's ALREADY-BUILT domain. build_reduction_system needs a `field` with
parse_entry / make_domain; _DomainField makes both identity/return-dom so the
stored relation coefficients flow straight into A.domain (no re-parse, no second
domain). If build_reduction_system's field contract ever needs more, add exactly
the methods it calls (grep groebner/system.py: only parse_entry, make_domain)."""


class _DomainField:
    def __init__(self, dom):
        self._dom = dom

    def parse_entry(self, x):
        return x

    def make_domain(self, entries):
        return self._dom


def field_for_domain(dom):
    return _DomainField(dom)
