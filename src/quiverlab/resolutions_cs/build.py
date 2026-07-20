"""CS entry currency: re-derive the Plan-03 ReductionSystem from a stored Algebra.

reduction_system_of(A) re-runs the public groebner.build_reduction_system over
A's own domain (via the _fieldshim), because a groebner-lowered Algebra keeps its
quiver + parsed relations but NOT a reduction system (spec Pillar-4). The rs's
domain therefore IS A.domain (same instance), so CS's matrices and A's structure
constants live in one field."""


def reduction_system_of(A):
    from quiverlab.groebner import build_reduction_system
    from quiverlab.resolutions_cs._fieldshim import field_for_domain
    if A.quiver is None or A.relations is None:
        raise ValueError("CS needs an algebra built by Quiver.algebra (quiver+relations present)")
    return build_reduction_system(A.quiver, list(A.relations), field_for_domain(A.domain))
