"""Skew-gentle algebras (Plan 68 / R32).

The user-facing object is the **triple** ``(Q, I, Sp)`` (``triple.py``): a gentle
pair ``(Q, I)`` plus a set ``Sp`` of *special* vertices, each carrying a specialised
loop ``eps`` with ``eps^2 = eps`` (an idempotent, He-Zhou-Zhu Def 1.3).  The ideal
``<eps^2 - eps>`` is NOT admissible, so the internal object is always the
**idempotent-split algebra** ``kQ_hat/I_hat`` (``split.py`` / Chen sec 3): each special
vertex splits into two copies ``i+, i-``, every ordinary arrow splits over its
endpoints' copies, the special loops become the vertex splitting, and ``I`` becomes
zero + commutative mesh relations.  The split algebra is admissible, ``dim``-certified
against the associated gentle algebra (HZZ Lemma 1.5), and characteristic-free, so
every existing engine runs on it verbatim.

Public surface:

- ``SkewGentleTriple`` / ``is_skew_gentle_triple`` / ``associated_gentle`` (``triple``)
- ``SkewGentleAlgebra`` / ``split_quiver`` (``split``)
- ``classify`` / ``skew_gentle_module`` / ``skew_gentle_indecomposables`` (``modules``)
- ``is_representation_finite`` / ``brick_finite_certificate`` /
  ``support_tau_tilting`` (``certificate``)
- ``skew_gentle_block`` (``block``)
"""
from quiverlab.skewgentle.block import skew_gentle_block
from quiverlab.skewgentle.certificate import (brick_finite_certificate,
                                              is_representation_finite,
                                              support_tau_tilting)
from quiverlab.skewgentle.modules import (SkewGentleString, classify,
                                          skew_gentle_indecomposables,
                                          skew_gentle_module)
from quiverlab.skewgentle.split import SkewGentleAlgebra, split_quiver
from quiverlab.skewgentle.triple import (SkewGentleTriple, associated_gentle,
                                         is_skew_gentle_triple)

__all__ = [
    "SkewGentleTriple",
    "associated_gentle",
    "is_skew_gentle_triple",
    "SkewGentleAlgebra",
    "split_quiver",
    "SkewGentleString",
    "classify",
    "skew_gentle_module",
    "skew_gentle_indecomposables",
    "is_representation_finite",
    "brick_finite_certificate",
    "support_tau_tilting",
    "skew_gentle_block",
]
