"""Builder registry: name -> callable(*args) -> quiverlab Algebra (quiverlab.batch).

Adapted from hanlab's ``labdb`` builder registry (Marco Armenta's HansConjecture,
MIT): the same catalog of named families, but each name now resolves to the
corresponding quiverlab family constructor.  Every builder threads an optional
``field`` (default: the family's own default, exact CC) so ``scan.analyze`` can
rebuild an algebra over GF(p) for the exact F_p Hochschild computation.

The six core builders (mapped per the Plan-06 brief):
  truncated_polynomial(a)        -> families.basic.truncated_polynomial(a)
  quantum_ci(c)                  -> families.QuantumCI(q=c)
  cyclic_nakayama(n, ell)        -> families.NakayamaAlgebra(n=n, l=ell, cyclic=True)
  linear_path_algebra(n)         -> families.basic.linear_path_algebra(n)
  dynkin(typ, n)                 -> families.PathAlgebra(f"{typ}{n}")
  reduction_system(ngen, rules, name)
                                 -> families.zoo.build_from_record({...})   (open zone)
"""
from quiverlab.families import basic
from quiverlab.families.nakayama import NakayamaAlgebra
from quiverlab.families.path_algebra import PathAlgebra
from quiverlab.families.quantum import QuantumCI
from quiverlab.families import zoo as _zoo


def _truncated_polynomial(a, field=None):
    return basic.truncated_polynomial(a, field=field)


def _quantum_ci(c, field=None):
    return QuantumCI(q=c, field=field)


def _cyclic_nakayama(n, ell, field=None):
    return NakayamaAlgebra(n=n, l=ell, cyclic=True, field=field)


def _linear_path_algebra(n, field=None):
    return basic.linear_path_algebra(n, field=field)


def _dynkin(typ, n, field=None):
    return PathAlgebra(f"{typ}{n}", field=field)


def _reduction_system(ngen, rules, name, field=None):
    """Open-zone builder: materialise A = k<g>/I from a confluent reduction system,
    routed through quiverlab's zoo reduction-system -> Algebra path."""
    return _zoo.build_from_record(
        {"ngen": ngen, "rules": rules, "name": name, "dim": None}, field=field)


BUILDERS = {
    "truncated_polynomial": _truncated_polynomial,   # (a,)
    "quantum_ci": _quantum_ci,                        # (c,)
    "cyclic_nakayama": _cyclic_nakayama,              # (n, ell)
    "linear_path_algebra": _linear_path_algebra,      # (n,)
    "dynkin": _dynkin,                                # (typ, n)
    "reduction_system": _reduction_system,            # (ngen, rules, name) -- open zone
}


def build_algebra(spec):
    """Construct the quiverlab Algebra named by a spec dict (builder + args)."""
    return BUILDERS[spec["builder"]](*spec.get("args", []))
