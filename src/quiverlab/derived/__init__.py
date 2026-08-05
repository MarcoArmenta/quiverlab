"""The derived-category surface (Plan 43 / C8): reified hyper-Hom classes, the
derived AR translate ``tau_{D^b} = nu[-1]`` on perfect complexes, a tilting-complex
verifier with ``End(T)`` as a structure-constant algebra, and a necessary-condition
derived fingerprint. A thin exact-linear-algebra layer over P37 (End/ModuleHom),
P38 (Cartan/Coxeter), P39 (complexes/hyper-Hom) and P41 (Nakayama/corner-transpose);
no new math engine. Public surface only -- engine internals stay internal."""
from quiverlab.derived.homs import hyper_hom_basis

__all__ = ["hyper_hom_basis"]
