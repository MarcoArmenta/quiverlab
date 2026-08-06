"""Euler bilinear form, Tits quadratic form, definiteness classification
(Plan 38 / C2). All exact: sympy Rational off the integer Cartan matrix.

The Euler form matrix is E = C^{-1} (with cartan_matrix's convention
``C[i][j]`` = dim e_i A e_j = #paths i->j), so <d, e> = d E e^T; for finite
global dimension <dim M, dim N> = sum (-1)^i dim Ext^i(M, N). On a
hereditary path algebra this is exactly the arrow formula
<d, e> = sum_v d_v e_v - sum_{a: s->t} d_s e_t (E = I - N, N the arrow
adjacency), and it reproduces dim Hom(M, N) - dim Ext^1(M, N) with THIS
package's ext convention (verified: Ext^1(S_i, S_j) counts arrows i->j).

The Tits form is q(d) = <d, d>; its symmetrized matrix is E + E^T. The
finite/tame/wild reading of definiteness is a THEOREM only for hereditary
algebras (Gabriel; Donovan-Freislich/Nazarova) -- form_type computes the
signature for any invertible Cartan and the docstring says exactly this.
"""
from __future__ import annotations

import sympy as sp

from quiverlab.errors import QuiverlabError
from quiverlab.invariants.cartan import cartan_matrix


def euler_form_matrix(A) -> sp.Matrix:
    """The Euler form matrix E = C^{-1} (exact Rational). Raises loudly when the
    Cartan matrix is not unimodular (|det C| != 1): the homological Euler form
    <M, N> = sum (-1)^i dim Ext^i needs finite global dimension, whose necessary
    condition is a unimodular Cartan. A singular Cartan (det C = 0) and a
    non-unimodular one (e.g. k[x]/(x^2), det C = 2, gl.dim = infinity) are both
    refused -- in the latter C^{-1} has non-integer entries and is not the true
    Euler form."""
    C = sp.Matrix(cartan_matrix(A))
    det = C.det()
    if det not in (1, -1):
        reason = ("singular (det C = 0)" if det == 0
                  else f"not unimodular (det C = {det})")
        raise QuiverlabError(
            "Euler form needs a unimodular Cartan matrix (a necessary condition "
            f"for finite global dimension); the Cartan matrix here is {reason}.",
            hint="an algebra of infinite global dimension (e.g. k[x]/(x^2)) has "
                 "no homological Euler form")
    return C.inv()


def euler_form(A, d, e):
    """<d, e> = d E e^T for integer dimension vectors d, e (vertex order)."""
    E = euler_form_matrix(A)
    dv = sp.Matrix(1, E.rows, [sp.Integer(x) for x in d])
    ev = sp.Matrix(E.rows, 1, [sp.Integer(x) for x in e])
    return sp.nsimplify((dv * E * ev)[0, 0])


def tits_matrix(A) -> sp.Matrix:
    """The symmetrized (Tits) matrix E + E^T of the Euler form."""
    E = euler_form_matrix(A)
    return E + E.T


def tits_form(A, d):
    """The Tits quadratic form q(d) = <d, d>."""
    return euler_form(A, d, d)


def form_type(A) -> str:
    """'finite' / 'tame' / 'wild' by exact definiteness of the Tits form
    (positive definite / positive semidefinite-not-definite / neither). The
    representation-type meaning (finite/tame/wild representation type) is a
    theorem for HEREDITARY algebras only (Gabriel; Donovan-Freislich/Nazarova);
    for others this is the form's signature, computed exactly over Rationals."""
    Q = tits_matrix(A)
    if Q.is_positive_definite:
        return "finite"
    if Q.is_positive_semidefinite:
        return "tame"
    return "wild"
