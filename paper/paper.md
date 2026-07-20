---
title: "quiverlab: exact Hochschild theory for quivers with relations in Python"
tags:
  - Python
  - representation theory
  - quivers with relations
  - Hochschild cohomology
  - Gerstenhaber bracket
  - homological algebra
authors:
  - name: Marco Armenta
    orcid: 0000-0000-0000-0000        # TODO(Marco): replace with the real ORCID before submission
    corresponding: true
    email: drmarcoarmenta@gmail.com
    affiliation: 1
affiliations:
  - index: 1
    name: "Affiliation TODO(Marco)"   # TODO(Marco): set the real affiliation before submission
date: 18 July 2026
bibliography: paper.bib
---

# Summary

`quiverlab` is a pure-Python library for computing with finite-dimensional
associative algebras presented as **quivers with relations**, `A = kQ/I`, over
exact fields. Given a quiver and a list of relation strings, it certifies that the
algebra is finite-dimensional, builds an exact multiplication table, and computes
**Hochschild cohomology and homology** together with their algebraic operations —
the cup product, the Gerstenhaber bracket, and the cap action — as well as cyclic
homology, module Ext, and Cartan/Coxeter invariants. Every number is exact: the
library works over the rationals, over exact subfields of the complex numbers
(algebraic number fields `Q(α)`), and over every finite field `GF(p^n)`, and it
fails loudly on any floating-point input rather than returning an approximation.
It is aimed at research algebraists who barely program: three lines take a user
from a presentation to a certified Hochschild table, and every computation can emit
a human-readable worked-steps document.

# Statement of need

Hochschild cohomology `HH^\bullet(A)`, with its Gerstenhaber algebra structure, is
a central invariant in representation theory and deformation theory, yet it is
remarkably hard to compute by hand beyond the smallest examples. Researchers who
study quivers with relations — the standard presentation of finite-dimensional
algebras — have no installable tool that computes these invariants exactly and
without programming overhead. The existing options each stop short: they either do
not compute Hochschild cohomology at all, do not implement its operations, are not
exact, or require a substantial computer-algebra installation and scripting effort.
`quiverlab` fills this gap. It is `pip install quiverlab` with no external system
dependencies, exposes a flat, discoverable API, and returns certified exact answers
across characteristics — so an algebraist can, for instance, watch the
characteristic-`p` pathology of `k[x]/(x^2)` appear by changing a single argument.

# State of the field

The strongest existing system is **QPA** [@qpa], a mature GAP package for quivers
and path algebras: it constructs `kQ/I` by admissible ideals, computes minimal
module resolutions, module Ext, and Auslander–Reiten theory. QPA ships **no
Hochschild cohomology** — it must be assembled by hand via the enveloping algebra
and module Ext — and **no cup product or Gerstenhaber bracket**; installation
requires GAP, non-trivial for non-programmers, and historically awkward on Windows.
**SageMath** [@sagemath] provides only the *free* path algebra, with no
quotient-by-relations object, and an unreduced-bar Hochschild complex usable only at
toy sizes. **Magma**, **Macaulay2/Singular**, and **QuiverTools** address adjacent
problems (Ext algebras, noncommutative Gröbner bases, moduli of representations) but
none computes finite-dimensional Hochschild theory with its operations. On PyPI there
is **nothing**: no quivers-with-relations, no Hochschild cohomology, no bracket. In
short, no system on earth ships Hochschild cohomology *with its operations* for
finite-dimensional algebras, none ships the Chouhy–Solotar resolution, and nothing is
pip-installable for non-coders. `quiverlab` is built to be that system rather than a
contribution to any of the above, because its exact-only, non-coder-first design and
its resolution engine differ fundamentally from each existing architecture.

# Software design

`quiverlab` is layered. A `Domain` protocol carries all coefficient arithmetic
exactly (rationals with fraction-free elimination, number fields via `sympy`, finite
fields via a fast integer kernel with a pure-Python fallback), so no engine has a
floating-point code path; a static analysis gate forbids float literals in the source
tree. On top of it, a quiver-with-relations front end runs an exact noncommutative
**Gröbner (Buchberger–Mora overlap) completion** with a degree bound and an
**admissibility certificate**: it either returns a certified finite-dimensional
algebra with an irreducible-path basis, or fails loudly with the offending cycle —
never a hang and never a guess. The resolution layer offers four interchangeable
backends (normalized bar, minimal `A^e` syzygy, Bardzell for monomial algebras, and
the general **Chouhy–Solotar** resolution [@ChouhySolotar2015], the first full
implementation in any system, specializing exactly to Bardzell [@Bardzell1997] in the
monomial case). This general resolution is certified for quadratic-tip and all
monomial presentations; a non-quadratic non-monomial presentation raises at an
explicit boundary rather than risk a wrong answer. Deep Hochschild *dimensions* are
certified for closed-form families to arbitrary degree and, for general admissible
`kQ/I`, are computed and gated per instance by three independent checks (`d∘d = 0`, an
order condition, and degreewise agreement with the bar/minimal engines within their
window); the Gerstenhaber *operations* [@NegronWitherspoon2016; @Volkov2019] are
transported to bar cochains and are certified inside the bar-buildable degree window.
Module Ext uses minimal projective resolutions in the Green–Solberg–Zacharia style
[@GSZ2001]. An optional `pip install quiverlab[qpa]` backend runs
an **independent** recomputation of module Ext and Hochschild dimensions in QPA (via
the enveloping-algebra route) as a validation oracle in continuous integration.

# Research impact statement

`quiverlab` makes exact Hochschild computations routine where they were previously
manual, enabling systematic study of Gerstenhaber structure across families of
algebras and across characteristics — for example reproducing published Hochschild
dimensions of quantum complete intersections [@BGMS2005] and the vanishing behaviour
of hereditary algebras [@Happel1989], then extending them along parameter sweeps that
would be impractical by hand. Because it is exact and pip-installable, it lowers the
barrier to reproducible experiments in representation theory and provides a common,
citable reference implementation of the Chouhy–Solotar resolution against which future
work can be checked. The library unifies and generalizes the author's prior research
software into a tool independent of the application that produced it.

# AI usage disclosure

The library's design and mathematics are the author's. AI coding assistants were used
under close human supervision for implementation scaffolding, test authoring, and
documentation drafting; every mathematical claim, algorithm, and certified value was
verified by the author against hand computations, published literature, and an
independent computer-algebra cross-check (QPA). No result reported by the software or
this paper rests on unverified AI output.

# Acknowledgements

We thank the developers of QPA [@qpa], GAP [@gap4], SymPy, and NumPy, on whose ideas
and tools this work builds.
<!-- TODO(Marco): add financial support / funding acknowledgement before submission -->

# References
