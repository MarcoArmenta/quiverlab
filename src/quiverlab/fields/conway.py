"""Bundled modulus polynomials for GF(p^n), little-endian coefficient lists
[c0, ..., cn] with cn = 1 (monic). Sourced from Lübeck's Conway-polynomial
tables; each entry is machine-validated irreducible at field construction and
irreducible+primitive in the test suite — a failing entry must be replaced
from the table, never trusted. Exact Conway normalization (cross-system
embedding compatibility) is a non-goal in this phase."""

CONWAY = {
    (2, 2): [1, 1, 1],
    (2, 3): [1, 1, 0, 1],
    (2, 4): [1, 1, 0, 0, 1],
    (2, 5): [1, 0, 1, 0, 0, 1],
    (2, 6): [1, 1, 0, 1, 1, 0, 1],
    (2, 7): [1, 1, 0, 0, 0, 0, 0, 1],
    (2, 8): [1, 0, 1, 1, 1, 0, 0, 0, 1],
    (2, 9): [1, 0, 0, 0, 1, 0, 0, 0, 0, 1],
    (2, 10): [1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1],
    (3, 2): [2, 2, 1],
    (3, 3): [1, 2, 0, 1],
    (3, 4): [2, 0, 0, 2, 1],
    (3, 5): [1, 2, 0, 0, 0, 1],
    (5, 2): [2, 4, 1],
    (5, 3): [3, 3, 0, 1],
    (5, 4): [2, 4, 4, 0, 1],
    (7, 2): [3, 6, 1],
    (7, 3): [4, 0, 6, 1],
    (11, 2): [2, 7, 1],
    (13, 2): [2, 12, 1],
    (17, 2): [3, 16, 1],
    (19, 2): [2, 18, 1],
    (23, 2): [5, 21, 1],
}
