"""quiverlab.viz: draw() and tikz() for quivers with relations (spec §5 c.10,
§3.7). A layered layout is computed ONCE in exact int/Fraction coordinates
(layout.py) and rendered either with matplotlib (draw.py) or as TikZ (tikz.py)
from the identical coordinates.

FLOAT POLICY: this package writes NO float/complex literal and NO float() call.
All layout arithmetic is int/Fraction; matplotlib coerces to float internally
(outside src/), so tests/test_no_floats.py covers viz with no exemption (D3).

BACKEND POLICY: importing this package does NOT mutate matplotlib's global backend
(no matplotlib.use(...) call) -- that would hijack a user's interactive session.
draw.py builds a Figure with its own Agg canvas for file export and returns it; a
user who wants an interactive figure gets the returned Figure and shows it with
their own pyplot, untouched by us."""

__all__ = []  # draw/tikz are Algebra methods; layout() is imported directly when needed
