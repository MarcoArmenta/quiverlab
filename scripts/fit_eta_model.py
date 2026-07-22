"""Refit the GUI eta model (docs/gui/runner.py ETA_MODEL). Run manually:

    .venv/bin/python scripts/fit_eta_model.py

Blocks numba (Pyodide has none), benchmarks both HH routes on loop algebras
k[x]/(x^m) for m in 3..6 plus the scalar invariants, fits alpha/p per route by
log-space least squares on points >= 0.3 s, and prints a ready-to-paste
ETA_MODEL literal. Paste it into runner.py ONLY when it drifts >2x from the
baked values (bucket widths absorb less)."""
import math
import sys
import time

sys.modules["numba"] = None          # force the pure-Python kernels
import quiverlab as ql               # noqa: E402

ql.verbose = False
CAP = 4_000_000


def cells(m, n):
    return (m * (m - 1) ** (n + 1)) * (m * (m - 1) ** n)


def timed(f):
    s = time.monotonic()
    f()
    return time.monotonic() - s


def fit(points):
    best = None
    for p in [x / 100 for x in range(50, 201, 5)]:
        la = sum(math.log(s / sum(cells(m, n) ** p for n in range(top + 1)))
                 for m, top, s in points) / len(points)
        err = sum((math.log(s) - la - math.log(
            sum(cells(m, n) ** p for n in range(top + 1)))) ** 2
            for m, top, s in points)
        if best is None or err < best[0]:
            best = (err, p, math.exp(la))
    return best


def main():
    grid = {"bar": [], "fast": []}
    for m in (3, 4, 5, 6):
        for route, field in (("bar", ql.CC), ("fast", ql.GF(2))):
            A = ql.Quiver(vertices=[1], arrows={"x": (1, 1)}).algebra(
                relations=["*".join(["x"] * m)], field=field)
            for top in range(2, 11):
                if any(cells(m, n) > CAP for n in range(top + 1)):
                    break
                s = timed(lambda: A.hochschild_cohomology(top))
                grid[route].append((m, top, s))
                if s > 25:
                    break
    print("ETA_MODEL = {")
    for route in ("bar", "fast"):
        pts = [x for x in grid[route] if x[2] >= 0.3]
        err, p, alpha = fit(pts)
        worst = max(max(alpha * sum(cells(m, n) ** p for n in range(top + 1)) / s,
                        s / (alpha * sum(cells(m, n) ** p for n in range(top + 1))))
                    for m, top, s in pts)
        print('    "%s":  {"alpha": %.4e, "p": %s},   # worst-off %.2fx on %d pts'
              % (route, alpha, p, worst, len(pts)))
    A9 = next(iter(ql.zoo(dim_max=12)))
    A3 = ql.Quiver(vertices=[1], arrows={"x": (1, 1)}).algebra(
        relations=["x*x*x"], field=ql.CC)
    sc = {"cartan": 0.0, "coxeter_polynomial": 0.0, "center": 0.0,
          "global_dimension": 0.0}
    for A in (A3, A9):
        sc["cartan"] = max(sc["cartan"], timed(A.cartan_matrix))
        sc["coxeter_polynomial"] = max(sc["coxeter_polynomial"],
                                       timed(A.coxeter_polynomial))
        sc["center"] = max(sc["center"], timed(A.center))
        sc["global_dimension"] = max(sc["global_dimension"],
                                     timed(A.global_dimension))
    print('    "scalars": {"cartan": %.2g, "coxeter_polynomial": %.2g,\n'
          '                "center": %.2g, "global_dimension": %.2g},'
          % (max(sc["cartan"], 0.01), max(sc["coxeter_polynomial"], 0.2),
             max(sc["center"], 0.05), max(sc["global_dimension"], 0.5)))
    print("}")


if __name__ == "__main__":
    main()
