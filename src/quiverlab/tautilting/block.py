"""The algebra-level ``tau_tilting`` compute-kind payload (Plan 45 / C4, Task H). Assembles
the exchange graph (pairs + g-matrices + brick-labelled Hasse edges), the wall-and-chamber
fan (n = 2, 3), the maximal-green-sequence count, and the AIR four-way counts -- everything
budget-capped with the honest complete-iff contract. Consumed by all three tiers (hpc.spec,
the Pyodide runner, the GUI) and the worked-steps report."""
from __future__ import annotations

# Citation KEYS (registry names); each runner resolves them to (key, formatted) pairs via
# its own _citation_pairs, exactly as every other algebra-level block does.
_CITATIONS = ["air_tau_tilting", "demonet_iyama_jasso", "king_stability"]


def tau_tilting_block(A, budget=512):
    """The ``tau_tilting`` block for algebra ``A``: the support tau-tilting pairs, their
    g-matrices and brick-labelled Hasse edges, the fan (n=2,3), the green-sequence count,
    and the four-way counts (all honest / None when budget-capped)."""
    from quiverlab.tautilting.green import maximal_green_sequences
    from quiverlab.tautilting.mutation import exchange_graph
    from quiverlab.tautilting.silting import silting_count
    from quiverlab.tautilting.stability import wall_and_chamber_fan
    from quiverlab.tautilting.torsion import (hasse_orientation, semibricks,
                                              torsion_class_data)
    verts = list(A.quiver.vertices)
    n = len(verts)
    eg = exchange_graph(A, budget_pairs=budget)
    block = {
        "kind": "tau_tilting",
        "n": n,
        "complete": eg.is_complete,
        "status": eg.status,
        "num_pairs": len(eg.vertices),
        "references": list(_CITATIONS),        # citation KEYS; each runner adds "citations"
    }
    block["pairs"] = [
        {"id": i, "g_matrix": rec["g_matrix"], "label": rec["label"],
         "summand_dimvecs": [{str(k): int(v) for k, v in dv.items()}
                             for dv in rec["summand_dimvecs"]],
         "support": list(rec["support"]), "is_initial": rec["is_initial"]}
        for i, rec in enumerate(eg.vertices)]
    if not eg.is_complete:
        # honest budget cap: the fan/counts/green would be a partial lie -- omit them.
        block["hasse"] = []
        block["fan"] = None
        block["green_count"] = None
        block["counts"] = None
        return block
    orient = hasse_orientation(eg)
    hasse = []
    for (i, j), d in orient.items():
        a, b = (i, j) if d == "down" else (j, i)      # a -> b downward
        hasse.append({"from": a, "to": b, "brick_dimvec": eg.arrows[(i, j)]["brick"]})
    block["hasse"] = hasse
    block["fan"] = wall_and_chamber_fan(A, budget=budget) if n in (2, 3) else None
    block["green_count"] = maximal_green_sequences(A, cap=budget)["count"]
    n_torsion = len({tuple(torsion_class_data(rec["pair"])["gen_dimvecs"])
                     for rec in eg.vertices})
    block["counts"] = {
        "pairs": len(eg.vertices),
        "torsion": n_torsion,
        "silting": silting_count(A, budget=budget)["count"],
        "semibricks": len(semibricks(A, budget=budget)),
    }
    return block
