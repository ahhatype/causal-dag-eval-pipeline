"""Build fig2_portrait.json: the editor layout turned to portrait, in points, with labels placed by search.

Usage: python portrait_layout.py <fig2_layout.json> <edges.csv> <out.json> [pitch]
"""
import json
import sys

import pandas as pd

import dag_labels as dl

FS = 6.5
WIDTH = 288.0
L, R, TOP = 74.0, 214.0, 14.0
NR, NR_OUT = 2.6, 3.6

BREAKS = {
    "cumulative|mission duration": "cumulative|mission|duration",
    "detect long term|health outcomes": "detect long|term health|outcomes",
    "long term|health outcomes": "long term|health|outcomes",
    "loss of mission|objectives": "loss of|mission|objectives",
    "mineralized|renal material": "mineralized|renal|material",
}

FIXED = {  # hand-set placements where the search leaves a label ambiguous
    "individual_readiness": (("L", 0), "INDIVIDUAL|READINESS"),
    "crew_capability": (("R", 0), "CREW|CAPABILITY"),
    "long_term_health_outcomes": (("R", 0), "LONG TERM|HEALTH|OUTCOMES"),
    "detect_long_term_health_outcomes": (("R", 0), "DETECT LONG|TERM HEALTH|OUTCOMES"),
    "pre_flight_fitness": (("R", 0), "PRE-FLIGHT|FITNESS"),
    "percutaneous_nephrostomy": (("L", 0), "PERCUTANEOUS|NEPHROSTOMY"),
    "resistive_exercise": (("R", 0), "RESISTIVE|EXERCISE"),
    "urine_concentration": (("L", 0), "URINE|CONCENTRATION"),
    "bisphosphonates": (("R", 0), "BISPHOSPHONATES"),
    "urine_chemistry": (("L", 0), "URINE CHEMISTRY"),
}

NUDGE = {"percutaneous_nephrostomy": (0, 16), "bisphosphonates": (0, 16)}  # off the crowded third rank


def main(layout_path, edges_path, out_path, pitch=22.0):
    lay = json.load(open(layout_path))
    lay = lay.get("data", lay)
    nodes = lay["nodes"]
    edges = list(pd.read_csv(edges_path)[["parent", "child"]].itertuples(index=False, name=None))
    xs = sorted({v["x"] for v in nodes.values()})
    count = [sum(v["x"] == x for v in nodes.values()) for x in xs]
    ys = [TOP]
    for i in range(len(xs) - 1):
        ys.append(ys[-1] + pitch + 5 * max(0, max(count[i], count[i + 1]) - 2))
    P = {k: (L + (v["y"] - 120) / 800 * (R - L), ys[xs.index(v["x"])]) for k, v in nodes.items()}
    for k, (dx, dy) in NUDGE.items():
        P[k] = (P[k][0] + dx, P[k][1] + dy)
    height = ys[-1] + 14
    rad = {k: (NR_OUT if k == "nephrolithiasis" else NR) for k in nodes}
    options = {k: dl.variants(v["label"], [t.upper() for t in BREAKS.get(v["label"], v["label"]).split("|")], 72, FS)
               for k, v in nodes.items()}
    fixed = {k: (cd, ls.split("|")) for k, (cd, ls) in FIXED.items()}
    lab = dl.place(P, rad, list(nodes), options, [(P[a], P[b]) for a, b in edges], WIDTH, height, FS, fixed, gap=2.2)
    out = {"units": "pt", "width": WIDTH, "height": height,
           "font": {"size": FS, "tracking": dl.TRACK, "leading": dl.LEAD, "weight": dl.WEIGHT},
           "node_r": NR, "outcome_r": NR_OUT, "nodes": {}}
    for k in nodes:
        a = lab[k]
        out["nodes"][k] = {"x": round(P[k][0], 1), "y": round(P[k][1], 1), "label": "|".join(a["lines"]),
                           "anchor": a["anchor"], "lx": round(a["lx"], 1), "ly": round(a["ly"], 1), "place": a["place"]}
    json.dump(out, open(out_path, "w"), indent=1)
    for k, a in lab.items():
        if a["residual"] > 1:
            print(f"  {k}: {a['place']} {a['residual']}")


if __name__ == "__main__":
    main(*sys.argv[1:4], *(float(a) for a in sys.argv[4:5]))
