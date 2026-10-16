"""Write node_layout.csv: a layered left-to-right layout of the analysis graph.

Columns are longest-path layers from the roots ("mission era" is pinned to its own column); order within
a column comes from barycentre sweeps to reduce crossings. The figure script reads only node_layout.csv,
so hand-edit x and y there to fine-tune; rerunning this script resets the layout.
"""
import csv
from collections import defaultdict

from style import ROOT

edges = list(csv.DictReader(open(ROOT / "outputs" / "03_shape_validation" / "analysis_graph_edges.csv")))
roles = {r["id"]: r for r in csv.DictReader(open(ROOT / "02_data" / "feature_sets.csv"))}
ERA = "mission_era"
parents, children = defaultdict(set), defaultdict(set)
for e in edges:
    parents[e["child"]].add(e["parent"])
    children[e["parent"]].add(e["child"])
nodes = sorted(set(parents) | set(children))

layer = {ERA: 0}


def depth(n):
    if n not in layer:
        layer[n] = max((depth(p) + 1 for p in parents[n] if p != ERA), default=1)
    return layer[n]


for n in nodes:
    depth(n)
for _ in range(3):  # pull roots right, next to their earliest child
    for n in nodes:
        if n != ERA and not (parents[n] - {ERA}) and children[n]:
            layer[n] = min(layer[c] for c in children[n]) - 1
order = defaultdict(list)
for n in nodes:
    order[layer[n]].append(n)


def idx():
    return {n: i / max(len(v) - 1, 1) for v in order.values() for i, n in enumerate(v)}


for sweep in range(24):
    pos = idx()
    down = sweep % 2 == 0
    for k in sorted(order, reverse=not down):
        def bary(n):
            nb = [m for m in (parents[n] if down else children[n]) if m != ERA]
            return sum(pos[m] for m in nb) / len(nb) if nb else pos[n]
        order[k].sort(key=bary)

nmax = max(len(v) for v in order.values())
with open("node_layout.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["id", "label", "x", "y", "group"])
    for k in sorted(order):
        v = order[k]
        for i, n in enumerate(v):
            y = (nmax - 1) / 2 - i * (nmax - 1) / (len(v) - 1) if len(v) > 1 else 0
            group = "outcome" if n == "nephrolithiasis" else (
                "ancestor" if roles[n]["in_ancestor_set"] == "TRUE" else "non_ancestor")
            w.writerow([n, roles[n]["name"] if n in roles else "nephrolithiasis", k, round(y, 3), group])
print({k: order[k] for k in sorted(order)})
