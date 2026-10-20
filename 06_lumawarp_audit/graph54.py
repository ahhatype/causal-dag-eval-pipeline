"""Shared truth: the pipeline's analysis graph (36 observed variables, 63 edges) and node scales, read from the repository.
Paths: the repository root is two levels above this file; DATA_DIR (default <repo>/02_data) holds the replicates."""
import os, pandas as pd
P = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA_DIR = os.environ.get("DATA_DIR", os.path.join(P, "02_data"))
if not os.path.isabs(DATA_DIR): DATA_DIR = os.path.join(P, DATA_DIR.removeprefix("./"))
E = pd.read_csv(os.path.join(P, "outputs", "03_shape_validation", "analysis_graph_edges.csv"))
TRUE = {(a, b) for a, b in zip(E.parent, E.child)}; NODES = sorted(set(E.parent) | set(E.child))
PA = {v: {a for a, b in TRUE if b == v} for v in NODES}
SC = dict(zip(*[pd.read_csv(os.path.join(P, "01_data_generation", "config", "nodes.csv"))[c] for c in ("id", "scale")]))
CONT = {v for v in NODES if SC.get(v) == "continuous"}
UEDGES = {("nutrients", "bone_remodeling"), ("microbiome", "urine_chemistry"), ("hydration", "urine_concentration")}   # curved in the nonlinear variant
SAT = ("urine_chemistry", "mineralized_renal_material")   # saturating, monotone control
def topo():
    order = []; rem = {v: set(PA[v]) for v in NODES}
    while rem:
        ready = sorted(v for v, p in rem.items() if not p - set(order)); order += ready
        for v in ready: rem.pop(v)
    return order
TOPO = topo()
def replicate(rep, dataset, data_dir=None):
    return pd.read_csv(os.path.join(data_dir or DATA_DIR, "replicates", f"r{rep:03d}", f"{dataset}.csv"))[NODES]
def full_set(rep, data_dir=None): return replicate(rep, "full_set", data_dir)
