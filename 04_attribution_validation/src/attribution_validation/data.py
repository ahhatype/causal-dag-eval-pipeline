from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[3]
DG_CONFIG = REPO / "01_data_generation" / "config"
OUTCOME = "nephrolithiasis"


def env_dir(key: str, default: str) -> Path:
    """One directory setting from .env (nothing else is read), resolved against the repo."""
    value = default
    env = REPO / ".env"
    if env.exists():
        for line in env.read_text().splitlines():
            m = re.match(rf"^\s*{key}\s*=\s*(.*)$", line)
            if m:
                value = m.group(1).strip().strip("'\"")
                break
    path = Path(value)
    return path if path.is_absolute() else REPO / value.removeprefix("./")


@dataclass(frozen=True)
class GenerationConfig:
    nodes: pd.DataFrame
    edges: pd.DataFrame
    params: dict
    calibration: dict


def load_generation_config() -> GenerationConfig:
    return GenerationConfig(
        nodes=pd.read_csv(DG_CONFIG / "nodes.csv"),
        edges=pd.read_csv(DG_CONFIG / "edges.csv"),
        params=yaml.safe_load((DG_CONFIG / "params.yaml").read_text()),
        calibration=yaml.safe_load((DG_CONFIG / "calibration.yaml").read_text()),
    )


@dataclass(frozen=True)
class AnalysisGraph:
    nodes: list[str]
    edges: list[tuple[str, str]]
    scale: dict[str, str]

    @property
    def graph(self) -> nx.DiGraph:
        g = nx.DiGraph()
        g.add_nodes_from(self.nodes)
        g.add_edges_from(self.edges)
        return g

    def parents(self, node: str) -> list[str]:
        return [p for p, c in self.edges if c == node]


def analysis_graph(cfg: GenerationConfig) -> AnalysisGraph:
    """Observed variables only; era-fixed nodes collapse into mission era; edges from era-constant nodes dropped."""
    n = cfg.nodes
    constant = {k for k, v in cfg.params["era_values"].items() if v[0] == v[1]}
    observed = set(n.loc[n["observed"] & n["simulate"], "id"])
    collapsed = set(n.loc[(n["level"] == "era_fixed") & n["observed"], "id"])
    keep = observed - collapsed
    edges = set()
    for p, c in zip(cfg.edges["parent"], cfg.edges["child"]):
        if p not in observed or c not in observed or p in constant or c in constant:
            continue
        p2 = "mission_era" if p in collapsed else p
        c2 = "mission_era" if c in collapsed else c
        if p2 != c2:
            edges.add((p2, c2))
    g = nx.DiGraph()
    g.add_nodes_from(sorted(keep))
    g.add_edges_from(sorted(edges))
    order = list(nx.lexicographical_topological_sort(g))
    scale = dict(zip(n["id"], n["scale"]))
    return AnalysisGraph(nodes=order, edges=sorted(edges), scale={v: scale[v] for v in order})


def read_provenance(data_dir: Path) -> dict[str, str]:
    f = data_dir / "provenance.txt"
    if not f.exists():
        raise FileNotFoundError(f"no provenance.txt in {data_dir}: run `make data`")
    return dict(line.split(": ", 1) for line in f.read_text().splitlines() if ": " in line)


def check_data_current(data_dir: Path, cfg: GenerationConfig) -> int:
    prov = read_provenance(data_dir)
    if prov.get("config hash") != cfg.calibration["config_hash"]:
        raise RuntimeError(f"data in {data_dir} were generated from a different config: run `make data`")
    return int(prov["replicates"])


def load_replicate(data_dir: Path, replicate: int, dataset: str) -> pd.DataFrame:
    return pd.read_csv(data_dir / "replicates" / f"r{replicate:03d}" / f"{dataset}.csv")


def feature_sets(data_dir: Path) -> dict[str, list[str]]:
    f = pd.read_csv(data_dir / "feature_sets.csv")
    truthy = lambda s: s.astype(str).str.upper() == "TRUE"  # noqa: E731
    return {"ancestor": f.loc[truthy(f["in_ancestor_set"]), "id"].tolist(),
            "all": f.loc[truthy(f["in_all_features_set"]), "id"].tolist()}


def ground_truth(data_dir: Path, target: str) -> pd.Series:
    """True total effects (risk difference) for the source population ("full") or the selected one ("astronaut")."""
    g = pd.read_csv(data_dir / f"ground_truth_total_effects_{target}.csv")
    return g.set_index("feature")[f"effect_{target}"]


def ground_truth_se(data_dir: Path, target: str) -> pd.Series:
    """Paired Monte Carlo standard errors of the true total effects."""
    g = pd.read_csv(data_dir / f"ground_truth_total_effects_{target}.csv")
    return g.set_index("feature")[f"se_{target}"]


GROUND_TRUTH_TARGET = {"full_set": "full", "reference_subsample": "full", "astronaut_set": "astronaut"}


def true_outcome_probability(data: pd.DataFrame, cfg: GenerationConfig) -> np.ndarray:
    """P(outcome | its generating parents), from the frozen calibration: the bound on achievable discrimination."""
    e = cfg.edges[cfg.edges["child"] == OUTCOME]
    lin = cfg.calibration["logit"][OUTCOME] + sum(
        float(c) * data[p].to_numpy(dtype=float) for p, c in zip(e["parent"], e["coef"]))
    return 1.0 / (1.0 + np.exp(-lin))


def distance_to_outcome(ag: AnalysisGraph) -> dict[str, int]:
    """Shortest directed path length from each ancestor to the outcome."""
    rev = ag.graph.reverse()
    return {k: v for k, v in nx.single_source_shortest_path_length(rev, OUTCOME).items() if k != OUTCOME}
