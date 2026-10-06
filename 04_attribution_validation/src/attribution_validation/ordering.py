from __future__ import annotations

from typing import Callable, Sequence

import numpy as np
import pandas as pd

from .data import AnalysisGraph
from .interventional import AttributionResult, OrderSampler, ancestral_pairs


def ordering_only_shap(predict: Callable[[np.ndarray], np.ndarray], explain: pd.DataFrame, background: pd.DataFrame,
                       features: Sequence[str], ag: AnalysisGraph, orders: int, burn_in: int, thin: int,
                       seed: int) -> AttributionResult:
    """Asymmetric Shapley values: graph-consistent orders, features outside a coalition taken from the background.

    Differs from interventional SHAP only in the value function; the orders are sampled the same way.
    """
    features = list(features)
    rng = np.random.default_rng(seed)
    n_eval, n_bg = len(explain), len(background)
    bg = background[features].to_numpy(dtype=float)
    baseline = float(np.mean(predict(bg)))
    xv = {f: np.repeat(explain[f].to_numpy(dtype=float), n_bg) for f in features}
    col = {f: j for j, f in enumerate(features)}
    sampler = OrderSampler(features, ancestral_pairs(ag, features), rng, burn_in, thin)
    out = np.zeros((n_eval, len(features)))
    for _ in range(orders):
        batch = np.tile(bg, (n_eval, 1))
        prev = np.full(n_eval, baseline)
        for f in sampler.sample():
            batch[:, col[f]] = xv[f]
            cur = predict(batch).reshape(n_eval, n_bg).mean(axis=1)
            out[:, col[f]] += cur - prev
            prev = cur
    out /= orders
    eff = out.sum(axis=1) + baseline - predict(explain[features].to_numpy(dtype=float))
    return AttributionResult(pd.DataFrame(out, columns=features, index=explain.index), baseline, eff)
