from __future__ import annotations

from typing import Callable, Sequence

import numpy as np
import pandas as pd
import shap

from .interventional import AttributionResult


def permutation_shap(predict: Callable[[np.ndarray], np.ndarray], explain: pd.DataFrame, background: pd.DataFrame,
                     features: Sequence[str], permutations: int, seed: int) -> AttributionResult:
    """Model-agnostic permutation SHAP on the probability scale; `permutations` orders as antithetic pairs."""
    features = list(features)
    masker = shap.maskers.Independent(background[features].to_numpy(dtype=float), max_samples=len(background))
    explainer = shap.PermutationExplainer(predict, masker, seed=seed)
    k = len(features)
    res = explainer(explain[features].to_numpy(dtype=float), max_evals=(permutations // 2) * (2 * k + 1), silent=True)
    values = pd.DataFrame(res.values, columns=features, index=explain.index)
    base = float(np.mean(res.base_values))
    eff = values.sum(axis=1).to_numpy() + np.asarray(res.base_values) - predict(explain[features].to_numpy(dtype=float))
    return AttributionResult(values, base, eff)
