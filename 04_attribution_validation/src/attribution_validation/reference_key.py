"""Reference credits of the ConfoundingSHAP game on the true model.

The true propensity pi(x) and the potential-outcome probabilities mu_a(x) = P(outcome | do(drug = a)) are known for
every simulated record, so with X_S a covariate subset

    delta_S(x_S) = E[pi mu_1 | x_S] / E[pi | x_S] - E[(1 - pi) mu_0 | x_S] / E[1 - pi | x_S]

is a (propensity-weighted) smoothing of noise-free quantities, and the signed global game is v(S) = ATE - E[delta_S(X_S)] with
ATE = E[mu_1 - mu_0]. Shapley credits are estimated from random permutations (antithetic pairs), sharing the
coalition values between permutations. The smoother is gradient boosting; its error shows in v(all), which is
zero for the exact game when the covariates are sufficient.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.ensemble import HistGradientBoostingRegressor

from .truesim import TrueSCM


def potential_outcomes(scm: TrueSCM, u: dict, drug: str, n: int) -> tuple[np.ndarray, np.ndarray]:
    mu1 = scm.outcome_probability(scm.simulate(u, {drug: np.ones(n)}))
    mu0 = scm.outcome_probability(scm.simulate(u, {drug: np.zeros(n)}))
    return mu1, mu0


def _delta(X: np.ndarray, mu1: np.ndarray, mu0: np.ndarray, pi: np.ndarray, seed: int) -> float:
    """mean over records of delta_S(x_S) = m_1(x_S) - m_0(x_S), the propensity-weighted mean of mu_1 given x_S minus the
    (1 - pi)-weighted mean of mu_0. Weighted regressions avoid dividing two smoothed quantities, which explodes where
    the smoothed propensity is near zero."""
    if X.shape[1] == 0:
        return float(np.average(mu1, weights=pi) - np.average(mu0, weights=1 - pi))

    def fit(y: np.ndarray, w: np.ndarray) -> np.ndarray:
        m = HistGradientBoostingRegressor(max_iter=60, learning_rate=0.15, max_leaf_nodes=15, min_samples_leaf=100,
                                          random_state=seed)
        return np.clip(m.fit(X, y, sample_weight=w).predict(X), 0.0, 1.0)

    return float(np.mean(fit(mu1, pi) - fit(mu0, 1 - pi)))


def reference_credits(scm: TrueSCM, pop: pd.DataFrame, u: dict, drug: str, covariates: list[str], mask: np.ndarray,
                      permutations: int, seed: int, jobs: int = 7) -> dict:
    n = len(pop)
    mu1, mu0 = potential_outcomes(scm, u, drug, n)
    pi = scm.propensity(drug, {k: pop[k].to_numpy() for k in pop.columns if not k.startswith("U_")})
    mu1, mu0, pi = mu1[mask], mu0[mask], pi[mask]
    Xall = pop.loc[mask, covariates].to_numpy(dtype=float)
    ate = float((mu1 - mu0).mean())
    p = len(covariates)
    rng = np.random.default_rng(seed)
    orders = []
    for _ in range(permutations // 2):
        o = rng.permutation(p)
        orders += [o, o[::-1]]
    subsets = {frozenset(): None}
    for o in orders:
        for k in range(1, p + 1):
            subsets[frozenset(int(j) for j in o[:k])] = None
    keys = list(subsets)
    vals = Parallel(n_jobs=jobs)(delayed(_delta)(Xall[:, sorted(k)], mu1, mu0, pi, seed) for k in keys)
    v = {k: ate - d for k, d in zip(keys, vals)}
    credits = np.zeros(p)
    for o in orders:
        prev = v[frozenset()]
        seen: set[int] = set()
        for j in o:
            seen.add(int(j))
            cur = v[frozenset(seen)]
            credits[j] += cur - prev
            prev = cur
    credits /= len(orders)
    full = v[frozenset(range(p))]
    return {"credits": dict(zip(covariates, credits.tolist())), "ate": ate, "v_empty": v[frozenset()], "v_full": full,
            "crude": float(np.average(mu1, weights=pi) - np.average(mu0, weights=1 - pi)),
            "n": int(mask.sum()), "subsets": len(keys)}
