from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from sklearn.base import clone
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold


def default_library(seed: int) -> dict:
    return {
        "logistic": LogisticRegression(C=np.inf, max_iter=5000),
        "logistic_l2": LogisticRegression(C=1.0, max_iter=5000),
        "random_forest": RandomForestClassifier(n_estimators=200, min_samples_leaf=5, n_jobs=-1, random_state=seed),
        "gradient_boosting": HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, random_state=seed),
    }


def _log_loss(p: np.ndarray, y: np.ndarray) -> float:
    p = np.clip(p, 1e-9, 1 - 1e-9)
    return float(-np.mean(y * np.log(p) + (1 - y) * np.log(1 - p)))


class SuperLearner:
    """Convex combination of learners, weights chosen to minimize cross-validated log loss (van der Laan et al. 2007)."""

    def __init__(self, library: dict, folds: int, seed: int):
        self.library = library
        self.folds = folds
        self.seed = seed

    def fit(self, X: np.ndarray, y: np.ndarray) -> "SuperLearner":
        names = list(self.library)
        cv = StratifiedKFold(self.folds, shuffle=True, random_state=self.seed)
        z = np.zeros((len(y), len(names)))
        for train, test in cv.split(X, y):
            for j, name in enumerate(names):
                z[test, j] = clone(self.library[name]).fit(X[train], y[train]).predict_proba(X[test])[:, 1]
        k = len(names)
        res = minimize(lambda w: _log_loss(z @ w, y), np.full(k, 1 / k), method="SLSQP",
                       bounds=[(0, 1)] * k, constraints=({"type": "eq", "fun": lambda w: w.sum() - 1},))
        w = np.clip(res.x, 0, None)
        self.weights_ = dict(zip(names, w / w.sum()))
        self.cv_risk_ = {name: _log_loss(z[:, j], y) for j, name in enumerate(names)}
        self.fitted_ = {name: clone(self.library[name]).fit(X, y) for name in names if self.weights_[name] > 1e-6}
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Probability of the outcome."""
        X = np.asarray(X, dtype=float)
        return sum(self.weights_[n] * m.predict_proba(X)[:, 1] for n, m in self.fitted_.items())
