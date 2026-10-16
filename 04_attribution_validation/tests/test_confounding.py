import numpy as np
import pandas as pd
import pytest
from sklearn.linear_model import LinearRegression

from attribution_validation import data as D
from attribution_validation.confounding import (common_cause_sets, confounder_mass, confounder_recovery,
                                                confounding_credits, pretreatment_covariates)
from attribution_validation.truesim import TrueSCM

CFG = D.load_generation_config()
AG = D.analysis_graph(CFG)
DATA_DIR = D.env_dir("DATA_DIR", "./02_data")


def linear_learner(device="cpu", n_estimators=1):
    return LinearRegression()


def toy(n: int, seed: int = 3) -> pd.DataFrame:
    """The collaborators' tutorial DAG: C confounds, Z affects the treatment only, P the outcome only; effect 0.5."""
    rng = np.random.default_rng(seed)
    c, z, p = (rng.integers(0, 2, n) for _ in range(3))
    a = (rng.random(n) < 0.1 + 0.3 * c + 0.3 * z).astype(int)
    y = c + 0.7 * p + 0.5 * a + rng.normal(0, 0.5, n)
    return pd.DataFrame({"C": c, "Z": z, "P": p, "A": a, "nephrolithiasis": y})


def test_metrics():
    credits = {"a": 0.5, "b": -0.3, "c": 0.1, "d": 0.1}
    assert confounder_mass(credits, ["a", "b"]) == pytest.approx(0.8 / 1.0)
    assert confounder_recovery(credits, ["a", "b"]) == 1.0
    assert confounder_recovery(credits, ["a", "c"]) == pytest.approx(0.5)      # top-2 by |credit| is a, b
    assert np.isnan(confounder_mass({"a": 0.0}, ["a"]))


def test_covariates_exclude_treatment_descendants_and_other_drugs():
    anc = D.feature_sets(DATA_DIR)["ancestor"] if (DATA_DIR / "feature_sets.csv").exists() else None
    if anc is None:
        pytest.skip("generated data not present")
    drugs = ["k_citrate", "bisphosphonates"]
    k = pretreatment_covariates(AG, anc, "k_citrate", ["thiazides", *drugs])
    b = pretreatment_covariates(AG, anc, "bisphosphonates", ["thiazides", *drugs])
    assert len(k) == 11 and len(b) == 10
    assert not {"urine_chemistry", "mineralized_renal_material", "thiazides", "bisphosphonates", "k_citrate"} & set(k)
    assert "bone_remodeling" in k and "bone_remodeling" not in b          # a descendant of bisphosphonates only
    assert set(common_cause_sets(AG, k, "k_citrate")["direct"]) == {"individual_factors", "urine_concentration", "mission_era"}
    assert set(common_cause_sets(AG, b, "bisphosphonates")["direct"]) == {"cumulative_mission_duration", "pre_flight_fitness",
                                                                          "mission_era"}
    assert {"hydration", "water_intake"} <= set(common_cause_sets(AG, k, "k_citrate")["all"])


def test_exact_game_reproduces_the_tutorial_answer_key():
    """Tutorial key (exact): C 0.33201, Z -0.01951, P 0; they sum to the crude-minus-ATE gap 0.3125."""
    df = toy(40000)
    r = confounding_credits(df, ["C", "Z", "P"], "A", linear_learner, budget=8, seed=1)
    assert r["exact"]
    assert r["credits"]["C"] == pytest.approx(0.33201, abs=0.03)
    assert r["credits"]["Z"] == pytest.approx(-0.01951, abs=0.03)
    assert abs(r["credits"]["P"]) < 0.03
    assert abs(sum(r["credits"].values()) - (r["v_full"] - r["v_empty"])) < 1e-9
    assert r["v_full"] == pytest.approx(0, abs=1e-9)
    assert r["adjusted_ate"] == pytest.approx(0.5, abs=0.03)


def test_randomized_treatment_gives_near_zero_credits():
    df = toy(40000, seed=5)
    rng = np.random.default_rng(9)
    df["A"] = (rng.random(len(df)) < 0.4).astype(int)
    df["nephrolithiasis"] = df.C + 0.7 * df.P + 0.5 * df.A + rng.normal(0, 0.5, len(df))
    r = confounding_credits(df, ["C", "Z", "P"], "A", linear_learner, budget=8, seed=1)
    assert max(abs(v) for v in r["credits"].values()) < 0.03


def test_approximator_path_is_efficient_and_picks_the_confounder():
    df = toy(20000, seed=7)
    for j in range(5):                                  # extra noise covariates make 2^p exceed the budget
        df[f"N{j}"] = np.random.default_rng(j).integers(0, 2, len(df))
    cov = ["C", "Z", "P"] + [f"N{j}" for j in range(5)]
    r = confounding_credits(df, cov, "A", linear_learner, budget=96, seed=2)
    assert not r["exact"] and r["coalitions_fitted"] <= 96
    assert abs(sum(r["credits"].values()) - (r["v_full"] - r["v_empty"])) < 1e-6
    assert max(r["credits"], key=lambda f: abs(r["credits"][f])) == "C"


def test_true_propensity_matches_realized_uptake():
    scm = TrueSCM(CFG)
    pop = scm.population(30000, np.random.default_rng(4))
    for drug in ("k_citrate", "bisphosphonates"):
        p = scm.propensity(drug, {k: pop[k].to_numpy() for k in scm.order})
        assert p.mean() == pytest.approx(pop[drug].mean(), abs=0.01)
        assert np.corrcoef(p, pop[drug])[0, 1] > 0.1
    with pytest.raises(ValueError):
        scm.propensity("hydration", {k: pop[k].to_numpy() for k in scm.order})


def test_reference_key_has_a_zero_total_for_sufficient_covariates():
    from attribution_validation.reference_key import reference_credits
    scm = TrueSCM(CFG)
    rng = np.random.default_rng(6)
    u = scm.noise(20000, rng)
    v = scm.simulate(u)
    pop = pd.DataFrame({**v, **u})
    cov = ["mission_era", "individual_factors", "urine_concentration", "hydration", "pre_flight_fitness"]
    r = reference_credits(scm, pop, u, "k_citrate", cov, np.ones(20000, dtype=bool), permutations=8, seed=1, jobs=2)
    assert abs(sum(r["credits"].values()) - (r["v_full"] - r["v_empty"])) < 1e-9
    assert r["v_empty"] == pytest.approx(-(r["crude"] - r["ate"]), abs=1e-9)         # v(empty) = -(crude - ATE)
    assert r["credits"]["individual_factors"] > 0 or r["credits"]["urine_concentration"] > 0


def test_secret_env_is_loaded_by_name_without_overriding(tmp_path, monkeypatch):
    f = tmp_path / ".env"
    f.write_text("DATA_DIR=./x\nTabPFNAPI=abc123\nHuggingFace=hf_demo\nANTHROPIC_API_KEY=should_not_be_exported\n")
    for var in ("TABPFN_TOKEN", "HF_TOKEN", "ANTHROPIC_API_KEY"):
        monkeypatch.delenv(var, raising=False)
    assert sorted(D.load_secret_env(f)) == ["HF_TOKEN", "TABPFN_TOKEN"]
    import os
    assert os.environ["TABPFN_TOKEN"] == "abc123" and os.environ["HF_TOKEN"] == "hf_demo"
    assert "ANTHROPIC_API_KEY" not in os.environ
    monkeypatch.setenv("TABPFN_TOKEN", "already")
    assert D.load_secret_env(f) == []                                      # an exported value is never replaced
    assert os.environ["TABPFN_TOKEN"] == "already"
    assert D.load_secret_env(tmp_path / "missing.env") == []
