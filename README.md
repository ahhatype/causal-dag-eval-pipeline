# causal-dag-eval-pipeline

## Overview

<!-- TO BE WRITTEN BY THE AUTHOR: purpose of the project, study design, and how the pieces fit together. -->

_(Overview to be written.)_

---

## Setup

R is needed for data generation. Package versions (`simcausal`, `data.table`, `yaml`, `dagitty`, `testthat`) are pinned in `renv.lock`.

```bash
Rscript -e 'renv::restore()'   # install the pinned R packages
cp .env.example .env           # then add your API key; .env is git-ignored
make data                      # generate the data sets
make test                      # run the generation tests
```

`DATA_DIR` in `.env` (default `./data`) sets where generated files are written; it is git-ignored.

## Data generation

All code lives in [data_generation/](data_generation/). Three synthetic data sets are drawn from one prespecified DAG with [simcausal](https://cran.r-project.org/package=simcausal), so the true causal structure and total effects are known. The DAG is the NASA Human System Risk Board renal-stone DAG (SA-07566; 51 nodes, 75 edges) with three added nodes (**cumulative mission duration** as the exposure, **mission era** as a root, **pre-flight fitness** as a root), giving 54 nodes and 90 edges. The outcome is **nephrolithiasis** (binary). NASA's DAG specifies structure only: all equations, coefficients, era values and detection parameters here are study design parameters, not estimates of real-world effects.

| | Full-set | Reference subsample | Astronaut-set |
|---|---|---|---|
| n | 10,000 | 900 | 900 |
| Source | Simulated from the DAG | Random subsample of the full-set | Same source population, selected by sampling |
| Selection | None | None | Retained with probability expit(α − 1.0·individual factors + 1.0·pre-flight fitness), α set for ~30% retention |
| Era mix | 25% 1960s / 75% 2000s | as full-set | not era-dependent |
| Seed | 20261004 | 20261007 | 20261005 |
| Role | Ideal data, no selection | Separates sample-size from selection effects | Small, selected, astronaut-like cohort |

Ground-truth total effects use seed 20261006. Output files are the same set of 46 observed columns in every data set (7 latent nodes withheld; latent values go to separate `*_latent.csv` debug files).

### How it works

1. **Config** ([data_generation/config/](data_generation/config/)): the single source of truth.
   - `nodes.csv`, `edges.csv`: node roles/scales and every edge with its type and coefficient. Built from the DAG appendix (Supplementary Tables 1–2) by `tools/build_config_tables.py`.
   - `params.yaml`: seeds, sample sizes, era values, uptake probabilities, detection sensitivity/specificity, selection parameters, base rates. Entries marked PROVISIONAL are pending domain-expert confirmation.
   - `dag.txt`: the adapted DAGitty code used to validate the config.
2. **Validation** (`R/validate_dag.R`): the config must reproduce `dag.txt` exactly: 54 nodes, 90 edges, same edge set, acyclic, 16 outcome ancestors, 19 other observed features, 7 latent nodes.
3. **Calibration** (`R/calibrate.R`): a pilot run (n = 200,000) freezes, in topological order, each continuous node's standardization and residual SD, the logit intercept of each binary node (e.g. 10% nephrolithiasis prevalence), and the selection intercept α. Results are stored in `config/calibration.yaml` with a hash of the config files. Generation stops if the config has changed since calibration; re-run with `make recalibrate`.
4. **Model** (`R/model_spec.R`, `R/dag.R`): each node becomes a formula built from `edges.csv` and `params.yaml`.
   - Mission era is Bernoulli (0.75 for the 2000s). The 10 observed and 2 latent mission-level (context) nodes are fixed by era, so edges among them are kept in the graph but not separately parameterized. Cumulative mission duration is drawn from truncated era-specific distributions (1960s; 2000s Shuttle/ISS mixture).
   - Continuous children: linear combination of parents plus noise, in SD units of the child. Continuous parents enter standardized; binary parents (drugs, procedures) enter as 0/1, so their coefficient is the shift for users vs non-users. Binary children: logistic. Detection nodes use sensitivity/specificity; near-deterministic edges use a 5% flip.
   - Every node has its own uniform noise node, so interventions change a node's formula without shifting any random draws (**common random numbers**).
5. **Data sets** (`R/pipeline.R`): the full-set is simulated; the reference subsample is a random 900 of it; the astronaut-set simulates a larger source population, retains records by the selection rule and keeps the first 900 retained (`R/select.R`).
6. **Ground truth** (`R/ground_truth.R`): total effect of each of the 35 features on nephrolithiasis, as a risk difference in the outcome's expected probability under intervention (75th vs 25th percentile for continuous nodes, 1 vs 0 for binary), on 50,000 records, with a Monte Carlo standard error. For the astronaut target, effects are averaged over records retained at baseline, with selection held fixed.
7. **Checks** (`R/checks.R`): sizes, prevalence, era share, retention, selection-induced association, 35 features and no latent columns, non-ancestors have exactly zero effect. A failed check stops generation. The exception is the ancestor ranking: adjacent ancestors whose effects lie within 2 standard errors of each other are reported as a warning (near-tie). The current parameters give one: water intake and thiazides.

### Run

```bash
make data                                   # or: Rscript data_generation/generate_all.R
```

Takes about 2 minutes. Outputs in `data/`:

| File | Contents |
|---|---|
| `full_set.csv`, `reference_subsample.csv`, `astronaut_set.csv` | Observed columns (46 incl. outcome) |
| `*_latent.csv` | Latent nodes, for debugging only |
| `ground_truth_total_effects_full.csv`, `…_astronaut.csv` | Total effect, SE and rank per feature |
| `feature_sets.csv` | The 35-feature (all) and 16-feature (ancestor) sets |
| `data_summary.csv` | n, seed, era share, realized prevalence, retention |
| `provenance.txt` | Git commit, config hash, R and package versions |

### Tests

`make test` runs `data_generation/tests/`: DAG fidelity, topological order, seed reproducibility, common-random-number behavior, calibration targets and staleness, duration truncation, astronaut selection, ground truth, the checks themselves, an end-to-end run at reduced size, and (when `references/` is present locally) that the config tables rebuild exactly from the appendix.

## License

MIT, see [LICENSE](LICENSE).
