# Results, arm A: the 54-node generator with asymmetric curved mechanisms (13 October 2026)

Arm A = PREREGISTRATION.md as written: inverted U (Nutrients to Bone Remodeling), U (Microbiome to Urine Chemistry),
U (Hydration to Urine Concentration), all centred at 0.2 on the standardized parent, plus a saturating Urine
Chemistry to Mineralized Renal Material. Data: 50 replicates of the three data sets, regenerated and recalibrated
on branch `nonmonotone-variant` (commit 60d1ca3). Pipeline shape validation: `causal-dag-eval-pipeline/outputs_nonmono/`.
My hybrid and audit (20 replicates, full-set): `hybrid_summary.csv`, `hybrid_by_rep.csv`, `gaps.csv`, `audit_summary.csv`,
`implications_all.csv`.

## 1. What the pipeline's own discovery (MGM-PC-Stable) does with the curved edges

Adjacency recovery rate over 50 replicates, variant / linear generator:

| Edge | Full set (10,000) | Reference subsample (900) | Astronaut set (900) |
|---|---|---|---|
| Nutrients to Bone Remodeling (inverted U) | 1.00 / 1.00 | 0.46 / 1.00 | 0.30 / 1.00 |
| Microbiome to Urine Chemistry (U) | 1.00 / 1.00 | 0.42 / 0.98 | 0.42 / 0.96 |
| Hydration to Urine Concentration (U) | 1.00 / 1.00 | 0.44 / 1.00 | 0.48 / 1.00 |
| Urine Chemistry to Mineralized Renal Material (saturating control) | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 |

Orientation of the curved edges in the learned CPDAG on the full set: Hydration to Urine Concentration is left
undirected in every replicate (0.00 oriented, linear generator 0.36); the other two are oriented as often as before.
Overall full-set recall 0.78 and precision 0.995, the same as the linear generator; implication tests: 1.1
rejections per replicate on the full set (linear generator 0.04), all on the four flagged non-additive children.

A finding that predates the variant and is worth putting to Andy: the two countermeasure edges K+ Citrate to Urine
Chemistry and Thiazides to Urine Chemistry are never found (0.00 adjacency in every data set) on either generator.
MGM-PC does not recover the drug effects the paper is about.

## 2. Hybrid and audit on 20 full-set replicates

| Method | Edges | Precision | Recall | F1 | U edges found (of 3) | U edges oriented (of 3) | False edges |
|---|---|---|---|---|---|---|---|
| PC-stable alone (Fisher-z) | 55.2 | 0.98 | 0.85 | 0.91 | 3.0 | 1.95 | 1.4 |
| Hybrid, gap rule, no veto | 59.6 | 0.90 | 0.85 | 0.88 | 3.0 | 3.0 | 5.8 |
| Hybrid, chain veto | 57.9 | 0.93 | 0.85 | 0.89 | 3.0 | 3.0 | 4.1 |

The gap rule orients all three curved edges in every replicate (gaps 0.19 to 0.41, never near the band) but adds
2.75 false edges per replicate that the chain veto does not remove: three common-cause pairs of Cumulative Mission
Duration with Hydration, Resistive Exercise and Water Intake (20 of 20 replicates each), and the Hydration to Urine
Chemistry chain (20 of 20, as on nl_v2). Cumulative Mission Duration is the exposure, drawn from an era-dependent
truncated mixture: skewed and bimodal. The R2 asymmetry rule reads a skewed cause as "oriented", which folder 15's
pre-registration listed as its known weakness. A PC-adjacency veto would remove all three (none are PC-adjacent).

Audit (Shipley basis set, about 567 implications):

| Graph | Test | Rejections (mean) | Share | Deletion implications rejected (of 3) |
|---|---|---|---|---|
| True analysis graph | linear LRT / basis GCM | 1.1 / 1.9 | 0.002 / 0.003 | |
| Minus the 3 U edges | linear LRT / basis GCM | 5.2 / 8.7 | 0.009 / 0.015 | **3.0 / 3.0** |
| PC skeleton, ordered | linear LRT / basis GCM | 12.6 / 15.1 | 0.022 / 0.026 | |

## Reading
- With the U centred at 0.2 the curved edges keep marginal correlations of 0.12 to 0.23. At 10,000 rows that is
  enough for Fisher-z and for the linear audit, so on the paper's full set the "linearity hides" claim does not
  hold for the skeleton; it holds for orientation (one U edge never oriented) and, at the paper's astronaut scale of
  900 rows, for the skeleton too (curved edges found in under half the replicates against essentially always for
  their linear counterparts).
- The asymmetric centre is the physiologically honest choice (optima are rarely at the population mean), and it
  is also what makes the effect sample-size dependent. Arm B (symmetric U, centre 0, zero linear projection) shows
  the other end: see RESULTS_armB.md.
- The hybrid's value on this generator is orientation, not recall. Its false edges come from a skewed exposure
  distribution and need the skeleton veto, which the pre-registered rule did not include.
