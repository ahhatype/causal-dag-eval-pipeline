# Results, arm B: the 54-node generator with symmetric U-shapes (14 October 2026)

Arm B was added after arm A showed that U-shapes centred at 0.2 keep enough linear projection for Fisher-z to find
them at 10,000 rows. Same three curved edges and the saturating control, with the U centred at the parent's mean
(zero linear projection; commit 74a7e72). 50 replicates regenerated and recalibrated; pipeline shape validation in
`outputs_nonmono_sym/` (one rCausalMGM segmentation fault at unit 4, resumed without further incident). My hybrid
and audit on 20 full-set replicates: `armB/`.

## 1. The pipeline's discovery (MGM-PC-Stable, 50 replicates)

| Edge | Full set (10,000) | Reference subsample (900) | Astronaut set (900) |
|---|---|---|---|
| Nutrients to Bone Remodeling (symmetric inverted U) | 0.00 | 0.00 | 0.00 |
| Microbiome to Urine Chemistry (symmetric U) | 0.00 | 0.00 | 0.00 |
| Hydration to Urine Concentration (symmetric U) | 0.00 | 0.00 | 0.00 |
| Urine Chemistry to Mineralized Renal Material (saturating control) | 1.00 | 1.00 | 1.00 |

Full-set adjacency recall 0.74 (0.78 on the linear generator and on arm A), precision 0.99. The pipeline's linear
implication tests reject 1.1 implications per replicate on the full set, all on the flagged non-additive children,
and their calibration on the implications expected to hold stays uniform (KS rejected in 6% of replicates).

## 2. Hybrid and audit, 20 full-set replicates

| Method | Edges | Precision | Recall | F1 | U edges found (of 3) | U edges oriented (of 3) | False edges |
|---|---|---|---|---|---|---|---|
| PC-stable alone (Fisher-z) | 53.0 | 0.96 | 0.81 | 0.88 | 0.1 | 0.1 | 2.0 |
| Hybrid, gap rule, no veto | 60.0 | 0.90 | 0.86 | 0.88 | 3.0 | 3.0 | 6.1 |
| Hybrid, chain veto | 59.2 | 0.91 | 0.86 | 0.88 | 3.0 | 3.0 | 5.2 |

Gap-rule margins on the curved edges 0.19 to 0.44, never near the 0.10 band. False edges as in arm A: the three
common-cause pairs of the skewed exposure (Cumulative Mission Duration with Hydration, Resistive Exercise, Water
Intake; 20 of 20 each) and the Hydration to Urine Chemistry chain (20 of 20). None is PC-adjacent.

Audit (Shipley basis set):

| Graph | Test | Rejections (mean) | Share | Deletion implications rejected (of 3) |
|---|---|---|---|---|
| True analysis graph | linear LRT / basis GCM | 1.1 / 1.9 | 0.002 / 0.003 | |
| Minus the 3 U edges | linear LRT / basis GCM | 1.1 / 9.0 | 0.002 / 0.016 | **0.0 / 3.0** |
| PC skeleton, ordered | linear LRT / basis GCM | 12.5 / 20.6 | 0.022 / 0.036 | 0.15 / 2.9 |

## Reading
- This is the "linearity hides" demonstration inside the paper's own experiment. With symmetric curvature the
  pipeline's discovery step never finds the three curved edges in 150 data sets, its linear audit certifies the
  graph that lacks them (0 of 3 deletion implications rejected), and both the flexible two-way fit test (hybrid)
  and the basis-expanded GCM audit recover them in every replicate. The saturating edge, monotone, is found every
  time by every method, which isolates non-monotonicity as the cause.
- Together the two arms give the paper a clean statement with its boundary: a linear pipeline misses a curved edge
  when the curve's linear projection is small relative to the sample; at 900 rows that covers physiologically
  asymmetric U-shapes (arm A, found in 30 to 48% of replicates); at 10,000 rows it covers symmetric ones (arm B,
  never found) but not asymmetric ones (arm A, always found).
- Two things to carry into the method: the gap rule needs the PC-adjacency veto for skewed causes (the exposure
  here is an era-dependent mixture), and the audit needs a dependence-sensitive test, not a residual-correlation one.
- Not run: `make attribution` on either arm (about 10 hours per arm on this machine); the attribution results for
  the curved generator are the natural next step if Andy adopts the variant.
