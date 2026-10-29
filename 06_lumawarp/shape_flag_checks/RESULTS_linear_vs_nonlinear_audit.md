# Results: a linear audit passes a wrong graph; which nonlinear audit rejects it (13 October 2026)

Bars fixed in PREREGISTRATION.md. nl_v2, 8 seeds, n = 10,000; Shipley basis set (about 1,200 implications per
graph, X independent of Y given the parents of both), BH at 0.05. Tables: `implications_all.csv`,
`summary_by_seed.csv`, `summary.csv`; post hoc test in `implications_basis.csv`, `summary_basis.csv`.

## Pre-registered tests

| Graph audited | Test | Implications rejected (mean) | Share rejected | Graphs rejected (of 8) | Deletion implications rejected (of 4) |
|---|---|---|---|---|---|
| True graph | Linear/logistic LRT | 0.5 | 0.000 | 4 | |
| True graph | GCM, residual correlation | 0.5 | 0.000 | 4 | |
| True minus the 4 non-monotone edges | Linear/logistic LRT | 0.9 | 0.001 | 5 | **0.0** |
| True minus the 4 non-monotone edges | GCM, residual correlation | 1.6 | 0.001 | 7 | **0.0** |
| PC-stable skeleton, oriented by the ordering | Linear/logistic LRT | 0.1 | 0.000 | 1 | |
| PC-stable skeleton, oriented by the ordering | GCM, residual correlation | 0.9 | 0.001 | 3 | |

- **A1 (the linear audit is blind): PASS.** 0 of 4 deletion implications rejected in every seed.
- **A2 (the GCM sees them): FAIL.** Also 0 of 4. The GCM as specified tests the correlation of the two residuals,
  and a U-shaped dependence between residuals has zero correlation; conditioning on the parent's own parents
  reduces the parent to its innovation, which is exactly where the linear projection vanishes. The sine edge is
  hidden the same way (its phase varies with the conditioning set, so the within-set linear component averages out).
- **A3 (no false alarms): PASS.** Both tests reject 0.0% of the true graph's implications on average (half a
  rejection per 1,200).

## Post hoc test, labelled as such: GCM with a vector-valued Y
Regress X on Z and a cubic-spline basis of Y on Z (flexible regressions), test every residual product, take the
largest statistic with a Bonferroni correction. This is Shah and Peters' measure with Y expanded to a basis, so it
is sensitive to dependence, not only correlation.

| Graph audited | Implications rejected (mean) | Share rejected | Graphs rejected (of 8) | Deletion implications rejected (of 4) |
|---|---|---|---|---|
| True graph | 4.6 | 0.004 | 7 | |
| True minus the 4 non-monotone edges | 16.3 | 0.013 | 8 | **4.0** |
| PC-stable skeleton | 14.0 | 0.012 | 8 | (31 rejections on the 4 non-monotone pairs over 8 seeds; 81 elsewhere) |

All four deleted edges are flagged in all eight seeds; the false-rejection rate on the true graph is 0.4% of
implications (about 4.6 of 1,200), within the A3 bar but above the correlation tests' 0.0%.

## Reading
- The figure the paper's second pitch needs is here, with one correction to the rework note: swapping the linear
  test for "the GCM" is not enough. A residual-correlation test inherits the linear test's blindness to even
  dependence. The audit must test dependence between residuals (a basis-expanded GCM as here, or a kernel test
  such as HSIC or KCI), at the cost of a few more false rejections (0.4% of implications on the true graph).
- Stated for the manuscript: "A linear implication audit certified a graph with its non-monotone edges removed in
  every replicate; so did a residual-correlation GCM. A GCM applied to a spline basis of the tested variable
  rejected that graph in every replicate, flagging all four missing edges, at a false-rejection rate of 0.4%."
- On the discovery output itself (PC skeleton), the basis test rejects in every seed and about a quarter of its
  rejections sit on the four non-monotone pairs; the rest point at PC's other errors (its recall is 0.87). An audit
  that only counts rejections cannot tell the two apart; that is the job of the discrepancy trigger (folders 27, 29),
  which names the edge.
- The test is cheap (two least-squares fits per implication) and civilian; nothing here uses the engine.
