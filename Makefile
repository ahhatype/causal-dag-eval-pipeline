.PHONY: data data-nonlinear recalibrate-nonlinear recalibrate structure attribution convergence confounding table1 credit conditions figures test

data:            ## generate all replicate data sets + ground truth (uses frozen calibration)
	Rscript 01_data_generation/generate_all.R

recalibrate:     ## re-run the pilot calibration, then generate
	Rscript 01_data_generation/generate_all.R --recalibrate

data-nonlinear:  ## supplement: the same DAG with three curved edges and one saturating edge, written to $(DATA_DIR)_nonlinear_symmetric_u
	Rscript 01_data_generation/generate_all.R --variant=nonlinear_symmetric_u

recalibrate-nonlinear:  ## re-run the pilot calibration for the nonlinear variant, then generate
	Rscript 01_data_generation/generate_all.R --recalibrate --variant=nonlinear_symmetric_u

structure:       ## structural validation: implication tests, power control and causal discovery
	Rscript 03_structure_validation/run_all.R

attribution:     ## standard and causal SHAP against simulated total effects
	cd 04_attribution_validation && uv run python run_all.py

convergence:     ## causal predictive SHAP convergence on phi: two independent halves at increasing order budgets
	cd 04_attribution_validation && uv run python run_ng_convergence.py

confounding:     ## ConfoundingSHAP: reference credits from the true model, then the TabPFN runs
	cd 04_attribution_validation && uv run python run_confounding_shap.py key && uv run python run_confounding_shap.py run

table1:          ## Table 1 data, ancestor displacement and the rendered table, from the stored results
	cd 04_attribution_validation && uv run python make_table1.py

credit:          ## share of attribution by feature role, from the stored importance
	cd 04_attribution_validation && uv run python make_credit_by_role.py

conditions:      ## sample-size and selection contrasts for both modules, from the stored results
	Rscript 03_structure_validation/run_conditions.R
	cd 04_attribution_validation && uv run python make_data_condition_contrasts.py

figures:         ## manuscript figures and LaTeX tables, from the stored results
	cd 05_figures && python3 fig2_difference_dag.py && python3 fig3_attribution.py && python3 table1_latex.py && python3 table1_latex.py --full && python3 supp_tables.py

test:
	Rscript -e 'testthat::test_dir("01_data_generation/tests"); testthat::test_dir("03_structure_validation/tests")'
	cd 04_attribution_validation && uv run pytest -q
	cd 05_figures && python3 -m pytest -q tests
