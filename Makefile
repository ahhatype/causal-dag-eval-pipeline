.PHONY: data recalibrate shape attribution convergence confounding test

data:            ## generate all replicate data sets + ground truth (uses frozen calibration)
	Rscript 01_data_generation/generate_all.R

recalibrate:     ## re-run the pilot calibration, then generate
	Rscript 01_data_generation/generate_all.R --recalibrate

shape:           ## structural validation: implication tests, power control and causal discovery
	Rscript 03_shape_validation/run_all.R

attribution:     ## standard and causal SHAP against simulated total effects
	cd 04_attribution_validation && uv run python run_all.py

convergence:     ## causal predictive SHAP convergence on phi: two independent halves at increasing order budgets
	cd 04_attribution_validation && uv run python run_ng_convergence.py

confounding:     ## ConfoundingSHAP: reference credits from the true model, then the TabPFN runs
	cd 04_attribution_validation && uv run python run_confounding_shap.py key && uv run python run_confounding_shap.py run

test:
	Rscript -e 'testthat::test_dir("01_data_generation/tests"); testthat::test_dir("03_shape_validation/tests")'
	cd 04_attribution_validation && uv run pytest -q
