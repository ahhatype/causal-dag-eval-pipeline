.PHONY: data recalibrate shape attribution test

data:            ## generate all replicate data sets + ground truth (uses frozen calibration)
	Rscript 01_data_generation/generate_all.R

recalibrate:     ## re-run the pilot calibration, then generate
	Rscript 01_data_generation/generate_all.R --recalibrate

shape:           ## structural validation: implication tests, power control and causal discovery
	Rscript 03_shape_validation/run_all.R

attribution:     ## standard and causal SHAP against simulated total effects
	cd 04_attribution_validation && uv run python run_all.py

test:
	Rscript -e 'testthat::test_dir("01_data_generation/tests"); testthat::test_dir("03_shape_validation/tests")'
	cd 04_attribution_validation && uv run pytest -q
