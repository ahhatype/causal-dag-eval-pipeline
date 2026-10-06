.PHONY: data recalibrate shape test

data:            ## generate all three data sets + ground truth (uses frozen calibration)
	Rscript 01_data_generation/generate_all.R

recalibrate:     ## re-run the pilot calibration, then generate
	Rscript 01_data_generation/generate_all.R --recalibrate

shape:           ## structural validation: implication tests and causal discovery on the generated data
	Rscript 03_shape_validation/run_all.R

test:
	Rscript -e 'testthat::test_dir("01_data_generation/tests"); testthat::test_dir("03_shape_validation/tests")'
