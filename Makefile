.PHONY: data recalibrate test

data:
	Rscript data_generation/generate_all.R

recalibrate:
	Rscript data_generation/generate_all.R --recalibrate

test:
	Rscript -e 'testthat::test_dir("data_generation/tests")'
