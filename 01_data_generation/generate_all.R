#!/usr/bin/env Rscript
# Usage: Rscript 01_data_generation/generate_all.R [--recalibrate]

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
dg_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else normalizePath(".")
repo_dir <- dirname(dg_dir)
for (f in c("paths", "model_spec", "calibrate", "dag", "select", "ground_truth", "finalize", "validate_dag", "checks", "pipeline"))
  source(file.path(dg_dir, "R", paste0(f, ".R")))

out_dir <- env_dir(repo_dir, "DATA_DIR", "./02_data")

cfg_dir <- file.path(dg_dir, "config")
cfg <- read_config(cfg_dir)
message("[1/4] validating DAG against config/dag.txt")
validate_dag(cfg_dir, cfg)

cal_file <- file.path(cfg_dir, "calibration.yaml")
if ("--recalibrate" %in% args || !file.exists(cal_file)) {
  message("[2/4] calibrating (pilot n = ", cfg$params$sizes$calibration, ")")
  calibrate(cfg, cal_file)
} else {
  message("[2/4] using frozen calibration: ", cal_file)
}
cal <- read_calibration(cal_file)
check_calibration_current(cfg, cal)

out <- run_pipeline(cfg, cal, out_dir, repo_dir)
agg <- aggregate(cbind(n, share_2000s, prevalence, if_fitness_cor, retention) ~ dataset, out$summary,
                 FUN = mean, na.action = na.pass)
print(agg, row.names = FALSE)
print(out$checks, row.names = FALSE)
failed <- out$checks[!out$checks$pass, ]
for (i in which(failed$severity == "warning")) warning(failed$check[i], ": ", failed$detail[i], call. = FALSE, immediate. = TRUE)
if (any(failed$severity == "error"))
  stop("checks failed: ", paste(failed$check[failed$severity == "error"], collapse = ", "), call. = FALSE)
message("Wrote data sets to ", out_dir)
