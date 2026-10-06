#!/usr/bin/env Rscript
# Usage: Rscript 03_shape_validation/run_calibration.R [--replicates=N]
# Re-runs only the implication tests on every replicate (no discovery) to check their calibration.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
sv_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else normalizePath(".")
repo_dir <- dirname(sv_dir)
for (f in c("paths", "model_spec", "calibrate")) source(file.path(repo_dir, "01_data_generation", "R", paste0(f, ".R")))
for (f in c("analysis_graph", "implications", "dsep", "discovery", "compare", "power", "pipeline"))
  source(file.path(sv_dir, "R", paste0(f, ".R")))

params <- yaml::read_yaml(file.path(sv_dir, "config", "params.yaml"))
arg <- function(name) {
  a <- grep(sprintf("^--%s=", name), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(a)) as.integer(sub("^[^=]*=", "", a[1])) else NULL
}
data_dir <- env_dir(repo_dir, "DATA_DIR", "./02_data")
out_dir <- file.path(env_dir(repo_dir, "OUTPUT_DIR", "./outputs"), "03_shape_validation")
cfg_dir <- file.path(repo_dir, "01_data_generation", "config")
cfg <- read_config(cfg_dir)
cal <- read_calibration(file.path(cfg_dir, "calibration.yaml"))

summ <- run_calibration(cfg, cal, params, data_dir, out_dir, arg("replicates"))
print(summ)
message("Wrote implication calibration to ", out_dir)
