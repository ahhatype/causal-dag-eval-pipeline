#!/usr/bin/env Rscript
# Usage: Rscript data_generation/generate_all.R [--recalibrate]

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
dg_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else normalizePath(".")
repo_dir <- dirname(dg_dir)
for (f in c("model_spec", "calibrate", "dag", "select", "ground_truth", "finalize", "validate_dag", "checks", "pipeline"))
  source(file.path(dg_dir, "R", paste0(f, ".R")))

data_dir <- "./data"
env_file <- file.path(repo_dir, ".env")
if (file.exists(env_file)) {
  line <- grep("^\\s*DATA_DIR\\s*=", readLines(env_file, warn = FALSE), value = TRUE)
  if (length(line)) data_dir <- trimws(gsub("^['\"]|['\"]$", "", trimws(sub("^[^=]*=", "", line[1]))))
}
out_dir <- if (grepl("^/", data_dir)) data_dir else file.path(repo_dir, sub("^\\./", "", data_dir))

cfg_dir <- file.path(dg_dir, "config")
cfg <- read_config(cfg_dir)
message("[1/6] validating DAG against config/dag.txt")
validate_dag(cfg_dir, cfg)

cal_file <- file.path(cfg_dir, "calibration.yaml")
if ("--recalibrate" %in% args || !file.exists(cal_file)) {
  message("[2/6] calibrating (pilot n = ", cfg$params$sizes$calibration, ")")
  calibrate(cfg, cal_file)
} else {
  message("[2/6] using frozen calibration: ", cal_file)
}
cal <- read_calibration(cal_file)
check_calibration_current(cfg, cal)

out <- run_pipeline(cfg, cal, out_dir, repo_dir)
print(out$summary, row.names = FALSE)
print(out$checks, row.names = FALSE)
failed <- out$checks[!out$checks$pass, ]
for (i in which(failed$severity == "warning")) warning(failed$check[i], ": ", failed$detail[i], call. = FALSE, immediate. = TRUE)
if (any(failed$severity == "error"))
  stop("checks failed: ", paste(failed$check[failed$severity == "error"], collapse = ", "), call. = FALSE)
message("Wrote data sets to ", out_dir)
