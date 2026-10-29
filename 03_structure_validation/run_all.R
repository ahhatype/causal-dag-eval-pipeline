#!/usr/bin/env Rscript
# Usage: Rscript 03_structure_validation/run_all.R [--replicates=N] [--subsamples=N] [--resume | --fresh] [--variant=<name>]
# --variant=<name> validates the data of config/variants/<name> (read from <DATA_DIR>_<name>, written to outputs/03_structure_validation_<name>).

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
if (!is.null(arg("subsamples"))) params$stability$subsamples <- arg("subsamples")
flags <- commandArgs(trailingOnly = TRUE)

variant <- sub("^--variant=", "", grep("^--variant=", flags, value = TRUE))
data_dir <- env_dir(repo_dir, "DATA_DIR", "./02_data")
out_dir <- file.path(env_dir(repo_dir, "OUTPUT_DIR", "./outputs"), "03_structure_validation")
cfg_dir <- file.path(repo_dir, "01_data_generation", "config")
if (length(variant)) {   # a generator variant: its own configuration and frozen calibration, its own data and output directories
  cfg_dir <- file.path(cfg_dir, "variants", variant)
  if (!dir.exists(cfg_dir)) stop("no such configuration: ", cfg_dir)
  data_dir <- file.path(dirname(data_dir), paste0(basename(data_dir), "_", variant))
  out_dir <- file.path(dirname(out_dir), paste0("03_structure_validation_", variant))
}
cfg <- read_config(cfg_dir)
cal <- read_calibration(file.path(cfg_dir, "calibration.yaml"))

out <- run_structure(cfg, cal, params, data_dir, out_dir, arg("replicates"), repo_dir,
                 resume = "--resume" %in% flags, fresh = "--fresh" %in% flags, sv_dir = sv_dir)
print(t(out$summary[, c("dataset", "variant", "replicates", "selection_violation_detection_rate",
                        "any_unexpected_rejection_rate", "edges_detectable_mean", "adjacency_precision_mean",
                        "adjacency_recall_mean", "orientation_agreement_mean")]))
message("Wrote structure validation outputs to ", out_dir)
