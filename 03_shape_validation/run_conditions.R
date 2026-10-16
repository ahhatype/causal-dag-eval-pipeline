#!/usr/bin/env Rscript
# Usage: Rscript 03_shape_validation/run_conditions.R
# Sample-size and selection contrasts from the stored per_replicate.csv (nothing is rerun): data_condition_contrasts.csv.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
sv_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else normalizePath(".")
repo_dir <- dirname(sv_dir)
source(file.path(repo_dir, "01_data_generation", "R", "paths.R"))
source(file.path(sv_dir, "data_conditions.R"))

params <- yaml::read_yaml(file.path(sv_dir, "config", "params.yaml"))
out_dir <- file.path(env_dir(repo_dir, "OUTPUT_DIR", "./outputs"), "03_shape_validation")
per_rep <- data.table::fread(file.path(out_dir, "per_replicate.csv"), data.table = FALSE)
res <- data_condition_contrasts(per_rep, params$variants[[1]])
data.table::fwrite(res, file.path(out_dir, "data_condition_contrasts.csv"))
print(res, digits = 3)
message("Wrote data_condition_contrasts.csv to ", out_dir)
