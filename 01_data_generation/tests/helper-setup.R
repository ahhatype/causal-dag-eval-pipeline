dg_dir <- normalizePath(file.path(testthat::test_path(), ".."))
for (f in c("model_spec", "calibrate", "dag", "select", "ground_truth", "finalize", "validate_dag", "checks", "pipeline"))
  source(file.path(dg_dir, "R", paste0(f, ".R")))
cfg_dir <- file.path(dg_dir, "config")
cfg <- read_config(cfg_dir)
cal <- read_calibration(file.path(cfg_dir, "calibration.yaml"))
D <- build_dag(cfg, cal)

small_cfg <- function() {
  c2 <- cfg
  c2$params$sizes[c("full_set", "reference_subsample", "astronaut_set", "ground_truth")] <- list(3000, 300, 300, 3000)
  c2$params$replicates$count <- 2
  c2
}
