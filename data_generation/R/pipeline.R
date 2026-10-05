run_pipeline <- function(cfg, cal, out_dir, repo_dir = NULL) {
  p <- cfg$params
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  D <- build_dag(cfg, cal)

  message("[3/6] full-set")
  full <- sim_source(D, p$sizes$full_set, p$seeds$full_set)
  message("[4/6] reference subsample")
  set.seed(p$seeds$reference_subsample)
  reference <- full[sort(sample(nrow(full), p$sizes$reference_subsample))]
  message("[5/6] astronaut-set")
  astro <- draw_astronaut_set(D, cfg, cal, p$sizes$astronaut_set, p$seeds$astronaut_source)
  message("[6/6] ground-truth total effects")
  gt <- ground_truth(D, cfg, cal)

  write_dataset(full, "full_set", cfg, out_dir)
  write_dataset(reference, "reference_subsample", cfg, out_dir)
  write_dataset(astro$data, "astronaut_set", cfg, out_dir)
  write_feature_sets(cfg, out_dir)
  gt$name <- cfg$nodes$name[match(gt$feature, cfg$nodes$id)]
  for (t in c("full", "astronaut")) {
    cols <- c("feature", "name", "ancestor_set", paste0(c("low_", "high_", "effect_", "se_", "rank_"), t))
    data.table::fwrite(gt[, cols], file.path(out_dir, sprintf("ground_truth_total_effects_%s.csv", t)))
  }

  summ <- data.frame(
    dataset = c("full_set", "reference_subsample", "astronaut_set"),
    n = c(nrow(full), nrow(reference), nrow(astro$data)),
    seed = c(p$seeds$full_set, p$seeds$reference_subsample, p$seeds$astronaut_source),
    share_2000s = c(mean(full$mission_era), mean(reference$mission_era), mean(astro$data$mission_era)),
    prevalence = c(mean(full$nephrolithiasis), mean(reference$nephrolithiasis), mean(astro$data$nephrolithiasis)),
    source_n = c(nrow(full), nrow(full), astro$source_n),
    retention = c(NA, NA, astro$source_retention))
  data.table::fwrite(summ, file.path(out_dir, "data_summary.csv"))
  write_provenance(cfg, cal, out_dir, repo_dir)

  res <- list(full = full, reference = reference, astronaut = astro, ground_truth = gt)
  list(results = res, summary = summ, checks = run_checks(res, cfg))
}

write_provenance <- function(cfg, cal, out_dir, repo_dir = NULL) {
  pkgs <- c("simcausal", "data.table", "yaml", "dagitty")
  commit <- if (!is.null(repo_dir)) {
    h <- suppressWarnings(tryCatch(system2("git", c("-C", repo_dir, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
                                   error = function(e) NA_character_))
    dirty <- suppressWarnings(tryCatch(length(system2("git", c("-C", repo_dir, "status", "--porcelain"), stdout = TRUE)) > 0,
                                       error = function(e) NA))
    sprintf("%s%s", h[1], if (isTRUE(dirty)) " (uncommitted changes)" else "")
  } else NA_character_
  writeLines(c(
    paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("git commit:", commit),
    paste("config hash:", cal$config_hash),
    paste("R:", R.version.string),
    sprintf("%s: %s", pkgs, vapply(pkgs, function(x) as.character(utils::packageVersion(x)), ""))
  ), file.path(out_dir, "provenance.txt"))
}
