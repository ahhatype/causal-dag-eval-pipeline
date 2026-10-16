replicate_seeds <- function(p, r) {
  shift <- p$replicates$seed_stride * (r - 1)
  list(full_set = p$seeds$full_set + shift, reference_subsample = p$seeds$reference_subsample + shift)
}

replicate_dir <- function(out_dir, r) file.path(out_dir, "replicates", sprintf("r%03d", r))

run_pipeline <- function(cfg, cal, out_dir, repo_dir = NULL) {
  p <- cfg$params
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  D <- build_dag(cfg, cal)

  message("[3/4] ground truth (n = ", p$sizes$ground_truth, ")")
  t0 <- Sys.time()
  gt <- ground_truth(cfg, cal)
  write_feature_sets(cfg, out_dir)
  write_truth(gt, cfg, out_dir)
  write_scm_check(D, cfg, out_dir)
  message(sprintf("      %d of %d truth records retained by selection; %.1f min",
                  gt$retained, gt$n, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  message(sprintf("[4/4] %d replicates of each data set", p$replicates$count))
  rows <- lapply(seq_len(p$replicates$count), function(r) {
    s <- replicate_seeds(p, r)
    full <- sim_source(D, p$sizes$full_set, s$full_set)
    set.seed(s$reference_subsample)
    reference <- full[sort(sample(nrow(full), p$sizes$reference_subsample))]
    astro <- astronaut_from_full(full, cfg, cal, p$sizes$astronaut_set)
    rd <- replicate_dir(out_dir, r)
    dir.create(rd, showWarnings = FALSE, recursive = TRUE)
    write_dataset(full, "full_set", cfg, rd)
    write_dataset(reference, "reference_subsample", cfg, rd)
    write_dataset(astro$data, "astronaut_set", cfg, rd)
    sets <- list(full_set = full, reference_subsample = reference, astronaut_set = astro$data)
    data.frame(
      replicate = r, dataset = names(sets),
      n = vapply(sets, nrow, integer(1)),
      seed = c(s$full_set, s$reference_subsample, s$full_set),
      share_2000s = vapply(sets, function(d) mean(d$mission_era), numeric(1)),
      prevalence = vapply(sets, function(d) mean(d$nephrolithiasis), numeric(1)),
      if_fitness_cor = vapply(sets, function(d) cor(d$individual_factors, d$pre_flight_fitness), numeric(1)),
      source_n = c(nrow(full), nrow(full), nrow(full)),
      retained = c(NA, NA, astro$retained),
      retention = c(NA, NA, astro$source_retention))
  })
  summ <- do.call(rbind, rows)
  rownames(summ) <- NULL
  data.table::fwrite(summ, file.path(out_dir, "data_summary.csv"))
  write_provenance(cfg, cal, out_dir, repo_dir)

  list(summary = summ, ground_truth = gt, checks = run_checks(summ, gt$values, cfg))
}

write_truth <- function(gt, cfg, out_dir) {
  v <- gt$values
  v$name <- cfg$nodes$name[match(v$feature, cfg$nodes$id)]
  cols <- c("feature", "name", "ancestor_set", "truth_type", "population", "value", "mcse", "low", "high")
  data.table::fwrite(v[, cols], file.path(out_dir, "truth_values.csv"))
  data.table::fwrite(gt$binary_checks, file.path(out_dir, "truth_binary_checks.csv"))
}

# A few simcausal records with their noise draws, so other modules can check their own copy of the equations.
write_scm_check <- function(D, cfg, out_dir, n = 500, seed = 1) {
  s <- sim_source(D, n, seed)
  data.table::fwrite(s[, c(cfg$order, paste0("U_", noise_ids(cfg))), with = FALSE], file.path(out_dir, "scm_check.csv"))
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
    paste("replicates:", cfg$params$replicates$count),
    paste("R:", R.version.string),
    sprintf("%s: %s", pkgs, vapply(pkgs, function(x) as.character(utils::packageVersion(x)), ""))
  ), file.path(out_dir, "provenance.txt"))
}
