retention_mask <- function(df, cfg, cal) {
  s <- cfg$params$selection
  df$U_selection < plogis(cal$selection_alpha +
                            s$b_individual_factors * df$individual_factors +
                            s$b_fitness * df$pre_flight_fitness)
}

# The astronaut-set is the first n_target retained records of the same replicate's full-set, so the two are paired.
astronaut_from_full <- function(full, cfg, cal, n_target) {
  keep <- retention_mask(full, cfg, cal)
  if (sum(keep) < n_target) stop(sprintf("only %d of %d full-set records retained; need %d", sum(keep), nrow(full), n_target))
  list(data = full[keep][seq_len(n_target)], retained = sum(keep), source_retention = mean(keep))
}
