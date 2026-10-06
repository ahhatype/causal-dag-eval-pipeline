retention_mask <- function(df, cfg, cal) {
  s <- cfg$params$selection
  df$U_selection < plogis(cal$selection_alpha +
                            s$b_individual_factors * df$individual_factors +
                            s$b_fitness * df$pre_flight_fitness)
}

draw_astronaut_set <- function(D, cfg, cal, n_target, seed) {
  n_src <- ceiling(n_target / cfg$params$selection$target_retention * 1.2)
  repeat {
    src <- sim_source(D, n_src, seed)
    keep <- retention_mask(src, cfg, cal)
    if (sum(keep) >= n_target) break
    n_src <- ceiling(n_src * 1.5)
  }
  list(data = src[keep][seq_len(n_target)], source_n = n_src, source_retention = mean(keep))
}
