run_checks <- function(res, cfg) {
  p <- cfg$params
  obs <- observed_ids(cfg)
  feats <- feature_ids(cfg, "all")
  gt <- res$ground_truth
  non_anc <- gt$feature[!gt$ancestor_set]
  anc_eff <- gt[gt$ancestor_set, ]
  ast <- res$astronaut
  ck <- function(name, ok, detail, severity = "error") {
    data.frame(check = name, pass = isTRUE(ok), severity = severity, detail = detail)
  }
  rbind(
    ck("full_set_dims", nrow(res$full) == p$sizes$full_set, sprintf("n=%d", nrow(res$full))),
    ck("reference_dims", nrow(res$reference) == p$sizes$reference_subsample, sprintf("n=%d", nrow(res$reference))),
    ck("astronaut_dims", nrow(ast$data) == p$sizes$astronaut_set, sprintf("n=%d", nrow(ast$data))),
    ck("full_prevalence", abs(mean(res$full$nephrolithiasis) - p$base_rates$nephrolithiasis) < 0.01,
       sprintf("%.4f", mean(res$full$nephrolithiasis))),
    ck("era_share", abs(mean(res$full$mission_era) - p$era$share_2000s) < 0.015,
       sprintf("%.4f", mean(res$full$mission_era))),
    ck("astronaut_retention", abs(ast$source_retention - p$selection$target_retention) < 0.03,
       sprintf("%.4f of %d source records", ast$source_retention, ast$source_n)),
    ck("astronaut_prevalence_lower", mean(ast$data$nephrolithiasis) < mean(res$full$nephrolithiasis),
       sprintf("%.4f vs %.4f", mean(ast$data$nephrolithiasis), mean(res$full$nephrolithiasis))),
    ck("selection_induces_if_fitness_association",
       cor(ast$data$individual_factors, ast$data$pre_flight_fitness) > 0.05,
       sprintf("cor=%.3f", cor(ast$data$individual_factors, ast$data$pre_flight_fitness))),
    ck("feature_count_35", length(feats) == 35 && all(feats %in% obs), sprintf("%d", length(feats))),
    ck("no_latent_columns", !any(cfg$nodes$id[!cfg$nodes$observed] %in% obs), "latent nodes excluded"),
    ck("non_ancestors_zero_effect", all(abs(gt$effect_full[!gt$ancestor_set]) < 1e-12) &&
         all(abs(gt$effect_astronaut[!gt$ancestor_set]) < 1e-12),
       sprintf("%d non-ancestors", length(non_anc))),
    ck("ancestors_nonzero_effect", all(abs(anc_eff$effect_full) > 0), sprintf("%d ancestors", nrow(anc_eff))),
    {
      u <- c(unresolved_pairs(gt, "full"), unresolved_pairs(gt, "astronaut"))
      ck("ancestor_ranks_resolved", !length(u),
         if (length(u)) paste("within 2 SE:", paste(unique(u), collapse = "; ")) else "all adjacent gaps > 2 SE",
         severity = "warning")
    }
  )
}
