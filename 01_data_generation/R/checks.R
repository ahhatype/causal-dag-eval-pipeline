run_checks <- function(summ, truth, cfg) {
  p <- cfg$params
  obs <- observed_ids(cfg)
  feats <- feature_ids(cfg, "all")
  pu <- truth[truth$truth_type == "per_unit", ]
  ck <- function(name, ok, detail, severity = "error") {
    data.frame(check = name, pass = isTRUE(ok), severity = severity, detail = detail)
  }
  by <- function(ds) summ[summ$dataset == ds, ]
  full <- by("full_set"); ref <- by("reference_subsample"); ast <- by("astronaut_set")
  target_n <- c(full_set = p$sizes$full_set, reference_subsample = p$sizes$reference_subsample,
                astronaut_set = p$sizes$astronaut_set)
  rng <- function(x) sprintf("mean %.4f, range %.4f-%.4f", mean(x), min(x), max(x))
  rbind(
    ck("replicate_count", nrow(full) == p$replicates$count && nrow(ref) == nrow(full) && nrow(ast) == nrow(full),
       sprintf("%d replicates", nrow(full))),
    ck("dataset_sizes", all(summ$n == target_n[summ$dataset]), "every replicate of every data set"),
    ck("full_prevalence", abs(mean(full$prevalence) - p$base_rates$nephrolithiasis) < 0.005 &&
         all(abs(full$prevalence - p$base_rates$nephrolithiasis) < 0.02), rng(full$prevalence)),
    ck("era_share", all(abs(full$share_2000s - p$era$share_2000s) < 0.02), rng(full$share_2000s)),
    ck("astronaut_retention", all(abs(ast$retention - p$selection$target_retention) < 0.035), rng(ast$retention)),
    ck("astronaut_prevalence_lower", mean(ast$prevalence) < mean(full$prevalence),
       sprintf("%.4f vs %.4f", mean(ast$prevalence), mean(full$prevalence))),
    ck("selection_induces_if_fitness_association", mean(ast$if_fitness_cor) > 0.05 && abs(mean(full$if_fitness_cor)) < 0.02,
       sprintf("astronaut %s; full %s", rng(ast$if_fitness_cor), rng(full$if_fitness_cor))),
    ck("feature_count_35", length(feats) == 35 && all(feats %in% obs), sprintf("%d", length(feats))),
    ck("no_latent_columns", !any(cfg$nodes$id[!cfg$nodes$observed] %in% obs), "latent nodes excluded"),
    ck("non_ancestors_zero_effect", all(abs(truth$value[!truth$ancestor_set]) < 1e-12),
       sprintf("%d non-ancestors", length(unique(truth$feature[!truth$ancestor_set])))),
    ck("ancestors_nonzero_effect", all(abs(truth$value[truth$ancestor_set]) > 0),
       sprintf("%d ancestors", length(unique(truth$feature[truth$ancestor_set])))),
    {
      u <- c(unresolved_pairs(truth, "source"), unresolved_pairs(truth, "selected"))
      ck("ancestor_ranks_resolved", !length(u),
         if (length(u)) paste("within 2 SE:", paste(unique(u), collapse = "; ")) else "all adjacent gaps > 2 SE",
         severity = "warning")
    }
  )
}
