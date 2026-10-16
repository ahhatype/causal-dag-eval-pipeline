# Sample-size and selection contrasts of the shape-validation metrics. Kept outside R/ so that the code hash of the
# stored replicate units is unchanged.
DATA_CONDITION_METRICS <- c("edges_detectable", "adjacency_recall", "adjacency_precision", "orientation_agreement",
                            "raw_rejection_rate")
DATA_CONDITIONS <- list(sample_size_effect = c("reference_subsample", "full_set"),
                        selection_effect = c("astronaut_set", "reference_subsample"))

# Per metric and contrast: the first data set minus the second, paired by replicate. Replicates missing either
# value (for example no orientation to compare) drop out of that pair; n counts the pairs used.
data_condition_contrasts <- function(per_rep, variant, metrics = DATA_CONDITION_METRICS, z = 1.96) {
  x <- per_rep[per_rep$variant == variant, ]
  by_replicate <- function(m, ds) {
    s <- x[x$dataset == ds, ]
    setNames(s[[m]], as.character(s$replicate))
  }
  rows <- list()
  for (m in metrics) for (k in names(DATA_CONDITIONS)) {
    a <- by_replicate(m, DATA_CONDITIONS[[k]][1])
    d <- a - by_replicate(m, DATA_CONDITIONS[[k]][2])[names(a)]
    d <- d[!is.na(d)]
    n <- length(d)
    se <- if (n > 1) sd(d) / sqrt(n) else NA_real_
    rows[[length(rows) + 1]] <- data.frame(variant = variant, metric = m, contrast = k, mean = mean(d), mcse = se,
                                           ci_low = mean(d) - z * se, ci_high = mean(d) + z * se, n = n)
  }
  do.call(rbind, rows)
}
