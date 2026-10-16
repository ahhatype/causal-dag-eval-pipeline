source(file.path(sv_dir, "data_conditions.R"))

toy_replicates <- function(reps = 10, variant = "MGM-PC-Stable") {
  noise <- 0.3 * (-1)^seq_len(reps) + 0.01 * seq_len(reps)    # shared across data sets: cancels in the pairing
  at <- c(full_set = 0, reference_subsample = 0.1, astronaut_set = 0.1 - 0.25)
  do.call(rbind, lapply(names(at), function(ds)
    data.frame(replicate = seq_len(reps), dataset = ds, variant = variant, edges_detectable = 50 + noise + at[[ds]],
               adjacency_recall = 0.6 + noise + at[[ds]], adjacency_precision = 0.9 + noise,
               orientation_agreement = 0.5 + noise + at[[ds]], raw_rejection_rate = 0.05 + noise - at[[ds]])))
}

test_that("contrasts recover known paired differences and ignore the other variants", {
  x <- toy_replicates()
  x <- rbind(x, transform(toy_replicates(variant = "MGM-CPC-Stable"), adjacency_recall = 0))
  res <- data_condition_contrasts(x, "MGM-PC-Stable")
  expect_equal(nrow(res), length(DATA_CONDITION_METRICS) * length(DATA_CONDITIONS))
  get <- function(m, k) res[res$metric == m & res$contrast == k, ]
  expect_equal(get("adjacency_recall", "sample_size_effect")$mean, 0.1)
  expect_equal(get("adjacency_recall", "selection_effect")$mean, -0.25)
  expect_equal(get("raw_rejection_rate", "selection_effect")$mean, 0.25)
  expect_equal(get("adjacency_precision", "sample_size_effect")$mean, 0)
  expect_equal(get("edges_detectable", "sample_size_effect")$mcse, 0, tolerance = 1e-12)
  expect_true(all(res$n == 10))
})

test_that("intervals are mean +/- 1.96 MCSE and pairs missing a value drop out", {
  x <- toy_replicates()
  x$orientation_agreement <- x$orientation_agreement + 0.02 * ifelse(x$replicate %% 2 == 0, 1, -1) *
    ifelse(x$dataset == "reference_subsample", 1, 0)
  x$orientation_agreement[x$dataset == "full_set" & x$replicate %in% c(2, 5)] <- NA
  r <- data_condition_contrasts(x, "MGM-PC-Stable")
  r <- r[r$metric == "orientation_agreement" & r$contrast == "sample_size_effect", ]
  expect_equal(r$n, 8)
  expect_equal(r$ci_high - r$mean, 1.96 * r$mcse)
  expect_equal(r$mean - r$ci_low, 1.96 * r$mcse)
  expect_gt(r$mcse, 0)
})
