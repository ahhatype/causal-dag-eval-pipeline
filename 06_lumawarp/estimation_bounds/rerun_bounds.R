# Re-targets one folder-34 unit under several propensity bounds. Usage: Rscript rerun_bounds.R <tmle_inputs.csv> <out.csv> <treatment> <outcome> <truth> <results_csv_from_34>
suppressPackageStartupMessages(library(tmle))
args <- commandArgs(trailingOnly = TRUE); inp <- args[1]; outp <- args[2]; treatment <- args[3]; outcome <- args[4]; truth <- as.numeric(args[5]); r34 <- args[6]
d <- read.csv(inp, check.names = FALSE); prev <- read.csv(r34)
Wn <- setdiff(names(d), c("id", "fold", treatment, outcome, "g_lw", "Q0_lw", "Q1_lw")); W <- d[, Wn, drop = FALSE]; A <- d[[treatment]]; Y <- d[[outcome]]
seed <- as.integer(sub("^rep", "", basename(dirname(inp))))
# the Super Learner nuisances were not saved in folder 34's inputs; refit them identically (same seeds, same library) so both learners are compared
suppressPackageStartupMessages(library(SuperLearner))
d$g_sl <- NA_real_; d$Q0_sl <- NA_real_; d$Q1_sl <- NA_real_
for (k in sort(unique(d$fold))) {
  tr <- d$fold != k; te <- d$fold == k
  set.seed(seed * 10 + 100 * k + 50)
  gf <- SuperLearner(Y = A[tr], X = W[tr, , drop = FALSE], family = binomial(), SL.library = c("SL.glm", "SL.glmnet", "SL.earth"), cvControl = list(V = 5L))
  d$g_sl[te] <- as.numeric(predict(gf, newdata = W[te, , drop = FALSE], onlySL = TRUE)$pred)
  set.seed(seed * 10 + 100 * k + 51); XA <- cbind(setNames(data.frame(A), treatment), W)
  qf <- SuperLearner(Y = Y[tr], X = XA[tr, , drop = FALSE], family = binomial(), SL.library = c("SL.glm", "SL.glmnet", "SL.earth"), cvControl = list(V = 5L))
  x0 <- XA[te, , drop = FALSE]; x0[[treatment]] <- 0; x1 <- x0; x1[[treatment]] <- 1
  d$Q0_sl[te] <- as.numeric(predict(qf, newdata = x0, onlySL = TRUE)$pred); d$Q1_sl[te] <- as.numeric(predict(qf, newdata = x1, onlySL = TRUE)$pred)
}
rows <- list()
for (lo in c(0.025, 0.05, 0.01, 0.001)) for (suffix in c("lw", "sl")) {
  gb <- c(lo, 1 - lo)
  fit <- tmle(Y = Y, A = A, W = W, Q = cbind(Q0W = d[[paste0("Q0_", suffix)]], Q1W = d[[paste0("Q1_", suffix)]]), g1W = d[[paste0("g_", suffix)]],
              family = "binomial", fluctuation = "logistic", gbound = gb, prescreenW.g = FALSE, target.gwt = FALSE, evalATT = FALSE)
  e <- fit$estimates$ATE
  rows[[length(rows) + 1]] <- data.frame(treatment = treatment, seed = seed, truth = truth, learner = suffix, bound = lo, est = unname(e$psi), se = sqrt(unname(e$var.psi)),
                                         lower = unname(e$CI[1]), upper = unname(e$CI[2]), covered = (e$CI[1] <= truth & truth <= e$CI[2]),
                                         clipped = mean(d[[paste0("g_", suffix)]] < gb[1] | d[[paste0("g_", suffix)]] > gb[2]))
}
out <- do.call(rbind, rows)
# consistency check against folder 34 at the original bound
chk <- out[out$bound == 0.025 & out$learner == "lw", ]; out$check_vs_34 <- abs(chk$est - prev$lw_est[1])
dir.create(dirname(outp), recursive = TRUE, showWarnings = FALSE); tmp <- paste0(outp, ".tmp"); write.csv(out, tmp, row.names = FALSE); file.rename(tmp, outp)
