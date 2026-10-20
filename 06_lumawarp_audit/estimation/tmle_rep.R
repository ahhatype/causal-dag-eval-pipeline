# One replicate: Super Learner on the same folds, then tmle() with the LumaWarp and the Super Learner components,
# plus the naive difference and a logistic g-computation. Usage: Rscript tmle_rep.R <tmle_inputs.csv> <out.csv> <treatment> <outcome> <truth>
suppressPackageStartupMessages({library(tmle); library(SuperLearner)})
args <- commandArgs(trailingOnly = TRUE); inp <- args[1]; outp <- args[2]; treatment <- args[3]; outcome <- args[4]; truth <- as.numeric(args[5])
d <- read.csv(inp, check.names = FALSE); gbound <- c(0.025, 0.975)
Wn <- setdiff(names(d), c("id", "fold", treatment, outcome, "g_lw", "Q0_lw", "Q1_lw")); W <- d[, Wn, drop = FALSE]; A <- d[[treatment]]; Y <- d[[outcome]]
seed <- as.integer(sub("^rep", "", basename(dirname(inp))))
d$g_sl <- NA_real_; d$Q0_sl <- NA_real_; d$Q1_sl <- NA_real_; t0 <- proc.time()[["elapsed"]]
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
sl_seconds <- proc.time()[["elapsed"]] - t0
one <- function(suffix) {
  fit <- tmle(Y = Y, A = A, W = W, Q = cbind(Q0W = d[[paste0("Q0_", suffix)]], Q1W = d[[paste0("Q1_", suffix)]]), g1W = d[[paste0("g_", suffix)]],
              family = "binomial", fluctuation = "logistic", gbound = gbound, prescreenW.g = FALSE, target.gwt = FALSE, evalATT = FALSE)
  e <- fit$estimates$ATE; c(est = unname(e$psi), se = sqrt(unname(e$var.psi)), lower = unname(e$CI[1]), upper = unname(e$CI[2]))
}
lw <- one("lw"); sl <- one("sl")
naive <- mean(Y[A == 1]) - mean(Y[A == 0])
gl <- glm(Y ~ ., data = cbind(Y = Y, XA), family = binomial()); x0 <- XA; x0[[treatment]] <- 0; x1 <- XA; x1[[treatment]] <- 1
gcomp <- mean(predict(gl, x1, type = "response")) - mean(predict(gl, x0, type = "response"))
nuis <- function(suffix) c(g_auc_proxy = cor(d[[paste0("g_", suffix)]], A), g_clipped = mean(d[[paste0("g_", suffix)]] < gbound[1] | d[[paste0("g_", suffix)]] > gbound[2]))
row <- data.frame(treatment = treatment, seed = seed, truth = truth, n = nrow(d), events = sum(Y), treated = sum(A),
                  lw_est = lw["est"], lw_se = lw["se"], lw_lower = lw["lower"], lw_upper = lw["upper"], lw_covered = (lw["lower"] <= truth & truth <= lw["upper"]),
                  sl_est = sl["est"], sl_se = sl["se"], sl_lower = sl["lower"], sl_upper = sl["upper"], sl_covered = (sl["lower"] <= truth & truth <= sl["upper"]),
                  naive_est = naive, gcomp_logistic_est = gcomp, g_lw_clipped = nuis("lw")["g_clipped"], g_sl_clipped = nuis("sl")["g_clipped"], sl_seconds = sl_seconds)
dir.create(dirname(outp), recursive = TRUE, showWarnings = FALSE); tmp <- paste0(outp, ".tmp"); write.csv(row, tmp, row.names = FALSE); file.rename(tmp, outp)
