outcome_prob <- function(df, cfg, cal) {
  plogis(cal$logit[[outcome_id]] + eval(parse(text = lin_expr(outcome_id, cfg)), envir = df))
}

ground_truth <- function(D, cfg, cal) {
  p <- cfg$params
  n <- p$sizes$ground_truth; seed <- p$seeds$ground_truth
  base <- sim_source(D, n, seed)
  keep <- retention_mask(base, cfg, cal)
  q <- c(p$ground_truth$low_quantile, p$ground_truth$high_quantile)
  anc <- feature_ids(cfg, "ancestor")
  scale <- setNames(cfg$nodes$scale, cfg$nodes$id)
  eff <- function(r, lo, hi, m) {
    d <- outcome_prob(r[[hi]], cfg, cal)[m] - outcome_prob(r[[lo]], cfg, cal)[m]
    c(mean(d), sd(d) / sqrt(length(d)))
  }
  rows <- lapply(feature_ids(cfg, "all"), function(f) {
    if (scale[[f]] == "binary") {
      vf <- va <- c(0, 1)
      r <- sim_actions(D, cfg, n, seed, list(lo = list(node = f, value = 0), hi = list(node = f, value = 1)))
      ef <- eff(r, "lo", "hi", rep(TRUE, n)); ea <- eff(r, "lo", "hi", keep)
    } else {
      vf <- unname(quantile(base[[f]], q)); va <- unname(quantile(base[[f]][keep], q))
      r <- sim_actions(D, cfg, n, seed, list(
        lo_full = list(node = f, value = vf[1]), hi_full = list(node = f, value = vf[2]),
        lo_ast = list(node = f, value = va[1]), hi_ast = list(node = f, value = va[2])))
      ef <- eff(r, "lo_full", "hi_full", rep(TRUE, n)); ea <- eff(r, "lo_ast", "hi_ast", keep)
    }
    data.frame(feature = f, ancestor_set = f %in% anc,
               low_full = vf[1], high_full = vf[2], low_astronaut = va[1], high_astronaut = va[2],
               effect_full = ef[1], se_full = ef[2], effect_astronaut = ea[1], se_astronaut = ea[2])
  })
  out <- do.call(rbind, rows)
  out$rank_full <- rank(-abs(out$effect_full), ties.method = "min")
  out$rank_astronaut <- rank(-abs(out$effect_astronaut), ties.method = "min")
  out
}

unresolved_pairs <- function(gt, target = c("full", "astronaut")) {
  target <- match.arg(target)
  e <- gt[[paste0("effect_", target)]]; s <- gt[[paste0("se_", target)]]
  a <- gt[gt$ancestor_set, ]; ea <- abs(e[gt$ancestor_set]); sa <- s[gt$ancestor_set]
  o <- order(-ea); a <- a[o, ]; ea <- ea[o]; sa <- sa[o]
  i <- which(diff(-ea) <= 2 * sqrt(sa[-1]^2 + sa[-length(sa)]^2))
  if (!length(i)) return(character(0))
  paste(a$feature[i], a$feature[i + 1], sep = " ~ ")
}
