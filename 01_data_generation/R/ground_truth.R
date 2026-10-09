# Ground truth evaluates the node formulas directly (the same strings simcausal uses), so an intervention
# re-evaluates only the target's descendants and every other draw is shared (common random numbers).

outcome_prob <- function(df, cfg, cal) {
  plogis(cal$logit[[outcome_id]] + eval(parse(text = lin_expr(outcome_id, cfg)), envir = df))
}

compiled_formulas <- function(cfg, cal) {
  setNames(lapply(cfg$order, function(id) parse(text = node_formula(id, cfg, cal))[[1]]), cfg$order)
}

# Nodes are evaluated in topological order into an environment holding the noise draws.
truth_population <- function(cfg, cal, n, seed, forms = compiled_formulas(cfg, cal)) {
  set.seed(seed)
  env <- new.env(parent = globalenv())
  for (id in noise_ids(cfg)) assign(paste0("U_", id), runif(n), envir = env)
  for (id in cfg$order) assign(id, eval(forms[[id]], envir = env), envir = env)
  env
}

descendants <- function(id, edges) {
  out <- character(0); frontier <- id
  while (length(frontier)) {
    nxt <- setdiff(edges$child[edges$parent %in% frontier], out)
    out <- c(out, nxt); frontier <- nxt
  }
  out
}

ancestors <- function(id, edges) {
  out <- character(0); frontier <- id
  while (length(frontier)) {
    nxt <- setdiff(edges$parent[edges$child %in% frontier], out)
    out <- c(out, nxt); frontier <- nxt
  }
  out
}

# Outcome probability per record under do(node = value).
do_outcome_prob <- function(env, cfg, cal, forms, node, value) {
  need <- intersect(cfg$order, intersect(descendants(node, cfg$edges), c(ancestors(outcome_id, cfg$edges), outcome_id)))
  child <- new.env(parent = env)
  assign(node, rep_len(value, length(env[[node]])), envir = child)
  for (id in setdiff(need, outcome_id)) assign(id, eval(forms[[id]], envir = child), envir = child)
  outcome_prob(child, cfg, cal)
}

# Weighted mean absolute deviation of mu around its weighted centre; zero-weight points contribute nothing.
pop_spread <- function(mu, w) {
  m <- sum(w * mu)
  sum(w * abs(mu - m))
}

# Grid for the population-scaled truth within one era: binary values weighted by within-era prevalence,
# continuous values at the midpoints of equal-probability bins, weighted equally.
truth_grid <- function(x, scale, points) {
  if (scale == "binary") {
    p <- mean(x)
    list(x = c(0, 1), w = c(1 - p, p))
  } else {
    list(x = unname(quantile(x, (seq_len(points) - 0.5) / points)), w = rep(1 / points, points))
  }
}

batch_mcse <- function(values) sd(values) / sqrt(length(values))

# Population-scaled (pop) and record-level (rec) truth for one feature in one population.
# P is a list of per-record outcome probabilities, one per grid value, named by grid key.
population_truths <- function(f, scale, base_x, era, idx, batch, points, eval_at) {
  eras <- if (f == "mission_era") list(all = rep(TRUE, length(era))) else list(`0` = era == 0, `1` = era == 1)
  grids <- lapply(eras, function(m) {
    if (f == "mission_era") list(x = c(0, 1), w = c(1 - mean(era[idx]), mean(era[idx])))
    else truth_grid(base_x[idx & m], scale, points)
  })
  P <- lapply(grids, function(g) lapply(g$x, eval_at))
  share <- vapply(eras, function(m) mean(m[idx]), numeric(1))
  pop_value <- function(keep) {
    sum(vapply(seq_along(eras), function(k) {
      mu <- vapply(P[[k]], function(p) mean(p[keep & eras[[k]]]), numeric(1))
      if (!any(keep & eras[[k]])) return(0)
      mean(eras[[k]][keep]) * pop_spread(mu, grids[[k]]$w)
    }, numeric(1)))
  }
  centre <- numeric(length(era))
  for (k in seq_along(eras)) {
    ck <- Reduce(`+`, Map(function(p, w) w * p, P[[k]], grids[[k]]$w))
    centre[eras[[k]]] <- ck[eras[[k]]]
  }
  list(grids = grids, P = P, share = share, eras = eras,
       pop = c(pop_value(idx), batch_mcse(vapply(sort(unique(batch)), function(b) pop_value(idx & batch == b), numeric(1)))),
       centre = centre)
}

ground_truth <- function(cfg, cal) {
  p <- cfg$params
  g <- p$ground_truth
  n <- p$sizes$ground_truth
  forms <- compiled_formulas(cfg, cal)
  env <- truth_population(cfg, cal, n, p$seeds$ground_truth, forms)
  keep <- retention_mask(env, cfg, cal)
  pops <- list(source = rep(TRUE, n), selected = keep)
  batch <- (seq_len(n) - 1L) %% g$batches + 1L
  era <- env$mission_era
  p_obs <- outcome_prob(env, cfg, cal)
  anc <- feature_ids(cfg, "ancestor")
  scale <- setNames(cfg$nodes$scale, cfg$nodes$id)
  cache <- new.env()

  rows <- list(); checks <- list()
  for (f in feature_ids(cfg, "all")) {
    base_x <- env[[f]]
    eval_at <- function(x) {
      key <- paste(f, format(x, digits = 17))
      if (is.null(cache[[key]])) cache[[key]] <- do_outcome_prob(env, cfg, cal, forms, f, x)
      cache[[key]]
    }
    v <- if (scale[[f]] == "binary") c(0, 1) else unname(quantile(base_x, c(g$low_quantile, g$high_quantile)))
    d <- eval_at(v[2]) - eval_at(v[1])
    for (pn in names(pops)) {
      idx <- pops[[pn]]
      rows[[length(rows) + 1]] <- data.frame(feature = f, truth_type = "per_unit", population = pn,
                                             value = mean(d[idx]), mcse = sd(d[idx]) / sqrt(sum(idx)), low = v[1], high = v[2])
      if (!(f %in% anc)) {
        for (tt in c("pop", "rec"))
          rows[[length(rows) + 1]] <- data.frame(feature = f, truth_type = tt, population = pn, value = 0, mcse = 0,
                                                 low = NA, high = NA)
        next
      }
      pt <- population_truths(f, scale[[f]], base_x, era, idx, batch, g$grid_points, eval_at)
      r <- abs(p_obs - pt$centre)
      rows[[length(rows) + 1]] <- data.frame(feature = f, truth_type = c("pop", "rec"), population = pn,
                                             value = c(pt$pop[1], mean(r[idx])), mcse = c(pt$pop[2], sd(r[idx]) / sqrt(sum(idx))),
                                             low = NA, high = NA)
      if (scale[[f]] == "binary" && f != "mission_era")
        checks[[length(checks) + 1]] <- binary_check(f, pn, pt, base_x, p_obs, idx, batch)
    }
    rm(list = ls(cache), envir = cache)
  }
  out <- do.call(rbind, rows)
  out$ancestor_set <- out$feature %in% anc
  list(values = out, binary_checks = do.call(rbind, checks), retained = sum(keep), n = n)
}

# Within era, T_pop equals 2p(1 - p)|interventional contrast| exactly. Replacing that contrast with the observational
# one shows how far confounding moves the naive value; a gap beyond 3 MCSE flags a confounded feature.
binary_check <- function(f, pn, pt, base_x, p_obs, idx, batch) {
  do.call(rbind, lapply(seq_along(pt$eras), function(k) {
    m <- idx & pt$eras[[k]]
    prev <- mean(base_x[m])
    t_era <- pop_spread(vapply(pt$P[[k]], function(q) mean(q[m]), numeric(1)), pt$grids[[k]]$w)
    naive <- function(mm) {
      if (!any(mm & base_x == 1) || !any(mm & base_x == 0)) return(NA_real_)
      pr <- mean(base_x[mm])
      2 * pr * (1 - pr) * abs(mean(p_obs[mm & base_x == 1]) - mean(p_obs[mm & base_x == 0]))
    }
    tb <- vapply(sort(unique(batch)), function(b) {
      mb <- m & batch == b
      naive(mb) - pop_spread(vapply(pt$P[[k]], function(q) mean(q[mb]), numeric(1)), c(1 - mean(base_x[mb]), mean(base_x[mb])))
    }, numeric(1))
    diff <- naive(m) - t_era
    mc <- batch_mcse(tb)
    data.frame(feature = f, population = pn, era = names(pt$eras)[k], prevalence = prev, t_pop_era = t_era,
               observational_form = naive(m), difference = diff, mcse = mc,
               flagged = isTRUE(abs(diff) > 3 * mc))
  }))
}

# Adjacent per-unit truths (ancestor set, |value| ordered) within 2 combined MCSE.
unresolved_pairs <- function(truth, population = c("source", "selected")) {
  population <- match.arg(population)
  a <- truth[truth$ancestor_set & truth$truth_type == "per_unit" & truth$population == population, ]
  a <- a[order(-abs(a$value)), ]
  e <- abs(a$value); s <- a$mcse
  i <- which(diff(-e) <= 2 * sqrt(s[-1]^2 + s[-length(s)]^2))
  if (!length(i)) return(character(0))
  paste(a$feature[i], a$feature[i + 1], sep = " ~ ")
}
