suppressMessages(library(rCausalMGM))

# rCausalMGM needs at least min_category records in each category of a discrete variable;
# sparser binary variables are left out of discovery for that data set.
mgm_data <- function(d, ag, min_category = 5) {
  d <- as.data.frame(d)[, ag$nodes]
  binary <- ag$nodes[ag$scale[ag$nodes] == "binary"]
  sparse <- binary[vapply(binary, function(v) min(table(factor(d[[v]], levels = 0:1))) < min_category, logical(1))]
  d <- d[, setdiff(ag$nodes, sparse), drop = FALSE]
  for (v in intersect(binary, names(d))) d[[v]] <- factor(d[[v]])
  attr(d, "dropped") <- sparse
  d
}

# MGM skeleton (lambda chosen by StEPS, or fixed), then PC-Stable with each collider-orientation rule.
run_mgm_pc <- function(d, params, seed) {
  m <- params$mgm
  if (identical(m$lambda_selection, "steps")) {
    set.seed(seed)
    s <- NULL
    invisible(capture.output(s <- steps(d, numSub = m$steps_subsamples, threads = params$threads)))
    g0 <- s$graph.steps
  } else {
    g0 <- mgm(d, lambda = as.numeric(unlist(m$lambda)))
  }
  rules <- names(params$variants)
  g <- pcStable(d, initialGraph = g0, orientRule = rules, alpha = params$alpha, threads = params$threads)
  graphs <- if (length(rules) == 1) list(g) else g
  list(graphs = setNames(graphs, rules), lambda = as.numeric(g0$lambda))
}

stability_edges <- function(d, g, params) {
  s <- params$stability
  set.seed(s$seed)
  b <- tryCatch(bootstrap(d, g, numBoots = s$subsamples, replace = s$replace, threads = params$threads),
                error = function(e) e)
  if (inherits(b, "error")) return(data.frame(error = conditionMessage(b)))
  out <- as.data.frame(b$stabilities)
  out$adjacency <- 1 - out$none
  out
}

learned_edges <- function(g) {
  t <- as.data.frame(graphTable(g))
  names(t) <- c("from", "to", "type")
  t
}
