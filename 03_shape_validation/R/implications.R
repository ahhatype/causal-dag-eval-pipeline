# Shipley's d-separation basis set: one implication per non-adjacent pair,
# X _||_ Y | pa(X) u pa(Y), with X the later node in topological order.
shipley_basis <- function(nodes, edges) {
  ord <- topo_order(nodes, edges)
  pa <- split(edges$parent, factor(edges$child, levels = ord))
  adj <- paste(edges$parent, edges$child)
  rows <- list()
  for (i in seq_along(ord)[-1]) for (j in seq_len(i - 1)) {
    x <- ord[i]; y <- ord[j]
    if (paste(y, x) %in% adj || paste(x, y) %in% adj) next
    z <- setdiff(union(pa[[x]], pa[[y]]), c(x, y))
    rows[[length(rows) + 1]] <- data.frame(X = x, Y = y, Z = paste(ord[ord %in% z], collapse = "; "))
  }
  if (!length(rows)) return(data.frame(X = character(0), Y = character(0), Z = character(0)))
  do.call(rbind, rows)
}

split_z <- function(z) if (is.na(z) || z == "") character(0) else strsplit(z, "; ", fixed = TRUE)[[1]]

flag_implications <- function(imp, ag, params) {
  g_sel <- dagitty::dagitty(sprintf("dag { %s %s %s }", paste(ag$nodes, collapse = " "),
                                    paste(ag$edges$parent, "->", ag$edges$child, collapse = " "),
                                    paste(params$selection_on, "-> S", collapse = " ")))
  imp$nonadditive_response <- imp$X %in% params$nonadditive_nodes
  imp$expected_under_selection <- vapply(seq_len(nrow(imp)), function(i)
    !dagitty::dseparated(g_sel, imp$X[i], imp$Y[i], c(split_z(imp$Z[i]), "S")), logical(1))
  imp
}

standardize <- function(d, scale) {
  for (v in names(d)) if (scale[[v]] == "continuous") d[[v]] <- as.numeric(scale(d[[v]]))
  d
}

# Likelihood-ratio test of Y in the regression of X on Z (+ Y). Effect: coefficient of Y per SD.
test_implication <- function(d, x, y, z, binary_x) {
  warn <- character(0)
  fit <- function(rhs) withCallingHandlers(
    if (binary_x) glm(reformulate(rhs, x), binomial, data = d) else lm(reformulate(rhs, x), data = d),
    warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (length(unique(d[[x]])) < 2 || length(unique(d[[y]])) < 2)
    return(data.frame(stat = NA, p = NA, effect = NA, fit_warning = "constant variable"))
  d[[y]] <- as.numeric(scale(d[[y]]))
  m0 <- fit(if (length(z)) z else "1")
  m1 <- fit(c(z, y))
  stat <- max(0, 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0))))
  data.frame(stat = stat, p = pchisq(stat, df = 1, lower.tail = FALSE), effect = unname(coef(m1)[y]),
             fit_warning = paste(unique(warn), collapse = "; "))
}

test_implications <- function(d, imp, ag, fdr) {
  d <- standardize(as.data.frame(d)[, ag$nodes], ag$scale)
  res <- do.call(rbind, lapply(seq_len(nrow(imp)), function(i)
    test_implication(d, imp$X[i], imp$Y[i], split_z(imp$Z[i]), ag$scale[[imp$X[i]]] == "binary")))
  out <- cbind(imp, res)
  out$p_bh <- p.adjust(out$p, method = "BH")
  out$rejected <- !is.na(out$p_bh) & out$p_bh < fdr
  out
}
