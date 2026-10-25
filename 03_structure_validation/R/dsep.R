# d-separation by the moralized ancestral graph criterion, on adjacency matrices.
# Equivalent to dagitty::dseparated, without a call into its JavaScript engine per query.

dag_matrix <- function(nodes, edges) {
  a <- matrix(FALSE, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  a[cbind(edges$parent, edges$child)] <- TRUE
  reach <- a
  repeat {
    nxt <- reach | (reach %*% a > 0)
    if (identical(nxt, reach)) break
    reach <- nxt
  }
  list(adj = a, anc = t(reach))            # anc[i, j]: j is an ancestor of i
}

d_separated <- function(g, x, y, z = character(0)) {
  keep <- c(x, y, z)
  an <- union(keep, colnames(g$anc)[colSums(g$anc[keep, , drop = FALSE]) > 0])
  a <- g$adj[an, an, drop = FALSE]
  m <- a | t(a)
  co <- crossprod(t(a) * 1) > 0             # co[i, j]: i and j share a child
  m <- (m | co) & !diag(length(an))
  open <- setdiff(an, z)
  m <- m[open, open, drop = FALSE]
  seen <- setNames(open == x, open)
  frontier <- seen
  repeat {
    nxt <- (as.numeric(frontier) %*% m > 0)[1, ] & !seen
    if (!any(nxt)) break
    seen <- seen | nxt
    if (seen[[y]]) return(FALSE)
    frontier <- nxt
  }
  !seen[[y]]
}

holds_in <- function(imp, g) {
  vapply(seq_len(nrow(imp)), function(i) d_separated(g, imp$X[i], imp$Y[i], split_z(imp$Z[i])), logical(1))
}
