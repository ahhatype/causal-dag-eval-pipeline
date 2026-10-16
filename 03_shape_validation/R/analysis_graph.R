analysis_graph <- function(cfg) {
  n <- cfg$nodes
  ev <- cfg$params$era_values
  constant <- names(ev)[vapply(ev, function(v) v[1] == v[2], logical(1))]
  observed <- n$id[n$observed & n$simulate]
  collapsed <- n$id[n$level == "era_fixed" & n$observed]
  nodes <- setdiff(observed, collapsed)

  e <- cfg$edges[cfg$edges$parent %in% observed & cfg$edges$child %in% observed &
                   !cfg$edges$parent %in% constant & !cfg$edges$child %in% constant, ]
  to_era <- function(x) ifelse(x %in% collapsed, "mission_era", x)
  e$source_edge <- paste(e$parent, "->", e$child)
  e$parent <- to_era(e$parent)
  e$child <- to_era(e$child)
  e <- e[e$parent != e$child, ]
  edges <- aggregate(source_edge ~ parent + child, data = e, FUN = function(x) paste(sort(x), collapse = "; "))
  edges <- edges[order(edges$parent, edges$child), ]
  rownames(edges) <- NULL

  list(nodes = topo_order(nodes, edges), edges = edges, collapsed = collapsed, constant = constant,
       scale = setNames(n$scale, n$id)[nodes])
}

as_dagitty <- function(nodes, edges) {
  dagitty::dagitty(sprintf("dag { %s %s }", paste(nodes, collapse = " "),
                           paste(edges$parent, "->", edges$child, collapse = " ")))
}
