# For each true edge, the false independence implied by the graph without it.
deletion_implications <- function(ag) {
  rows <- lapply(seq_len(nrow(ag$edges)), function(i) {
    e <- ag$edges[i, ]
    imp <- shipley_basis(ag$nodes, ag$edges[-i, ])
    hit <- imp[edge_key(imp$X, imp$Y) == edge_key(e$parent, e$child), ]
    cbind(edge = paste(e$parent, "->", e$child), hit)
  })
  do.call(rbind, rows)
}

# Positive control: is each deleted edge's false independence rejected (Bonferroni over edges)?
edge_power <- function(d, del, ag, alpha) {
  r <- test_implications(d, del[, c("X", "Y", "Z")], ag, alpha)
  data.frame(edge = del$edge, p = r$p, detected = !is.na(r$p) & r$p < alpha / nrow(del))
}
