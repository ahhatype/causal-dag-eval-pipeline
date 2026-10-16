validate_dag <- function(cfg_dir, cfg) {
  g <- dagitty::dagitty(paste(readLines(file.path(cfg_dir, "dag.txt")), collapse = "\n"))
  e <- dagitty::edges(g)
  ref <- sort(paste(slug_id(e$v), slug_id(e$w)))
  all_edges <- read.csv(file.path(cfg_dir, "edges.csv"), stringsAsFactors = FALSE)
  mine <- sort(paste(all_edges$parent, all_edges$child))
  nm <- sort(slug_id(names(g)))
  checks <- list(
    n_nodes = length(nm) == 54 && identical(nm, sort(cfg$nodes$id)),
    n_edges = length(ref) == 90,
    same_edges = identical(ref, mine),
    acyclic = dagitty::isAcyclic(g)
  )
  anc <- slug_id(dagitty::ancestors(g, "Nephrolithiasis"))
  n <- cfg$nodes
  keep <- n$observed & n$simulate & n$level != "era_fixed"
  anc_set <- setdiff(n$id[keep & n$id %in% anc], "nephrolithiasis")
  checks$ancestor_set_16 <- setequal(anc_set, n$id[n$feature_set == "ancestor"]) && length(anc_set) == 16
  other <- setdiff(n$id[keep], c(anc_set, "nephrolithiasis"))
  checks$non_ancestor_19 <- setequal(other, n$id[n$feature_set == "non_ancestor"]) && length(other) == 19
  checks$latent_7 <- sum(!n$observed) == 7
  bad <- names(checks)[!unlist(checks)]
  if (length(bad)) stop("DAG validation failed: ", paste(bad, collapse = ", "))
  invisible(checks)
}
