v_structures <- function(edges, adjacent) {
  out <- character(0)
  for (ch in unique(edges$child)) {
    pa <- sort(edges$parent[edges$child == ch])
    if (length(pa) < 2) next
    for (pr in combn(pa, 2, simplify = FALSE))
      if (!paste(sort(pr), collapse = "|") %in% adjacent) out <- c(out, paste(pr[1], ch, pr[2]))
  }
  sort(out)
}

# Dor & Tarsi (1992): a consistent extension of a partially directed graph, or NULL if none exists.
consistent_extension <- function(nodes, learned) {
  dir <- learned[learned$type == "dir", c("from", "to")]
  und <- learned[learned$type != "dir", c("from", "to")]
  out <- data.frame(parent = dir$from, child = dir$to)
  left <- nodes
  while (length(left)) {
    nb_und <- function(x) c(und$to[und$from == x], und$from[und$to == x])
    adj <- function(x) unique(c(nb_und(x), dir$to[dir$from == x], dir$from[dir$to == x]))
    pick <- NULL
    for (x in left) {
      if (any(dir$from == x)) next
      u <- nb_und(x)
      others <- adj(x)
      if (all(vapply(u, function(y) all(setdiff(others, y) %in% c(adj(y), y)), logical(1)))) { pick <- x; break }
    }
    if (is.null(pick)) return(NULL)
    u <- nb_und(pick)
    out <- rbind(out, data.frame(parent = u, child = rep(pick, length(u))))
    dir <- dir[dir$to != pick, ]
    und <- und[und$from != pick & und$to != pick, ]
    left <- setdiff(left, pick)
  }
  out
}

# All acyclic orientations of a CPDAG's undirected edges (sampled beyond max_enumerated), with a flag
# for consistent extensions: those adding no v-structure beyond the CPDAG's directed edges.
dag_orientations <- function(nodes, learned, max_enumerated = 4096, seed = 1) {
  dir <- learned[learned$type == "dir", ]
  und <- learned[learned$type != "dir", ]
  base <- data.frame(parent = dir$from, child = dir$to)
  adjacent <- apply(learned[, c("from", "to")], 1, function(r) paste(sort(r), collapse = "|"))
  v0 <- v_structures(base, adjacent)
  k <- nrow(und)
  bits <- if (2^k <= max_enumerated) {
    as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), k)))
  } else {
    set.seed(seed)
    matrix(runif(max_enumerated * k) < 0.5, ncol = k)
  }
  dags <- list(); consistent <- logical(0)
  for (i in seq_len(max(1, nrow(bits)))) {
    flip <- if (k) bits[i, ] else logical(0)
    e <- rbind(base, data.frame(parent = ifelse(flip, und$to, und$from), child = ifelse(flip, und$from, und$to)))
    if (!dagitty::isAcyclic(as_dagitty(nodes, e))) next
    dags[[length(dags) + 1]] <- e
    consistent <- c(consistent, identical(v_structures(e, adjacent), v0))
  }
  list(dags = dags, consistent = consistent, undirected = k)
}

# Comparison over the consistent extension if one exists (all are Markov equivalent), otherwise
# summarized over all acyclic orientations, keeping the median orientation for the detail files.
compare_learned <- function(nodes, truth_edges, learned, params, imp_t = shipley_basis(nodes, truth_edges)) {
  ext <- consistent_extension(nodes, learned)
  ors <- if (!is.null(ext)) {
    list(dags = list(ext), consistent = TRUE, undirected = sum(learned$type != "dir"))
  } else {
    dag_orientations(nodes, learned, params$orientations$max, params$orientations$seed)
  }
  use <- if (any(ors$consistent)) which(ors$consistent)[1] else seq_along(ors$dags)
  if (!length(use)) {
    na <- data.frame(truth_implications = nrow(imp_t), learned_implications = NA, truth_held_by_learned = NA,
                     learned_held_by_truth = NA, truth_only = NA, learned_only = NA,
                     truth_held_by_learned_min = NA, truth_held_by_learned_max = NA,
                     learned_held_by_truth_min = NA, learned_held_by_truth_max = NA)
    return(list(summary = cbind(data.frame(consistent_extension = FALSE, undirected_edges = ors$undirected,
                                           orientations_compared = 0L), na), detail = NULL))
  }
  cmps <- lapply(ors$dags[use], function(e) compare_implications(nodes, truth_edges, e, imp_t))
  s <- do.call(rbind, lapply(cmps, `[[`, "summary"))
  med <- order(s$truth_held_by_learned)[ceiling(nrow(s) / 2)]
  summary <- cbind(
    data.frame(consistent_extension = any(ors$consistent), undirected_edges = ors$undirected,
               orientations_compared = length(use)),
    s[med, ],
    data.frame(truth_held_by_learned_min = min(s$truth_held_by_learned),
               truth_held_by_learned_max = max(s$truth_held_by_learned),
               learned_held_by_truth_min = min(s$learned_held_by_truth),
               learned_held_by_truth_max = max(s$learned_held_by_truth)))
  rownames(summary) <- NULL
  list(summary = summary, detail = cmps[[med]])
}

# Implications of each graph checked for d-separation in the other.
compare_implications <- function(nodes, truth_edges, learned_dag_edges, imp_t = shipley_basis(nodes, truth_edges)) {
  imp_l <- shipley_basis(nodes, learned_dag_edges)
  t_in_l <- holds_in(imp_t, dag_matrix(nodes, learned_dag_edges))
  l_in_t <- holds_in(imp_l, dag_matrix(nodes, truth_edges))
  list(summary = data.frame(
         truth_implications = nrow(imp_t), learned_implications = nrow(imp_l),
         truth_held_by_learned = sum(t_in_l), learned_held_by_truth = sum(l_in_t),
         truth_only = sum(!t_in_l), learned_only = sum(!l_in_t)),
       truth = cbind(imp_t, held_by_other = t_in_l), learned = cbind(imp_l, held_by_other = l_in_t))
}

# The true graph's equivalence class: compelled edges are directed, reversible edges undirected.
true_cpdag <- function(nodes, truth_edges) {
  e <- dagitty::edges(dagitty::equivalenceClass(as_dagitty(nodes, truth_edges)))
  data.frame(from = e$v, to = e$w, type = ifelse(e$e == "->", "dir", "undir"))
}

edge_key <- function(a, b) ifelse(a < b, paste(a, b), paste(b, a))

# Edge mark relative to the canonical (sorted) pair: "fwd", "rev" or "undir".
edge_mark <- function(from, to, type) ifelse(type != "dir", "undir", ifelse(from < to, "fwd", "rev"))

# Adjacency precision/recall, and orientation agreement with the true CPDAG on shared adjacencies.
structure_metrics <- function(cpdag, learned) {
  t <- edge_key(cpdag$from, cpdag$to)
  l <- edge_key(learned$from, learned$to)
  tm <- setNames(edge_mark(cpdag$from, cpdag$to, cpdag$type), t)
  lm <- setNames(edge_mark(learned$from, learned$to, learned$type), l)
  shared <- intersect(l, t)
  compelled <- shared[tm[shared] != "undir"]
  data.frame(learned_edges = length(l), true_edges = length(t),
             adjacency_precision = if (length(l)) mean(l %in% t) else NA, adjacency_recall = mean(t %in% l),
             shared_edges = length(shared),
             orientation_agreement = if (length(shared)) mean(lm[shared] == tm[shared]) else NA,
             compelled_shared = length(compelled),
             compelled_correct = sum(lm[compelled] == tm[compelled]),
             compelled_reversed = sum(lm[compelled] != "undir" & lm[compelled] != tm[compelled]),
             compelled_unoriented = sum(lm[compelled] == "undir"))
}
