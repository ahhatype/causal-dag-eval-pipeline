test_that("settings parse to the expected types", {
  expect_true(is.numeric(params$stability$subsamples))
  expect_false(params$stability$replace)
  expect_true(params$mgm$lambda_selection %in% c("steps", "fixed"))
  expect_length(unlist(params$mgm$lambda), 3)
  expect_true(all(params$nonadditive_nodes %in% cfg$nodes$id))
  expect_true(all(params$selection_on %in% cfg$nodes$id))
})

test_that("analysis graph: 36 observed variables, 59 edges, acyclic, era-fixed nodes collapsed", {
  expect_length(ag$nodes, 36)
  expect_equal(nrow(ag$edges), 59)
  expect_true(dagitty::isAcyclic(as_dagitty(ag$nodes, ag$edges)))
  expect_false(any(c(ag$collapsed, ag$constant) %in% c(ag$nodes, ag$edges$parent, ag$edges$child)))
  expect_setequal(ag$constant, c("altered_gravity", "humidity"))
  expect_true(all(c("ancestor", "non_ancestor") %in% cfg$nodes$feature_set[cfg$nodes$id %in% ag$nodes]))
})

test_that("Shipley basis set: one implication per non-adjacent pair, all implied by the graph", {
  imp <- shipley_basis(ag$nodes, ag$edges)
  expect_equal(nrow(imp), choose(36, 2) - 59)
  expect_false(anyDuplicated(paste(pmin(imp$X, imp$Y), pmax(imp$X, imp$Y))) > 0)
  expect_true(all(holds_in(imp, dag_matrix(ag$nodes, ag$edges))))
})

test_that("only individual factors vs pre-flight fitness is expected to fail under selection", {
  imp <- flag_implications(shipley_basis(ag$nodes, ag$edges), ag, params)
  sel <- imp[imp$expected_under_selection, ]
  expect_equal(nrow(sel), 1)
  expect_setequal(c(sel$X, sel$Y), params$selection_on)
})

test_that("likelihood-ratio test keeps a true independence and rejects a planted dependence", {
  set.seed(1)
  n <- 3000
  a <- rnorm(n); c <- rnorm(n)
  b <- 0.5 * a + rnorm(n)
  d <- data.frame(a = a, b = b, c = c, y = rbinom(n, 1, plogis(0.8 * a)))
  expect_gt(test_implication(d, "b", "c", "a", FALSE)$p, 0.01)
  d$b2 <- 0.5 * a + 0.3 * c + rnorm(n)
  expect_lt(test_implication(d, "b2", "c", "a", FALSE)$p, 1e-6)
  expect_gt(test_implication(d, "y", "c", "a", TRUE)$p, 0.01)
})

test_that("implication comparison detects an added and an omitted edge", {
  nodes <- c("a", "b", "c")
  truth <- data.frame(parent = c("a", "b"), child = c("b", "c"))
  added <- data.frame(parent = c("a", "b", "a"), child = c("b", "c", "c"))
  omitted <- data.frame(parent = "a", child = "b")
  s1 <- compare_implications(nodes, truth, added)$summary
  expect_equal(c(s1$truth_only, s1$learned_only), c(1, 0))
  s2 <- compare_implications(nodes, truth, omitted)$summary
  expect_equal(s2$truth_only, 0)
  expect_equal(s2$learned_only, 2)
})

test_that("fast d-separation agrees with dagitty", {
  imp <- shipley_basis(ag$nodes, ag$edges)
  g <- as_dagitty(ag$nodes, ag$edges)
  set.seed(3)
  probe <- data.frame(X = sample(ag$nodes, 200, TRUE), Y = sample(ag$nodes, 200, TRUE))
  probe <- probe[probe$X != probe$Y, ]
  probe$Z <- vapply(seq_len(nrow(probe)), function(i)
    paste(sample(setdiff(ag$nodes, c(probe$X[i], probe$Y[i])), sample(0:4, 1)), collapse = "; "), "")
  both <- rbind(imp, probe)
  slow <- vapply(seq_len(nrow(both)), function(i)
    dagitty::dseparated(g, both$X[i], both$Y[i], split_z(both$Z[i])), logical(1))
  expect_identical(holds_in(both, dag_matrix(ag$nodes, ag$edges)), slow)
})

test_that("consistent extension: found when one exists, NULL for an unchorded undirected 4-cycle", {
  chain <- data.frame(from = c("a", "b"), to = c("b", "c"), type = c("undir", "undir"))
  ext <- consistent_extension(c("a", "b", "c"), chain)
  expect_equal(nrow(ext), 2)
  expect_length(v_structures(ext, c("a|b", "b|c")), 0)
  forced <- data.frame(from = c("a", "b"), to = c("b", "c"), type = c("dir", "undir"))
  expect_true(all(consistent_extension(c("a", "b", "c"), forced) == data.frame(parent = c("a", "b"), child = c("b", "c"))))
  cycle <- data.frame(from = c("a", "b", "c", "d"), to = c("b", "c", "d", "a"), type = "undir")
  expect_null(consistent_extension(c("a", "b", "c", "d"), cycle))
})

test_that("comparison is deterministic and reports a range when no consistent extension exists", {
  nodes <- c("a", "b", "c", "d")
  truth <- data.frame(parent = c("a", "b", "c"), child = c("b", "c", "d"))
  cycle <- data.frame(from = c("a", "b", "c", "d"), to = c("b", "c", "d", "a"), type = "undir")
  p <- list(orientations = list(max = 4096, seed = 1))
  r1 <- compare_learned(nodes, truth, cycle, p)
  expect_false(r1$summary$consistent_extension)
  expect_identical(r1, compare_learned(nodes, truth, cycle, p))
  expect_true(r1$summary$truth_held_by_learned_min <= r1$summary$truth_held_by_learned)
  expect_true(r1$summary$truth_held_by_learned <= r1$summary$truth_held_by_learned_max)
  expect_gt(r1$summary$orientations_compared, 1)
})

test_that("orientation is scored against the true CPDAG", {
  nodes <- c("a", "b", "c")
  chain <- true_cpdag(nodes, data.frame(parent = c("a", "b"), child = c("b", "c")))
  expect_true(all(chain$type == "undir"))
  m1 <- structure_metrics(chain, data.frame(from = c("a", "b"), to = c("b", "c"), type = "dir"))
  expect_equal(c(m1$adjacency_recall, m1$orientation_agreement, m1$compelled_shared), c(1, 0, 0))
  collider <- true_cpdag(nodes, data.frame(parent = c("a", "c"), child = c("b", "b")))
  m2 <- structure_metrics(collider, data.frame(from = c("a", "b"), to = c("b", "c"), type = "dir"))
  expect_equal(c(m2$compelled_shared, m2$compelled_correct, m2$compelled_reversed), c(2, 1, 1))
})

test_that("edge-deletion positive control: one false independence per true edge, detected when the edge is strong", {
  del <- deletion_implications(ag)
  expect_equal(nrow(del), nrow(ag$edges))
  g <- dag_matrix(ag$nodes, ag$edges)
  expect_false(any(holds_in(del, g)))
  set.seed(4)
  n <- 2000
  a <- rnorm(n); b <- 0.5 * a + rnorm(n)
  toy <- list(nodes = c("a", "b"), edges = data.frame(parent = "a", child = "b"), scale = c(a = "continuous", b = "continuous"))
  pw <- edge_power(data.frame(a = a, b = b), deletion_implications(toy), toy, 0.05)
  expect_true(pw$detected)
})

test_that("MGM-PC-Stable wrapper recovers the skeleton of a small mixed-data chain", {
  set.seed(2)
  n <- 2000
  a <- rnorm(n); b <- 0.8 * a + rnorm(n)
  c <- rbinom(n, 1, plogis(1.5 * b))
  d <- data.frame(a = a, b = b, c = factor(c), e = rnorm(n))
  learned <- learned_edges(run_mgm_pc(d, modifyList(params, list(threads = 1)), seed = 1)$graphs[[1]])
  keys <- apply(learned[, 1:2], 1, function(r) paste(sort(r), collapse = "-"))
  expect_setequal(keys, c("a-b", "b-c"))
})

test_that("stale data are refused", {
  tmp <- file.path(tempdir(), "stale_data"); dir.create(tmp, showWarnings = FALSE)
  writeLines(c("config hash: not-the-current-hash", "replicates: 2"), file.path(tmp, "provenance.txt"))
  cal <- read_calibration(file.path(repo_dir, "01_data_generation", "config", "calibration.yaml"))
  expect_error(check_data_current(tmp, cal), "make data")
  writeLines(c(paste("config hash:", cal$config_hash), "replicates: 2"), file.path(tmp, "provenance.txt"))
  expect_equal(check_data_current(tmp, cal), 2L)
})

test_that("shape pipeline runs end to end on generated replicates", {
  dg <- new.env()
  for (f in c("dag", "select", "ground_truth", "finalize", "checks", "pipeline"))
    sys.source(file.path(repo_dir, "01_data_generation", "R", paste0(f, ".R")), envir = dg)
  cal <- read_calibration(file.path(repo_dir, "01_data_generation", "config", "calibration.yaml"))
  c2 <- cfg
  c2$params$sizes[c("full_set", "reference_subsample", "astronaut_set", "ground_truth")] <- list(2000, 300, 300, 2000)
  c2$params$replicates$count <- 2
  data_dir <- file.path(tempdir(), "e2e_data"); out_dir <- file.path(tempdir(), "e2e_out")
  suppressMessages(dg$run_pipeline(c2, cal, data_dir))
  p2 <- modifyList(params, list(mgm = list(lambda_selection = "fixed"), stability = list(subsamples = 2),
                                orientations = list(max = 8)))
  out <- suppressMessages(run_shape(cfg, cal, p2, data_dir, out_dir))
  nv <- length(params$variants)
  expect_equal(nrow(out$per_replicate), 6 * nv)
  expect_equal(nrow(out$summary), 3 * nv)
  expect_setequal(unique(out$summary$variant), unlist(params$variants))
  expect_true(all(file.exists(file.path(out_dir, c(
    "summary.csv", "per_replicate.csv", "implication_rejection_rates.csv", "power_by_edge.csv",
    "edge_recovery.csv", "analysis_graph_edges.csv", "analysis_graph_cpdag.csv", "provenance.txt",
    "edge_stability_r001_full_set.csv")))))
  expect_equal(nrow(read.csv(file.path(out_dir, "implication_rejection_rates.csv"))), 571)
  expect_true(all(out$per_replicate$tests == 571))
})
