read_provenance <- function(data_dir) {
  f <- file.path(data_dir, "provenance.txt")
  if (!file.exists(f)) stop("no provenance.txt in ", data_dir, ": run `make data`", call. = FALSE)
  l <- readLines(f)
  kv <- regmatches(l, regexpr(": ", l), invert = TRUE)
  setNames(vapply(kv, `[`, "", 2), vapply(kv, `[`, "", 1))
}

check_data_current <- function(data_dir, cal) {
  prov <- read_provenance(data_dir)
  if (!identical(unname(prov["config hash"]), cal$config_hash))
    stop("data in ", data_dir, " were generated from a different config: run `make data`", call. = FALSE)
  as.integer(prov["replicates"])
}

check_analysis_graph <- function(ag, imp) {
  ok <- c(nodes_36 = length(ag$nodes) == 36, edges_63 = nrow(ag$edges) == 63,
          implications_567 = nrow(imp) == 567, selection_flag_1 = sum(imp$expected_under_selection) == 1)
  if (!all(ok)) stop("analysis graph check failed: ", paste(names(ok)[!ok], collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

mc_se <- function(x) if (sum(!is.na(x)) > 1) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))) else NA_real_

shape_manifest <- function(params, cal, ag, sv_dir) {
  f <- tempfile()
  on.exit(unlink(f))
  code <- list.files(file.path(sv_dir, "R"), pattern = "\\.R$", full.names = TRUE)
  writeLines(c(yaml::as.yaml(params), cal$config_hash, paste(ag$nodes, collapse = ","),
               unlist(lapply(sort(code), function(x) paste(basename(x), unname(tools::md5sum(x)))))), f)
  unname(tools::md5sum(f))
}

# Each (replicate, data set) is one unit, saved when complete. `resume` keeps units written under the same
# settings, config and code; `fresh` discards them; with neither, existing units are an error.
prepare_units <- function(unit_dir, manifest, resume, fresh) {
  mf <- file.path(unit_dir, "manifest.txt")
  existing <- list.files(unit_dir, pattern = "\\.rds$", full.names = TRUE)
  if (fresh) unlink(existing)
  else if (length(existing)) {
    if (!resume) stop(unit_dir, " already holds results: pass --resume to continue them or --fresh to discard them", call. = FALSE)
    if (!file.exists(mf) || !identical(readLines(mf, n = 1), manifest))
      stop(unit_dir, " was written with different settings, config or code: pass --fresh to discard it", call. = FALSE)
  }
  dir.create(unit_dir, showWarnings = FALSE, recursive = TRUE)
  writeLines(manifest, mf)
}

save_unit <- function(x, file) {
  tmp <- paste0(file, ".tmp")
  saveRDS(x, tmp)
  file.rename(tmp, file)
}

run_unit <- function(r, ds, n_rep, first, data_dir, out_dir, params, ag, imp, del, cpdag) {
  d <- data.table::fread(file.path(data_dir, "replicates", sprintf("r%03d", r), paste0(ds, ".csv")))
  seed <- params$discovery_seed + r

  tests <- test_implications(d, imp, ag, params$fdr)
  pw <- edge_power(d, del, ag, params$power$alpha)
  md <- mgm_data(d, ag)
  fit <- run_mgm_pc(md, params, seed)
  primary <- names(params$variants)[1]

  data.table::fwrite(tests, file.path(out_dir, "detail", sprintf("implications_r%03d_%s.csv", r, ds)))
  ic <- implication_calibration(tests, params$calibration$alpha)
  if (first && params$stability$subsamples > 0) {
    message(sprintf("[replicate %d] %s: stability over %d subsamples", r, ds, params$stability$subsamples))
    data.table::fwrite(stability_edges(md, fit$graphs[[primary]], params),
                       file.path(out_dir, sprintf("edge_stability_r%03d_%s.csv", r, ds)))
  }

  rj <- tests$rejected
  base <- data.frame(replicate = r, dataset = ds, n = nrow(d), tests = sum(!is.na(tests$p)),
                     rejected = sum(rj),
                     selection_violation_detected = sum(rj & tests$expected_under_selection),
                     rejected_nonadditive = sum(rj & tests$nonadditive_response & !tests$expected_under_selection),
                     rejected_other = sum(rj & !tests$nonadditive_response & !tests$expected_under_selection),
                     edges_detectable = sum(pw$detected),
                     calibration_ks_rejected = ic$calibration_ks_rejected, raw_rejection_rate = ic$raw_rejection_rate,
                     discovery_dropped = paste(attr(md, "dropped"), collapse = "; "),
                     lambda_cc = fit$lambda[1], lambda_cd = fit$lambda[2], lambda_dd = fit$lambda[3])
  rows <- list(); edges <- list()
  for (v in names(params$variants)) {
    learned <- learned_edges(fit$graphs[[v]])
    cmp <- compare_learned(ag$nodes, ag$edges, learned, params, imp[, c("X", "Y", "Z")])
    data.table::fwrite(learned, file.path(out_dir, "detail", sprintf("learned_cpdag_r%03d_%s_%s.csv", r, ds, v)))
    rows[[length(rows) + 1]] <- cbind(base, variant = params$variants[[v]],
                                      structure_metrics(cpdag, learned), cmp$summary)
    edges[[length(edges) + 1]] <- data.frame(dataset = ds, variant = params$variants[[v]], replicate = r,
                                             key = edge_key(learned$from, learned$to),
                                             mark = edge_mark(learned$from, learned$to, learned$type))
  }
  list(rows = rows, edges = edges, calib = cbind(replicate = r, dataset = ds, ic),
       rej = data.frame(dataset = ds, implication = seq_len(nrow(tests)), rejected = rj),
       pow = cbind(dataset = ds, pw))
}

run_shape <- function(cfg, cal, params, data_dir, out_dir, replicates = NULL, repo_dir = NULL,
                      resume = FALSE, fresh = FALSE, sv_dir = file.path(repo_dir, "03_shape_validation")) {
  available <- check_data_current(data_dir, cal)
  reps <- seq_len(min(available, if (is.null(replicates)) available else replicates))
  dir.create(file.path(out_dir, "detail"), showWarnings = FALSE, recursive = TRUE)

  ag <- analysis_graph(cfg)
  imp <- flag_implications(shipley_basis(ag$nodes, ag$edges), ag, params)
  check_analysis_graph(ag, imp)
  cpdag <- true_cpdag(ag$nodes, ag$edges)
  del <- deletion_implications(ag)
  data.table::fwrite(ag$edges, file.path(out_dir, "analysis_graph_edges.csv"))
  data.table::fwrite(cpdag, file.path(out_dir, "analysis_graph_cpdag.csv"))

  unit_dir <- file.path(out_dir, "units")
  prepare_units(unit_dir, shape_manifest(params, cal, ag, sv_dir), resume, fresh)

  rows <- list(); rej <- list(); pow <- list(); edges <- list(); calib <- list()
  for (r in reps) {
    for (ds in params$datasets) {
      f <- file.path(unit_dir, sprintf("r%03d_%s.rds", r, ds))
      if (file.exists(f)) {
        message(sprintf("[replicate %d/%d] %s (kept from an earlier run)", r, length(reps), ds))
        u <- readRDS(f)
      } else {
        message(sprintf("[replicate %d/%d] %s", r, length(reps), ds))
        u <- run_unit(r, ds, length(reps), r == reps[1], data_dir, out_dir, params, ag, imp, del, cpdag)
        save_unit(u, f)
      }
      rows <- c(rows, u$rows); edges <- c(edges, u$edges)
      calib[[length(calib) + 1]] <- u$calib; rej[[length(rej) + 1]] <- u$rej; pow[[length(pow) + 1]] <- u$pow
    }
  }

  per_rep <- do.call(rbind, rows)
  data.table::fwrite(per_rep, file.path(out_dir, "per_replicate.csv"))
  summ <- summarize_replicates(per_rep)
  data.table::fwrite(summ, file.path(out_dir, "summary.csv"))

  rj <- do.call(rbind, rej)
  rates <- aggregate(rejected ~ dataset + implication, rj, mean)
  names(rates)[3] <- "rejection_rate"
  wide <- reshape(rates, idvar = "implication", timevar = "dataset", direction = "wide")
  data.table::fwrite(cbind(imp, wide[order(wide$implication), -1, drop = FALSE]),
                     file.path(out_dir, "implication_rejection_rates.csv"))

  write_calibration(do.call(rbind, calib), out_dir)

  pw <- do.call(rbind, pow)
  power <- merge(setNames(aggregate(detected ~ dataset + edge, pw, mean), c("dataset", "edge", "detection_rate")),
                 setNames(aggregate(p ~ dataset + edge, pw, median), c("dataset", "edge", "median_p")))
  data.table::fwrite(power, file.path(out_dir, "power_by_edge.csv"))

  recovery <- edge_recovery(do.call(rbind, edges), cpdag, length(reps))
  data.table::fwrite(recovery, file.path(out_dir, "edge_recovery.csv"))
  data.table::fwrite(difference_dag_table(recovery, cpdag, params$variants[[1]], params$difference_dag$spurious_threshold),
                     file.path(out_dir, "difference_dag_edges.csv"))
  write_shape_provenance(out_dir, data_dir, cal, length(reps), repo_dir,
                         extra = sprintf("difference DAG spurious-edge display threshold: %s (primary variant %s)",
                                         params$difference_dag$spurious_threshold, params$variants[[1]]))
  list(per_replicate = per_rep, summary = summ, power = power)
}

summarize_replicates <- function(per_rep) {
  metrics <- setdiff(names(per_rep)[vapply(per_rep, is.numeric, logical(1))], c("replicate", "n"))
  out <- lapply(split(per_rep, list(per_rep$dataset, per_rep$variant), drop = TRUE), function(x) {
    data.frame(dataset = x$dataset[1], variant = x$variant[1], replicates = nrow(x), n = x$n[1],
               selection_violation_detection_rate = mean(x$selection_violation_detected > 0),
               any_unexpected_rejection_rate = mean(x$rejected_other > 0),
               comparison_available_rate = mean(!is.na(x$truth_held_by_learned)),
               t(setNames(c(vapply(metrics, function(m) mean(x[[m]], na.rm = TRUE), 0), vapply(metrics, function(m) mc_se(x[[m]]), 0)),
                          c(paste0(metrics, "_mean"), paste0(metrics, "_se")))), check.names = FALSE)
  })
  do.call(rbind, out)
}

# Per true edge, data set and variant, what the learner did in each replicate, as rates that sum to 1:
#   correct      learned mark equals the true CPDAG mark
#   reversed     true edge directed, learned directed the other way
#   unoriented   true edge directed, learned undirected
#   overoriented true edge undirected, learned directed either way
#   absent       no adjacency
# cpdag_mark_rate repeats correct_rate. False adjacencies keep only adjacency_rate (the rate they were added).
recovery_rates <- function(marks, true_mark, n_rep) {
  directed <- true_mark != "undir"
  rate <- function(x) sum(x) / n_rep
  list(correct_rate = rate(marks == true_mark),
       reversed_rate = if (directed) rate(marks %in% c("fwd", "rev") & marks != true_mark) else 0,
       unoriented_rate = if (directed) rate(marks == "undir") else 0,
       overoriented_rate = if (directed) 0 else rate(marks %in% c("fwd", "rev")),
       absent_rate = 1 - length(marks) / n_rep)
}

edge_recovery <- function(learned, cpdag, n_rep) {
  tk <- edge_key(cpdag$from, cpdag$to)
  tm <- setNames(edge_mark(cpdag$from, cpdag$to, cpdag$type), tk)
  cats <- c("correct_rate", "reversed_rate", "unoriented_rate", "overoriented_rate", "absent_rate")
  out <- lapply(split(learned, list(learned$dataset, learned$variant), drop = TRUE), function(x) {
    keys <- union(tk, x$key)
    per <- lapply(keys, function(k) {
      marks <- x$mark[x$key == k]
      if (k %in% tk) recovery_rates(marks, tm[[k]], n_rep) else setNames(as.list(rep(NA_real_, 5)), cats)
    })
    cbind(data.frame(dataset = x$dataset[1], variant = x$variant[1], edge = keys, true_edge = keys %in% tk,
                     true_mark = ifelse(keys %in% tk, tm[keys], NA),
                     adjacency_rate = vapply(keys, function(k) sum(x$key == k) / n_rep, 0)),
          do.call(rbind.data.frame, per)[, cats], cpdag_mark_rate = vapply(per, `[[`, 0, "correct_rate"))
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res[order(res$dataset, res$variant, !res$true_edge, -res$adjacency_rate), ]
}

# Figure-ready table for the primary variant: one row per true edge, and per false adjacency added in at least
# `threshold` of the replicates. Direction is the true one for directed edges, alphabetical otherwise.
difference_dag_table <- function(recovery, cpdag, variant, threshold) {
  rec <- recovery[recovery$variant == variant, ]
  tk <- edge_key(cpdag$from, cpdag$to)
  cats <- c("correct", "reversed", "unoriented", "overoriented", "absent")
  rates <- paste0(cats, "_rate")
  parts <- strsplit(rec$edge, " ", fixed = TRUE)
  parent <- vapply(parts, `[`, "", 1)
  child <- vapply(parts, `[`, "", 2)
  true_dir <- match(rec$edge, tk)
  directed <- !is.na(true_dir) & cpdag$type[true_dir] == "dir"
  parent[directed] <- cpdag$from[true_dir[directed]]
  child[directed] <- cpdag$to[true_dir[directed]]
  modal <- rep("added", nrow(rec))
  modal[rec$true_edge] <- cats[apply(as.matrix(rec[rec$true_edge, rates]), 1, which.max)]
  out <- data.frame(dataset = rec$dataset, edge = rec$edge, parent = parent, child = child, true_edge = rec$true_edge,
                    true_mark = rec$true_mark, rec[, rates], added_rate = ifelse(rec$true_edge, NA, rec$adjacency_rate),
                    modal_category = modal, row.names = NULL)
  out[out$true_edge | out$added_rate >= threshold, ]
}

write_shape_provenance <- function(out_dir, data_dir, cal, n_rep, repo_dir, file = "provenance.txt", extra = NULL) {
  commit <- if (!is.null(repo_dir)) tryCatch(system2("git", c("-C", repo_dir, "rev-parse", "HEAD"), stdout = TRUE,
                                                     stderr = FALSE)[1], error = function(e) NA_character_) else NA
  pkgs <- c("rCausalMGM", "dagitty")
  writeLines(c(paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
               paste("git commit:", commit),
               paste("data config hash:", cal$config_hash),
               paste("replicates analysed:", n_rep),
               extra,
               paste("R:", R.version.string),
               sprintf("%s: %s (built under R %s)", pkgs,
                       vapply(pkgs, function(x) as.character(utils::packageVersion(x)), ""),
                       vapply(pkgs, function(x) utils::packageDescription(x)$Built, ""))),
             file.path(out_dir, file))
}

write_calibration <- function(per_rep, out_dir) {
  data.table::fwrite(per_rep, file.path(out_dir, "implication_calibration.csv"))
  summ <- calibration_summary(per_rep)
  data.table::fwrite(summ, file.path(out_dir, "implication_calibration_summary.csv"))
  summ
}

run_calibration <- function(cfg, cal, params, data_dir, out_dir, replicates = NULL, repo_dir = NULL) {
  available <- check_data_current(data_dir, cal)
  reps <- seq_len(min(available, if (is.null(replicates)) available else replicates))
  dir.create(file.path(out_dir, "detail"), showWarnings = FALSE, recursive = TRUE)
  ag <- analysis_graph(cfg)
  imp <- flag_implications(shipley_basis(ag$nodes, ag$edges), ag, params)
  check_analysis_graph(ag, imp)
  calib <- list()
  for (r in reps) for (ds in params$datasets) {
    message(sprintf("[replicate %d/%d] %s", r, length(reps), ds))
    d <- data.table::fread(file.path(data_dir, "replicates", sprintf("r%03d", r), paste0(ds, ".csv")))
    tests <- test_implications(d, imp, ag, params$fdr)
    data.table::fwrite(tests, file.path(out_dir, "detail", sprintf("implications_r%03d_%s.csv", r, ds)))
    calib[[length(calib) + 1]] <- cbind(replicate = r, dataset = ds, implication_calibration(tests, params$calibration$alpha))
  }
  write_shape_provenance(out_dir, data_dir, cal, length(reps), repo_dir, "implication_calibration_provenance.txt")
  write_calibration(do.call(rbind, calib), out_dir)
}
