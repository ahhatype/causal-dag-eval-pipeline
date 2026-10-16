test_that("duration is truncated without point masses at the bounds", {
  s <- sim_source(D, 20000, 11)
  k <- cal$cont$cumulative_mission_duration
  raw <- s$cumulative_mission_duration * k$sd + k$mu
  d <- cfg$params$duration
  e60 <- raw[s$mission_era == 0]
  expect_true(all(e60 >= d$era_1960s$min - 1e-9 & e60 <= d$era_1960s$max + 1e-9))
  expect_true(all(raw[s$mission_era == 1] >= d$era_2000s$shuttle$min - 1e-9))
  bounds <- c(d$era_1960s$min, d$era_1960s$max, d$era_2000s$shuttle$min, d$era_2000s$shuttle$max,
              d$era_2000s$iss$min, d$era_2000s$iss$max)
  expect_equal(sum(sapply(bounds, function(b) sum(abs(raw - b) < 1e-6))), 0)
  expect_equal(median(e60), d$era_1960s$median, tolerance = 0.1)
})

test_that("drug uptake follows its DAG parents, not mission era directly", {
  for (id in names(cfg$params$uptake)) {
    f <- node_formula(id, cfg, cal)
    expect_match(f, "medical_prevention_capability", fixed = TRUE)
    expect_false(grepl("mission_era", f, fixed = TRUE))
  }
  expect_match(node_formula("k_citrate", cfg, cal), "individual_factors", fixed = TRUE)
  expect_match(node_formula("bisphosphonates", cfg, cal), "pre_flight_fitness", fixed = TRUE)
})

test_that("confounded drug uptake is on target within each era", {
  s <- sim_source(D, 40000, 21)
  for (id in c("k_citrate", "bisphosphonates")) {
    u <- cfg$params$uptake[[id]]
    expect_equal(mean(s[[id]][s$mission_era == 0]), u[1], tolerance = 0.15)
    expect_equal(mean(s[[id]][s$mission_era == 1]), u[2], tolerance = 0.1)
  }
  expect_gt(cor(s$k_citrate, s$individual_factors), 0.05)
  expect_lt(cor(s$bisphosphonates, s$pre_flight_fitness), -0.05)
})

test_that("astronaut-set: first retained full-set records, every record passes selection", {
  full <- sim_source(D, 10000, 123)
  a <- astronaut_from_full(full, cfg, cal, 900)
  expect_equal(nrow(a$data), 900)
  expect_true(all(retention_mask(a$data, cfg, cal)))
  expect_identical(a$data, full[retention_mask(full, cfg, cal)][1:900])
  expect_equal(a$retained, sum(retention_mask(full, cfg, cal)))
  expect_error(astronaut_from_full(full[1:100], cfg, cal, 900), "retained")
})

test_that("direct evaluation of the node formulas reproduces simcausal from the same noise", {
  s <- sim_source(D, 2000, 5)
  env <- new.env(parent = globalenv())
  for (id in noise_ids(cfg)) assign(paste0("U_", id), s[[paste0("U_", id)]], envir = env)
  forms <- compiled_formulas(cfg, cal)
  for (id in cfg$order) assign(id, eval(forms[[id]], envir = env), envir = env)
  for (id in cfg$order) expect_equal(rep_len(env[[id]], 2000), s[[id]], tolerance = 1e-12, label = id)
})

test_that("do_outcome_prob matches simcausal interventions", {
  forms <- compiled_formulas(cfg, cal)
  s <- sim_source(D, 3000, 9)
  r <- sim_actions(D, cfg, 3000, 9, list(hi = list(node = "hydration", value = 1)))
  envs <- new.env(parent = globalenv())
  for (id in noise_ids(cfg)) assign(paste0("U_", id), s[[paste0("U_", id)]], envir = envs)
  for (id in cfg$order) assign(id, eval(forms[[id]], envir = envs), envir = envs)
  expect_equal(do_outcome_prob(envs, cfg, cal, forms, "hydration", 1), outcome_prob(r$hi, cfg, cal), tolerance = 1e-12)
})

test_that("pop_spread and truth_grid: binary spread equals 2p(1 - p)|contrast|", {
  expect_equal(pop_spread(c(0.1, 0.3), c(0.8, 0.2)), 2 * 0.8 * 0.2 * 0.2)
  expect_equal(pop_spread(c(0.5, 0.5, 0.5), rep(1 / 3, 3)), 0)
  expect_equal(pop_spread(c(0.1, 0.3), c(1, 0)), 0)
  g <- truth_grid(c(0, 0, 0, 1), "binary", 20)
  expect_equal(g$w, c(0.75, 0.25))
  g <- truth_grid(seq(0, 1, length.out = 10001), "continuous", 20)
  expect_equal(g$x, (1:20 - 0.5) / 20, tolerance = 1e-9)
  expect_equal(sum(g$w), 1)
})

test_that("ground truth: zero for non-ancestors, all three types for both populations, known directions", {
  gt <- ground_truth(small_cfg(), cal)
  v <- gt$values
  expect_equal(nrow(v), 35 * 3 * 2)
  expect_setequal(unique(v$truth_type), c("per_unit", "pop", "rec"))
  expect_true(all(v$value[!v$ancestor_set] == 0))
  expect_true(all(v$value[v$ancestor_set & v$truth_type != "per_unit"] > 0))
  expect_true(all(v$mcse >= 0))
  pu <- v[v$truth_type == "per_unit" & v$population == "source", ]
  expect_gt(pu$value[pu$feature == "individual_factors"], 0)
  expect_lt(pu$value[pu$feature == "hydration"], 0)
  expect_true(gt$retained > 0 && gt$retained < gt$n)
  expect_setequal(unique(gt$binary_checks$feature), c("bisphosphonates", "k_citrate", "thiazides"))
})

test_that("unresolved_pairs flags adjacent per-unit truths within 2 MCSE", {
  tv <- data.frame(feature = rep(c("a", "b", "c"), 2), ancestor_set = TRUE, truth_type = "per_unit",
                   population = rep(c("source", "selected"), each = 3),
                   value = c(0.10, 0.099, 0.05, 0.10, 0.05, 0.01), mcse = 0.001)
  expect_equal(unresolved_pairs(tv, "source"), "a ~ b")
  expect_length(unresolved_pairs(tv, "selected"), 0)
})

test_that("stale calibration is detected when config changes", {
  tmp <- file.path(tempdir(), "cfg_copy"); dir.create(tmp, showWarnings = FALSE)
  file.copy(list.files(cfg_dir, full.names = TRUE), tmp, overwrite = TRUE)
  c2 <- read_config(tmp)
  expect_silent(check_calibration_current(c2, cal))
  cat("\n# comment only\n", file = file.path(tmp, "params.yaml"), append = TRUE)
  expect_silent(check_calibration_current(read_config(tmp), cal))
  cat("unused_param: 1\n", file = file.path(tmp, "params.yaml"), append = TRUE)
  expect_error(check_calibration_current(read_config(tmp), cal), "make recalibrate")
})

test_that("pipeline writes every replicate with distinct seeds, and checks catch a broken result", {
  out_dir <- file.path(tempdir(), "pipeline_out")
  c2 <- small_cfg()
  out <- suppressMessages(run_pipeline(c2, cal, out_dir))
  expect_true(all(file.exists(file.path(out_dir, c(
    "feature_sets.csv", "truth_values.csv", "truth_binary_checks.csv", "scm_check.csv",
    "data_summary.csv", "provenance.txt")))))
  for (r in 1:2) {
    rd <- replicate_dir(out_dir, r)
    expect_true(all(file.exists(file.path(rd, c("full_set.csv", "reference_subsample.csv", "astronaut_set.csv")))))
    expect_equal(ncol(data.table::fread(file.path(rd, "full_set.csv"))), 46)
  }
  expect_false(identical(data.table::fread(file.path(replicate_dir(out_dir, 1), "full_set.csv")),
                         data.table::fread(file.path(replicate_dir(out_dir, 2), "full_set.csv"))))
  expect_equal(nrow(out$summary), 6)
  expect_equal(replicate_seeds(c2$params, 1)$full_set, c2$params$seeds$full_set)
  ck <- out$checks
  expect_true(ck$pass[ck$check == "non_ancestors_zero_effect"])
  expect_true(ck$pass[ck$check == "replicate_count"] && ck$pass[ck$check == "dataset_sizes"])
  broken <- out$summary; broken$n[1] <- 10
  ck2 <- run_checks(broken, out$ground_truth$values, c2)
  expect_false(ck2$pass[ck2$check == "dataset_sizes"])
  expect_equal(ck2$severity[ck2$check == "dataset_sizes"], "error")
  expect_equal(ck2$severity[ck2$check == "ancestor_ranks_resolved"], "warning")
})

test_that("config tables are reproducible from the appendix (skipped without references/)", {
  appendix <- file.path(dirname(dg_dir), "references", "appendix-renal-dag-SA-07566.md")
  skip_if_not(file.exists(appendix), "references/ not present")
  skip_if(Sys.which("python3") == "", "python3 not available")
  tmp <- file.path(tempdir(), "tables"); dir.create(tmp, showWarnings = FALSE)
  status <- system2("python3", c(file.path(dg_dir, "tools", "build_config_tables.py"), appendix, tmp))
  expect_equal(status, 0)
  for (f in c("nodes.csv", "edges.csv"))
    expect_identical(unname(tools::md5sum(file.path(tmp, f))), unname(tools::md5sum(file.path(cfg_dir, f))))
})
