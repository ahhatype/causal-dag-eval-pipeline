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

test_that("drug uptake follows its DAG parent, not mission era directly", {
  for (id in names(cfg$params$uptake)) {
    f <- node_formula(id, cfg, cal)
    expect_match(f, "medical_prevention_capability", fixed = TRUE)
    expect_false(grepl("mission_era", f, fixed = TRUE))
  }
})

test_that("astronaut-set: right size, every record passes selection, deterministic", {
  a <- draw_astronaut_set(D, cfg, cal, 900, 123)
  expect_equal(nrow(a$data), 900)
  expect_true(all(retention_mask(a$data, cfg, cal)))
  expect_identical(a$data, draw_astronaut_set(D, cfg, cal, 900, 123)$data)
})

test_that("ground truth: zero for non-ancestors, SEs present, known direction", {
  gt <- ground_truth(D, small_cfg(), cal)
  expect_equal(nrow(gt), 35)
  expect_true(all(gt$effect_full[!gt$ancestor_set] == 0))
  expect_true(all(gt$effect_astronaut[!gt$ancestor_set] == 0))
  expect_true(all(gt$se_full >= 0 & gt$se_astronaut >= 0))
  expect_gt(gt$effect_full[gt$feature == "individual_factors"], 0)
  expect_lt(gt$effect_full[gt$feature == "hydration"], 0)
})

test_that("unresolved_pairs flags adjacent effects within 2 SE", {
  gt <- data.frame(feature = c("a", "b", "c"), ancestor_set = TRUE,
                   effect_full = c(0.10, 0.099, 0.05), se_full = 0.001,
                   effect_astronaut = c(0.10, 0.05, 0.01), se_astronaut = 0.001)
  expect_equal(unresolved_pairs(gt, "full"), "a ~ b")
  expect_length(unresolved_pairs(gt, "astronaut"), 0)
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

test_that("pipeline runs end to end at reduced size and checks catch a broken result", {
  out_dir <- file.path(tempdir(), "pipeline_out")
  out <- suppressMessages(run_pipeline(small_cfg(), cal, out_dir))
  expect_true(all(file.exists(file.path(out_dir, c(
    "full_set.csv", "reference_subsample.csv", "astronaut_set.csv", "feature_sets.csv",
    "ground_truth_total_effects_full.csv", "ground_truth_total_effects_astronaut.csv",
    "data_summary.csv", "provenance.txt")))))
  expect_equal(ncol(data.table::fread(file.path(out_dir, "full_set.csv"))), 46)
  ck <- out$checks
  expect_true(ck$pass[ck$check == "non_ancestors_zero_effect"])
  broken <- out$results; broken$full <- broken$full[1:10]
  ck2 <- run_checks(broken, small_cfg())
  expect_false(ck2$pass[ck2$check == "full_set_dims"])
  expect_equal(ck2$severity[ck2$check == "full_set_dims"], "error")
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
