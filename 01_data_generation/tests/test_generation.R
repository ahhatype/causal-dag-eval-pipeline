test_that("config reproduces the adapted DAG (54 nodes, 90 edges, 16/19 feature split)", {
  expect_silent(validate_dag(cfg_dir, cfg))
})

test_that("every simulated node has a generation rule and the order is topological", {
  expect_length(cfg$order, 53)
  expect_no_error(sapply(cfg$order, node_rule, cfg = cfg))
  pos <- setNames(seq_along(cfg$order), cfg$order)
  expect_true(all(pos[cfg$edges$parent] < pos[cfg$edges$child]))
})

test_that("same seed gives identical data; different seed differs", {
  a <- sim_source(D, 2000, 1); b <- sim_source(D, 2000, 1); c2 <- sim_source(D, 2000, 2)
  expect_identical(a, b)
  expect_false(identical(a$nephrolithiasis, c2$nephrolithiasis))
})

test_that("common random numbers: intervening leaves non-descendants untouched and non-ancestors have zero effect", {
  base <- sim_source(D, 5000, 7)
  r <- sim_actions(D, cfg, 5000, 7, list(lo = list(node = "ultrasound", value = 0),
                                         hi = list(node = "ultrasound", value = 1)))
  expect_identical(r$lo$nephrolithiasis, r$hi$nephrolithiasis)
  expect_identical(r$lo$individual_factors, base$individual_factors)
  r2 <- sim_actions(D, cfg, 5000, 7, list(lo = list(node = "individual_factors", value = -1),
                                          hi = list(node = "individual_factors", value = 1)))
  expect_gt(mean(r2$hi$nephrolithiasis) - mean(r2$lo$nephrolithiasis), 0)
})

test_that("calibration matches the current config", {
  expect_silent(check_calibration_current(cfg, cal))
})

test_that("calibrated prevalence and selection retention are on target", {
  s <- sim_source(D, 20000, 3)
  expect_equal(mean(s$nephrolithiasis), 0.10, tolerance = 0.1)
  expect_equal(mean(retention_mask(s, cfg, cal)), 0.30, tolerance = 0.1)
  expect_equal(mean(s$mission_era), 0.75, tolerance = 0.03)
})

test_that("latent nodes are absent from analysis columns and 35 features are observed", {
  obs <- observed_ids(cfg)
  expect_false(any(cfg$nodes$id[!cfg$nodes$observed] %in% obs))
  expect_length(feature_ids(cfg, "all"), 35)
  expect_length(feature_ids(cfg, "ancestor"), 16)
  expect_true(all(feature_ids(cfg, "all") %in% obs))
})
