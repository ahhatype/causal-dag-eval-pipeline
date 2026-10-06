calibrate <- function(cfg, out_file = NULL) {
  p <- cfg$params
  n <- p$sizes$calibration
  set.seed(p$seeds$calibration)
  env <- new.env()
  for (id in noise_ids(cfg)) assign(paste0("U_", id), runif(n), envir = env)
  cal <- list(cont = list(), logit = list())
  ev <- function(f) eval(parse(text = f), envir = env)

  for (id in cfg$order) {
    rule <- node_rule(id, cfg)
    if (rule == "cont") {
      lin <- ev(lin_expr(id, cfg))
      rs <- if (!is.null(p$resid_sd[[id]])) p$resid_sd[[id]] else sqrt(max(1 - var(lin), p$resid_sd_floor^2))
      raw <- lin + rs * qnorm(get(paste0("U_", id), envir = env))
      cal$cont[[id]] <- list(rs = rs, mu = mean(raw), sd = sd(raw))
    } else if (rule == "duration") {
      raw <- ev(duration_raw_formula(cfg))
      cal$cont[[id]] <- list(mu = mean(raw), sd = sd(raw))
    } else if (rule == "logit") {
      lin <- ev(lin_expr(id, cfg))
      cal$logit[[id]] <- uniroot(function(b) mean(plogis(b + lin)) - p$base_rates[[id]],
                                 c(-30, 30), tol = 1e-10)$root
    } else if (rule == "lmo") {
      evac <- ev("evacuation") == 0
      lin <- coef_of("task_performance", id, cfg) * ev("task_performance")[evac]
      cal$logit[[id]] <- uniroot(function(b) mean(plogis(b + lin)) - p$base_rates$loss_of_mission_objectives_no_evac,
                                 c(-30, 30), tol = 1e-10)$root
    }
    assign(id, ev(node_formula(id, cfg, cal)), envir = env)
  }
  s <- p$selection
  lin_sel <- s$b_individual_factors * env$individual_factors + s$b_fitness * env$pre_flight_fitness
  cal$selection_alpha <- uniroot(function(a) mean(plogis(a + lin_sel)) - s$target_retention,
                                 c(-30, 30), tol = 1e-10)$root
  cal$config_hash <- config_hash(cfg$dir)
  if (!is.null(out_file)) yaml::write_yaml(cal, out_file, precision = 15)
  cal
}

read_calibration <- function(file) yaml::read_yaml(file)

config_hash <- function(cfg_dir) {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(c(yaml::as.yaml(yaml::read_yaml(file.path(cfg_dir, "params.yaml"))),
               readLines(file.path(cfg_dir, "nodes.csv")),
               readLines(file.path(cfg_dir, "edges.csv"))), f)
  unname(tools::md5sum(f))
}

check_calibration_current <- function(cfg, cal) {
  if (!identical(cal$config_hash, config_hash(cfg$dir)))
    stop("config changed since calibration: run `make recalibrate`", call. = FALSE)
  invisible(TRUE)
}
