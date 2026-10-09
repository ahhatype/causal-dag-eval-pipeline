# Node formulas are strings shared by calibration (plain R) and the simcausal DAG.

slug_id <- function(x) {
  x <- tolower(x)
  x <- gsub("k\\+", "k", x)
  gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", x))
}

num <- function(x) sprintf("%.15g", x)

read_config <- function(cfg_dir) {
  cfg <- list(
    dir    = cfg_dir,
    params = yaml::read_yaml(file.path(cfg_dir, "params.yaml")),
    nodes  = read.csv(file.path(cfg_dir, "nodes.csv"), stringsAsFactors = FALSE),
    edges  = rbind(read.csv(file.path(cfg_dir, "edges.csv"), stringsAsFactors = FALSE),
                   read.csv(file.path(cfg_dir, "simulation_edges.csv"), stringsAsFactors = FALSE))
  )
  sim_ids <- cfg$nodes$id[cfg$nodes$simulate]
  cfg$edges <- cfg$edges[cfg$edges$parent %in% sim_ids & cfg$edges$child %in% sim_ids, ]
  cfg$order <- topo_order(sim_ids, cfg$edges)
  cfg
}

topo_order <- function(ids, edges) {
  indeg <- setNames(rep(0L, length(ids)), ids)
  for (ch in edges$child) indeg[ch] <- indeg[ch] + 1L
  out <- character(0)
  ready <- ids[indeg[ids] == 0]
  while (length(ready)) {
    n <- ready[1]; ready <- ready[-1]; out <- c(out, n)
    for (ch in edges$child[edges$parent == n]) {
      indeg[ch] <- indeg[ch] - 1L
      if (indeg[ch] == 0L) ready <- c(ready, ch)
    }
  }
  if (length(out) != length(ids)) stop("DAG is cyclic")
  out
}

parents_of <- function(id, cfg) cfg$edges[cfg$edges$child == id, ]

lin_expr <- function(id, cfg) {
  e <- parents_of(id, cfg)
  e <- e[!is.na(e$coef), ]
  if (!nrow(e)) stop("no numeric coefficients into ", id)
  paste(sprintf("(%s) * %s", num(e$coef), e$parent), collapse = " + ")
}

coef_of <- function(parent, child, cfg) {
  e <- cfg$edges[cfg$edges$parent == parent & cfg$edges$child == child, ]
  stopifnot(nrow(e) == 1, !is.na(e$coef))
  e$coef
}

node_rule <- function(id, cfg) {
  p <- cfg$params
  special <- c(mission_era = "era_binary", cumulative_mission_duration = "duration",
               ultrasound = "ultrasound", detect_mrm_stone = "detect_mrm",
               detect_long_term_health_outcomes = "detect_ltho",
               loss_of_mission_objectives = "lmo", loss_of_mission = "lm")
  if (id %in% names(special)) return(special[[id]])
  if (id %in% names(p$era_values)) return("era_fixed")
  if (id %in% names(p$uptake)) {
    e <- parents_of(id, cfg)
    return(if (any(!is.na(e$coef))) "uptake_logit" else "era_prob")
  }
  if (nrow(parents_of(id, cfg)) == 0) return("root_normal")
  scale <- cfg$nodes$scale[cfg$nodes$id == id]
  if (scale == "continuous") return("cont")
  if (id %in% names(p$base_rates)) return("logit")
  stop("no generation rule for node: ", id)
}

duration_raw_formula <- function(cfg) {
  d <- cfg$params$duration; a <- d$era_1960s; s <- d$era_2000s
  trunc_norm <- function(mean, sd, lo, hi, U) {
    pa <- pnorm((lo - mean) / sd); pb <- pnorm((hi - mean) / sd)
    sprintf("(%s + %s * qnorm(%s + %s * %s))", num(mean), num(sd), num(pa), num(pb - pa), U)
  }
  U <- "U_cumulative_mission_duration"
  sprintf("ifelse(mission_era == 0, exp(%s), ifelse(U_duration_mix < %s, %s, %s))",
          trunc_norm(log(a$median), a$log_sd, log(a$min), log(a$max), U),
          num(s$shuttle_share),
          trunc_norm(s$shuttle$mean, s$shuttle$sd, s$shuttle$min, s$shuttle$max, U),
          trunc_norm(s$iss$mean, s$iss$sd, s$iss$min, s$iss$max, U))
}

node_formula <- function(id, cfg, cal) {
  p <- cfg$params
  U <- paste0("U_", id)
  switch(node_rule(id, cfg),
    era_binary = sprintf("as.numeric(%s < %s)", U, num(p$era$share_2000s)),
    root_normal = sprintf("qnorm(%s)", U),
    era_fixed = {
      v <- p$era_values[[id]]
      if (v[1] == v[2]) "0" else {
        sh <- p$era$share_2000s
        mu <- (1 - sh) * v[1] + sh * v[2]
        sdv <- sqrt(sh * (1 - sh)) * abs(v[2] - v[1])
        sprintf("((%s + (%s) * mission_era) - %s) / %s", num(v[1]), num(v[2] - v[1]), num(mu), num(sdv))
      }
    },
    duration = sprintf("((%s) - %s) / %s", duration_raw_formula(cfg),
                       num(cal$cont[[id]]$mu), num(cal$cont[[id]]$sd)),
    cont = {
      k <- cal$cont[[id]]
      sprintf("((%s + %s * qnorm(%s)) - %s) / %s", lin_expr(id, cfg), num(k$rs), U, num(k$mu), num(k$sd))
    },
    logit = sprintf("as.numeric(%s < plogis(%s + %s))", U, num(cal$logit[[id]]), lin_expr(id, cfg)),
    era_prob = {
      u <- p$uptake[[id]]
      sprintf("as.numeric(%s < (%s + (%s) * (medical_prevention_capability > 0)))", U, num(u[1]), num(u[2] - u[1]))
    },
    uptake_logit = {
      a <- cal$uptake[[id]]
      sprintf("as.numeric(%s < plogis(ifelse(medical_prevention_capability > 0, %s, %s) + %s))",
              U, num(a[2]), num(a[1]), lin_expr(id, cfg))
    },
    ultrasound = sprintf("as.numeric(%s < ifelse(medical_monitoring_capability > 0, %s, %s))",
                         U, num(1 - p$flip), num(p$flip)),
    detect_mrm = {
      d <- p$detection$mrm_stone; sl <- num(d$mrm_slope)
      sens <- sprintf("plogis(qlogis(ifelse(ultrasound == 1, %s, %s)) + %s * mineralized_renal_material)",
                      num(d$with_ultrasound$sensitivity), num(d$without_ultrasound$sensitivity), sl)
      fpr <- sprintf("ifelse(ultrasound == 1, %s, %s)",
                     num(1 - d$with_ultrasound$specificity), num(1 - d$without_ultrasound$specificity))
      sprintf("as.numeric(%s < ifelse(pmax(nephrolithiasis, as.numeric(mineralized_renal_material > %s)) == 1, %s, %s))",
              U, num(d$mrm_threshold), sens, fpr)
    },
    detect_ltho = {
      d <- p$detection$ltho
      sens <- sprintf("ifelse(surveillance > 0, %s, %s)", num(d$with_surveillance$sensitivity), num(d$without_surveillance$sensitivity))
      fpr <- sprintf("ifelse(surveillance > 0, %s, %s)", num(1 - d$with_surveillance$specificity), num(1 - d$without_surveillance$specificity))
      sprintf("as.numeric(%s < ifelse(long_term_health_outcomes > %s, %s, %s))", U, num(d$threshold), sens, fpr)
    },
    lmo = sprintf("as.numeric(%s < ifelse(evacuation == 1, %s, plogis(%s + (%s) * task_performance)))",
                  U, num(1 - p$flip), num(cal$logit[[id]]), num(coef_of("task_performance", id, cfg))),
    lm = sprintf("as.numeric(%s < ifelse(loss_of_mission_objectives == 1, %s, %s))", U, num(1 - p$flip), num(p$flip))
  )
}

noise_ids <- function(cfg) c(cfg$order, "duration_mix", "selection")

feature_ids <- function(cfg, set = c("all", "ancestor")) {
  set <- match.arg(set)
  n <- cfg$nodes
  n$id[n$feature_set %in% if (set == "ancestor") "ancestor" else c("ancestor", "non_ancestor")]
}

outcome_id <- "nephrolithiasis"

observed_ids <- function(cfg) {
  n <- cfg$nodes
  intersect(cfg$order, n$id[n$observed & n$simulate])
}
