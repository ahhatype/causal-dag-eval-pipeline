# Each node draws from its own uniform noise node (U_<id>), so an intervention replaces
# a formula without shifting any random draws (common random numbers across actions).

suppressMessages(library(simcausal))

VECFUN <- c("qnorm", "plogis", "qlogis", "ifelse", "pmin", "pmax", "exp", "sqrt")

# simcausal does not allow underscores in node names.
sc_name <- function(x) gsub("_", ".", x, fixed = TRUE)
sc_formula <- function(formula, ids) {
  for (id in ids[order(-nchar(ids))]) formula <- gsub(paste0("\\b", id, "\\b"), sc_name(id), formula)
  formula
}
unsc_names <- function(df) { data.table::setnames(df, names(df), gsub(".", "_", names(df), fixed = TRUE)); df }

dag_node <- function(name, formula, ids) {
  eval(substitute(node(NM, distr = "rconst", const = EXPR),
                  list(NM = sc_name(name), EXPR = parse(text = sc_formula(formula, ids))[[1]])))
}

build_dag <- function(cfg, cal) {
  D <- NULL
  invisible(capture.output(D <- build_dag_quiet(cfg, cal)))
  D
}

build_dag_quiet <- function(cfg, cal) {
  D <- DAG.empty()
  ids <- c(cfg$order, paste0("U_", noise_ids(cfg)))
  for (id in noise_ids(cfg)) {
    D <- D + eval(substitute(node(NM, distr = "runif", min = 0, max = 1), list(NM = sc_name(paste0("U_", id)))))
  }
  for (id in cfg$order) D <- D + dag_node(id, node_formula(id, cfg, cal), ids)
  suppressMessages(suppressWarnings(set.DAG(D, vecfun = VECFUN, verbose = FALSE)))
}

sim_source <- function(D, n, seed) {
  df <- NULL
  invisible(capture.output(df <- sim(D, n = n, rndseed = seed, verbose = FALSE)))
  unsc_names(data.table::as.data.table(df))
}

sim_actions <- function(D, cfg, n, seed, values) {
  res <- NULL
  invisible(capture.output({
    for (nm in names(values)) {
      D <- D + action(nm, nodes = eval(substitute(node(NM, distr = "rconst", const = V),
                                                  list(NM = sc_name(values[[nm]]$node), V = values[[nm]]$value))))
    }
    res <- sim(D, actions = names(values), n = n, rndseed = seed, verbose = FALSE)
  }))
  lapply(res, function(d) unsc_names(data.table::as.data.table(d)))
}
