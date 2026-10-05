write_dataset <- function(df, name, cfg, out_dir) {
  obs <- observed_ids(cfg)
  latent <- cfg$nodes$id[cfg$nodes$simulate & !cfg$nodes$observed]
  data.table::fwrite(df[, ..obs], file.path(out_dir, paste0(name, ".csv")))
  data.table::fwrite(df[, c("ID", latent), with = FALSE], file.path(out_dir, paste0(name, "_latent.csv")))
  invisible(file.path(out_dir, paste0(name, ".csv")))
}

write_feature_sets <- function(cfg, out_dir) {
  n <- cfg$nodes
  f <- n[n$feature_set %in% c("ancestor", "non_ancestor"), c("id", "name", "scale", "role", "feature_set")]
  f$in_ancestor_set <- f$feature_set == "ancestor"
  f$in_all_features_set <- TRUE
  data.table::fwrite(f[, c("id", "name", "scale", "role", "in_ancestor_set", "in_all_features_set")],
                     file.path(out_dir, "feature_sets.csv"))
}
