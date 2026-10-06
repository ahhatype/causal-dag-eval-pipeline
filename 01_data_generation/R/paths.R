# Reads one directory setting from .env (no other variables are loaded) and resolves it against the repo.
env_dir <- function(repo_dir, key, default) {
  value <- default
  env_file <- file.path(repo_dir, ".env")
  if (file.exists(env_file)) {
    line <- grep(sprintf("^\\s*%s\\s*=", key), readLines(env_file, warn = FALSE), value = TRUE)
    if (length(line)) value <- trimws(gsub("^['\"]|['\"]$", "", trimws(sub("^[^=]*=", "", line[1]))))
  }
  if (grepl("^/", value)) value else file.path(repo_dir, sub("^\\./", "", value))
}
