# R/project_utils.R
# ------------------------------------------------------------
# Project utilities (output folders, logging, dependency checks)
# ------------------------------------------------------------

init_project <- function(n_cores = 1L,
                         seed = 2026L,
                         out_dir = "output",
                         logs_subdir = "logs") {
  n_cores <- suppressWarnings(as.integer(n_cores))
  if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L

  seed <- suppressWarnings(as.integer(seed))
  if (!is.finite(seed)) seed <- 2026L

  out_dir <- as.character(out_dir)
  log_dir <- file.path(out_dir, logs_subdir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

  require_pkgs <- function(pkgs) {
    stopifnot(is.character(pkgs))
    miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(miss) > 0) {
      stop(
        sprintf(
          "Missing required R packages: %s\nInstall them first, e.g.: install.packages(c(%s))",
          paste(miss, collapse = ", "),
          paste(sprintf("'%s'", miss), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    invisible(TRUE)
  }

  log_line <- function(path, line) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    cat(line, file = path, append = TRUE)
    invisible(TRUE)
  }

  write_session_info <- function(tag = "sessionInfo") {
    p <- file.path(out_dir, sprintf("%s.txt", tag))
    capture.output(sessionInfo(), file = p)
    invisible(p)
  }

  list(
    n_cores = n_cores,
    seed = seed,
    out_dir = out_dir,
    log_dir = log_dir,
    require_pkgs = require_pkgs,
    log_line = log_line,
    write_session_info = write_session_info
  )
}
