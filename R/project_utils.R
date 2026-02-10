# R/project_utils.R
# ------------------------------------------------------------
# Project utilities (logging, dependency checks, data I/O helpers)
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

# -------------------------
# Small data helpers
# -------------------------

safe_posix <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXct")) return(x)
  as.POSIXct(x, tz = tz)
}

median_impute <- function(dt, cols) {
  for (cc in cols) {
    if (!cc %in% names(dt)) next
    med <- suppressWarnings(stats::median(dt[[cc]], na.rm = TRUE))
    if (!is.finite(med)) med <- 0
    idx <- !is.finite(dt[[cc]]) | is.na(dt[[cc]])
    dt[[cc]][idx] <- med
  }
  invisible(dt)
}

read_mimic <- function(data_dir, rel_path, ...) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for read_mimic().", call. = FALSE)
  }
  f <- file.path(data_dir, rel_path)
  if (!file.exists(f)) {
    stop(
      sprintf(
        "Could not find '%s' under data_dir='%s'.\n\nExpected path: %s\n\n" ,
        rel_path, data_dir, normalizePath(f, winslash = "/", mustWork = FALSE)
      ),
      call. = FALSE
    )
  }
  data.table::fread(f, ...)
}
