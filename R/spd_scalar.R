# R/spd_scalar.R
# ------------------------------------------------------------
# Scalar standardized prognostic dependence (SPD)
#
# SPD is defined as the Cox regression coefficient for deviation per 1 SD
# higher baseline prognostic score. This file provides:
# - estimate_spd_scalar(): estimate SPD and (optionally) a small invariance check
# - compute_reference_spd(): pseudo-true target via large-N Monte Carlo
# ------------------------------------------------------------

estimate_spd_scalar <- function(long_dt,
                                score_col = "z_obs",
                                time_col  = "tstop",
                                event_col = "dev",
                                ties = "efron",
                                do_invariance = TRUE) {
  stopifnot(is.data.frame(long_dt))
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.", call. = FALSE)
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.", call. = FALSE)
  }

  dt <- data.table::as.data.table(long_dt)

  if (!all(c(score_col, time_col, event_col) %in% names(dt))) {
    stop("Required columns missing in long_dt.", call. = FALSE)
  }

  x <- dt[[score_col]]
  t <- dt[[time_col]]
  d <- dt[[event_col]]

  n_dev <- suppressWarnings(sum(d == 1, na.rm = TRUE))

  fit_one <- function(x_vec) {
    out <- list(
      beta_raw = NA_real_, se_raw = NA_real_,
      beta_std = NA_real_, se_std = NA_real_,
      ok = FALSE, warn = "", err = ""
    )

    sx <- stats::sd(x_vec, na.rm = TRUE)
    if (!is.finite(sx) || sx <= 0) {
      out$err <- "sd_nonfinite_or_zero"
      return(out)
    }
    if (!is.finite(n_dev) || n_dev < 5) {
      out$err <- "too_few_events"
      return(out)
    }

    wtxt <- ""
    fit <- tryCatch(
      withCallingHandlers(
        survival::coxph(survival::Surv(t, d) ~ x_vec, ties = ties),
        warning = function(w) {
          wtxt <<- paste(wtxt, conditionMessage(w), sep = " | ")
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      out$err <- conditionMessage(fit)
      out$warn <- wtxt
      return(out)
    }

    b <- as.numeric(stats::coef(fit)[1])
    se <- suppressWarnings(sqrt(as.numeric(stats::vcov(fit)[1, 1])))

    if (!is.finite(b) || !is.finite(se)) {
      out$err <- "coef_or_se_nonfinite"
      out$warn <- wtxt
      return(out)
    }

    out$beta_raw <- b
    out$se_raw   <- se
    out$beta_std <- b * sx
    out$se_std   <- se * sx
    out$ok <- is.finite(out$beta_std) && is.finite(out$se_std)
    out$warn <- wtxt
    out
  }

  main <- fit_one(x)

  inv_dt <- data.table::data.table()
  if (isTRUE(do_invariance)) {
    # Moderate affine transforms (avoid extreme scaling that can induce separation)
    affine_list <- list(
      `1` = function(u) u,
      `2` = function(u) 2 * u + 1,
      `3` = function(u) 0.5 * u - 1
    )

    inv_dt <- data.table::rbindlist(
      lapply(names(affine_list), function(k) {
        xx <- affine_list[[k]](x)
        res <- fit_one(xx)
        data.table::data.table(
          scale_id = as.integer(k),
          beta_std = res$beta_std,
          se_std   = res$se_std,
          ok       = res$ok,
          warn     = res$warn,
          err      = res$err
        )
      }),
      fill = TRUE
    )
  }

  list(
    n_dev = n_dev,
    beta_raw = main$beta_raw,
    se_raw   = main$se_raw,
    ok_raw   = is.finite(main$beta_raw) && is.finite(main$se_raw),
    beta_std = main$beta_std,
    se_std   = main$se_std,
    ok_std   = main$ok,
    warn_std = main$warn,
    err_std  = main$err,
    inv_dt   = inv_dt
  )
}

compute_reference_spd <- function(sim_fun,
                                  N_ref,
                                  B_ref,
                                  seed,
                                  ...) {
  betas <- rep(NA_real_, B_ref)
  ok    <- rep(FALSE, B_ref)

  for (b in seq_len(B_ref)) {
    dat <- sim_fun(N = N_ref, seed = seed + 10000L * b, ...)
    est <- estimate_spd_scalar(dat$long_dt, score_col = "z_obs", do_invariance = FALSE)
    betas[b] <- est$beta_std
    ok[b] <- isTRUE(est$ok_std) && is.finite(betas[b])
  }

  n_ok <- sum(ok)
  beta_ref <- if (n_ok > 0) mean(betas[ok]) else NA_real_
  list(beta_ref = beta_ref, n_ok = n_ok)
}
