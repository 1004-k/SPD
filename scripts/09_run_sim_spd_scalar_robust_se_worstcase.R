# scripts/09_run_sim_spd_scalar_robust_se_worstcase.R
# ------------------------------------------------------------
# Robust-SE sensitivity (worst-case) for the scalar SPD simulation.
#
# Compare model-based vs robust (sandwich) standard errors for the Cox coefficient
# in the strongest dependence setting: gamma_true=1.0, rho_meas=1.0.
#
# Scenarios:
# - N in N_VEC (default: 200, 800, 2000)
# - Two rate settings:
#   * baseline: LAMBDA_DEV_BASE (default 0.25), LAMBDA_EVENT (default 0.10)
#   * fastdev:  LAMBDA_DEV_FAST (default 0.50), LAMBDA_EVENT (same)
#
# Outputs (under output/<OUT_TAG>/):
# - raw/spd_scalar_robustSE_replicates_<OUT_TAG>.csv
# - spd_scalar_robustSE_summary_<OUT_TAG>.csv
# - figures/Fig_SPD_scalar_robustSE_coverage_<setting>.pdf
# ------------------------------------------------------------

source("R/project_utils.R")

cfg <- init_project(
  n_cores = as.integer(Sys.getenv("N_CORES", "3")),
  seed    = as.integer(Sys.getenv("SEED", "2026")),
  out_dir = Sys.getenv("OUT_DIR", "output")
)
cfg$require_pkgs(c("data.table","survival","future.apply","future"))

# ---------- knobs ----------
out_tag <- Sys.getenv("OUT_TAG", "sim_spd_scalar_robustSE_worstcase")

B     <- as.integer(Sys.getenv("B", "500"))
N_vec <- as.integer(strsplit(Sys.getenv("N_VEC", "200,800,2000"), ",", fixed = TRUE)[[1]])

gamma_true <- as.numeric(Sys.getenv("GAMMA_TRUE", "1.0"))
rho_meas   <- as.numeric(Sys.getenv("RHO_MEAS", "1.0"))

t_max <- as.numeric(Sys.getenv("T_MAX", "5"))

lambda_event_base <- as.numeric(Sys.getenv("LAMBDA_EVENT", "0.10"))
lambda_dev_base   <- as.numeric(Sys.getenv("LAMBDA_DEV_BASE", "0.25"))
lambda_dev_fast   <- as.numeric(Sys.getenv("LAMBDA_DEV_FAST", "0.50"))

settings <- data.table::data.table(
  setting = c("baseline","fastdev"),
  lambda_dev = c(lambda_dev_base, lambda_dev_fast),
  lambda_event = c(lambda_event_base, lambda_event_base)
)

# output dirs
out_root <- file.path(cfg$out_dir, out_tag)
raw_dir  <- file.path(out_root, "raw")
fig_dir  <- file.path(out_root, "figures")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(cfg$log_dir, sprintf("robustSE_worstcase_%s.log", out_tag))
if (file.exists(log_file)) file.remove(log_file)

cfg$log_line(log_file, sprintf("[%s] Start robust-SE worst-case sensitivity
", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cfg$log_line(log_file, sprintf("OUT_TAG=%s | B=%d | N_VEC=%s | gamma=%.2f | rho=%.2f
",
                               out_tag, B, paste(N_vec, collapse=","), gamma_true, rho_meas))
cfg$write_session_info(tag = paste0("sessionInfo_robustSE_", out_tag))

future::plan(future::multisession, workers = cfg$n_cores)

# ---------- DGP (same as scalar SPD simulation) ----------
simulate_one_dataset_scalar <- function(
  N,
  t_max = 5,
  lambda_event_base = 0.10,
  lambda_dev_base   = 0.25,
  theta_z = 0.50,
  gamma_true = 0,
  rho_meas = 1.0,
  seed = 1L
) {
  set.seed(seed)
  id <- seq_len(N)

  z_true <- rnorm(N)
  rho_meas <- max(min(rho_meas, 1), 0)
  z_obs <- rho_meas * z_true + sqrt(1 - rho_meas^2) * rnorm(N)

  rate_dev <- lambda_dev_base * exp(gamma_true * z_true)
  rate_evt <- lambda_event_base * exp(theta_z * z_true)

  t_dev <- rexp(N, rate = rate_dev)
  t_evt <- rexp(N, rate = rate_evt)

  time <- pmin(t_dev, t_evt, t_max)
  dev  <- as.integer(t_dev <= t_evt & t_dev <= t_max)

  data.table::data.table(id = id, tstop = time, dev = dev, z_obs = z_obs)
}

# ---------- Cox fit helper: model vs robust SE ----------
fit_spd_se <- function(dt, se_type = c("model","robust"), ties = "efron") {
  se_type <- match.arg(se_type)

  n_dev <- sum(dt$dev == 1, na.rm = TRUE)
  sx <- stats::sd(dt$z_obs, na.rm = TRUE)

  out <- list(ok = FALSE, n_dev = n_dev,
              beta_std = NA_real_, se_std = NA_real_,
              warn = "", err = "")

  if (!is.finite(sx) || sx <= 0) {
    out$err <- "sd_nonfinite_or_zero"
    return(out)
  }
  if (!is.finite(n_dev) || n_dev < 5) {
    out$err <- "too_few_events"
    return(out)
  }

  wtxt <- ""
  dat <- data.frame(tstop = dt$tstop, dev = dt$dev, z_obs = dt$z_obs, id = dt$id)

  fit <- tryCatch(
    withCallingHandlers({
      if (se_type == "model") {
        survival::coxph(survival::Surv(tstop, dev) ~ z_obs, data = dat, ties = ties)
      } else {
        survival::coxph(
          survival::Surv(tstop, dev) ~ z_obs + survival::cluster(id),
          data = dat, ties = ties
        )
      }
    }, warning = function(w) {
      wtxt <<- paste(wtxt, conditionMessage(w), sep = " | ")
      invokeRestart("muffleWarning")
    }),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    out$err <- conditionMessage(fit)
    out$warn <- wtxt
    return(out)
  }

  b <- as.numeric(stats::coef(fit)[1])
  se <- suppressWarnings(sqrt(as.numeric(vcov(fit)[1, 1])))

  if (!is.finite(b) || !is.finite(se)) {
    out$err <- "coef_or_se_nonfinite"
    out$warn <- wtxt
    return(out)
  }

  out$beta_std <- b * sx
  out$se_std   <- se * sx
  out$ok <- is.finite(out$beta_std) && is.finite(out$se_std)
  out$warn <- wtxt
  out
}

# ---------- run ----------
rep_rows <- list()

for (k in seq_len(nrow(settings))) {
  set_name <- settings$setting[k]
  lam_dev  <- settings$lambda_dev[k]
  lam_evt  <- settings$lambda_event[k]

  for (N in N_vec) {
    cat(sprintf("[Setting=%s] N=%d (B=%d)
", set_name, N, B))
    cfg$log_line(log_file, sprintf("[%s] setting=%s N=%d lam_dev=%.3f lam_evt=%.3f
",
                                   format(Sys.time(), "%H:%M:%S"), set_name, N, lam_dev, lam_evt))

    reps <- future.apply::future_lapply(
      seq_len(B),
      function(b) {
        dt <- simulate_one_dataset_scalar(
          N = N,
          t_max = t_max,
          lambda_event_base = lam_evt,
          lambda_dev_base   = lam_dev,
          theta_z = 0.50,
          gamma_true = gamma_true,
          rho_meas = rho_meas,
          seed = cfg$seed + 100000L * b + 1000L * N + 10L * k
        )

        m  <- fit_spd_se(dt, se_type = "model")
        rb <- fit_spd_se(dt, se_type = "robust")

        data.table::data.table(
          setting = set_name,
          N = N,
          gamma_true = gamma_true,
          rho_meas = rho_meas,
          rep = b,
          n_dev = m$n_dev,
          beta_std_model = m$beta_std,
          se_std_model   = m$se_std,
          ok_model       = m$ok,
          beta_std_robust = rb$beta_std,
          se_std_robust   = rb$se_std,
          ok_robust       = rb$ok
        )
      },
      future.seed = TRUE
    )

    rep_rows[[paste(set_name, N, sep = "_")]] <- data.table::rbindlist(reps, fill = TRUE)
  }
}

rep_dt <- data.table::rbindlist(rep_rows, fill = TRUE)

# ---------- summarize ----------
z975 <- stats::qnorm(0.975)

summ <- rep_dt[, {
  ok_m <- ok_model == TRUE & is.finite(beta_std_model) & is.finite(se_std_model)
  ok_r <- ok_robust == TRUE & is.finite(beta_std_robust) & is.finite(se_std_robust)

  bt <- gamma_true[1]  # with rho=1.0 the pseudo-true target is gamma_true

  cov_m <- if (sum(ok_m) > 0) mean((beta_std_model[ok_m] - z975*se_std_model[ok_m] <= bt) &
                                  (beta_std_model[ok_m] + z975*se_std_model[ok_m] >= bt)) else NA_real_
  cov_r <- if (sum(ok_r) > 0) mean((beta_std_robust[ok_r] - z975*se_std_robust[ok_r] <= bt) &
                                  (beta_std_robust[ok_r] + z975*se_std_robust[ok_r] >= bt)) else NA_real_

  list(
    B_total = .N,
    n_dev_med = stats::median(n_dev, na.rm = TRUE),

    B_ok_model = sum(ok_m),
    fail_rate_model = 1 - (sum(ok_m) / .N),
    bias_model = if (sum(ok_m) > 0) mean(beta_std_model[ok_m] - bt) else NA_real_,
    mean_se_model = if (sum(ok_m) > 0) mean(se_std_model[ok_m]) else NA_real_,
    coverage_model = cov_m,

    B_ok_robust = sum(ok_r),
    fail_rate_robust = 1 - (sum(ok_r) / .N),
    bias_robust = if (sum(ok_r) > 0) mean(beta_std_robust[ok_r] - bt) else NA_real_,
    mean_se_robust = if (sum(ok_r) > 0) mean(se_std_robust[ok_r]) else NA_real_,
    coverage_robust = cov_r
  )
}, by = .(setting, N, gamma_true, rho_meas)]

# ---------- save ----------
rep_file <- file.path(raw_dir, sprintf("spd_scalar_robustSE_replicates_%s.csv", out_tag))
sum_file <- file.path(out_root, sprintf("spd_scalar_robustSE_summary_%s.csv", out_tag))

data.table::fwrite(rep_dt, rep_file)
data.table::fwrite(summ, sum_file)

cfg$log_line(log_file, sprintf("[%s] Wrote: %s
", format(Sys.time(), "%H:%M:%S"), rep_file))
cfg$log_line(log_file, sprintf("[%s] Wrote: %s
", format(Sys.time(), "%H:%M:%S"), sum_file))

cat("
Saved outputs to:
")
cat(sprintf(" - %s
 - %s
", rep_file, sum_file))

# ---------- simple B/W coverage plots ----------
make_cov_plot <- function(ss, setting_name) {
  pdf_path <- file.path(fig_dir, sprintf("Fig_SPD_scalar_robustSE_coverage_%s.pdf", setting_name))
  grDevices::pdf(pdf_path, width = 7.2, height = 4.8, useDingbats = FALSE)

  op <- par(mar = c(4.2, 4.2, 1.2, 0.8), mgp = c(2.4, 0.8, 0))

  ss <- ss[order(N)]
  plot(ss$N, ss$coverage_model, type = "b", pch = 1, lty = 1, col = 1,
       ylim = c(0.85, 1.00),
       xlab = "Sample size (N)",
       ylab = "Empirical 95% coverage",
       main = sprintf("Worst-case coverage: model vs robust SE (%s)", setting_name))
  lines(ss$N, ss$coverage_robust, type = "b", pch = 2, lty = 2, col = 1)
  abline(h = 0.95, lty = 3, col = 1)
  legend("bottomright",
         legend = c("Model-based SE", "Robust SE (cluster id)", "Nominal 0.95"),
         lty = c(1,2,3), pch = c(1,2,NA), bty = "n", col = 1)

  par(op)
  dev.off()
  cat(sprintf("Wrote %s
", pdf_path))
}

for (set_name in unique(summ$setting)) {
  ss <- summ[setting == set_name]
  make_cov_plot(ss, set_name)
}
