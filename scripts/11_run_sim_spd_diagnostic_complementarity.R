# scripts/11_run_sim_spd_diagnostic_complementarity.R
# ------------------------------------------------------------
# Diagnostic-complementarity simulation (Table 3, eTable S5, Figure 2).
#
# Aim: show what SPD adds to common deviation/censoring diagnostics and what
# it cannot establish. This is a compact diagnostic-complementarity analysis,
# not a replacement for the main finite-sample operating-characteristics study.
#
# The reported Table 3 was computed from the replicate results stored in
# results/complementarity (see scripts/16). A re-run with this script uses a
# different random-number stream, so its results agree with Table 3 within
# Monte Carlo error rather than digit for digit.
#
# Scenarios:
# A. Positive versus negative prognostic sorting
# B. Prognostic sorting plus treatment-related deviation
# C. Strong latent selection with a noisy observed score
# D. Near-zero baseline SPD with treatment-related deviation
#
# Required packages: data.table, survival
# Environment variables: N (default 2000), B (default 500), SEED, OUT_DIR
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

N <- as.integer(Sys.getenv("N", "2000"))
B <- as.integer(Sys.getenv("B", "500"))
SEED <- as.integer(Sys.getenv("SEED", "20260918"))
OUT_DIR <- Sys.getenv("OUT_DIR", "output/complementarity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cox_one <- function(time, event, x) {
  ok <- is.finite(time) & is.finite(event) & is.finite(x)
  if (sum(event[ok] == 1) < 2L || sd(x[ok]) <= 0) {
    return(list(beta = NA_real_, se = NA_real_, ok = FALSE))
  }
  fit <- tryCatch(
    suppressWarnings(coxph(Surv(time[ok], event[ok]) ~ x[ok], ties = "breslow")),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(beta = NA_real_, se = NA_real_, ok = FALSE))
  b <- as.numeric(coef(fit)[1])
  se <- sqrt(as.numeric(vcov(fit)[1, 1]))
  list(beta = b, se = se, ok = is.finite(b) && is.finite(se))
}

smd_one <- function(a, b) {
  if (length(a) < 2L || length(b) < 2L) return(NA_real_)
  den <- sqrt((var(a) + var(b)) / 2)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  (mean(a) - mean(b)) / den
}

km_before <- function(time, event) {
  o <- order(time, method = "radix")
  tt <- time[o]
  ee <- event[o]
  ans <- numeric(length(tt))
  surv <- 1
  i <- 1L
  while (i <= length(tt)) {
    j <- i
    while (j < length(tt) && tt[j + 1L] == tt[i]) j <- j + 1L
    ans[i:j] <- surv
    d <- sum(ee[i:j] == 1L, na.rm = TRUE)
    nrisk <- length(tt) - i + 1L
    if (d > 0 && nrisk > 0) surv <- surv * max(0, 1 - d / nrisk)
    i <- j + 1L
  }
  out <- numeric(length(time))
  out[o] <- ans
  out
}

simulate_one <- function(n, gamma, alpha, rho, seed,
                         t_max = 5, beta_true = -0.5,
                         theta_z = 0.8, lambda_dev = 0.15,
                         lambda_event = 0.04) {
  set.seed(seed)
  z <- rnorm(n)
  score <- rho * z + sqrt(max(0, 1 - rho^2)) * rnorm(n)
  x1 <- z + rnorm(n, sd = 0.75)
  x2 <- 0.5 * z + rnorm(n, sd = 1.0)
  treatment <- rbinom(n, 1, 0.5)
  t_dev <- rexp(n, rate = lambda_dev * exp(gamma * z + alpha * treatment))
  t_event <- rexp(n, rate = lambda_event * exp(theta_z * z + beta_true * treatment))
  time <- pmin(t_dev, t_event, t_max)
  dev <- as.integer(t_dev <= t_event & t_dev <= t_max)
  outcome <- as.integer(t_event <= t_dev & t_event <= t_max)
  data.table(time = time, dev = dev, outcome = outcome, score = score,
             x1 = x1, x2 = x2, treatment = treatment,
             beta_true = beta_true)
}

diagnose_one <- function(dat) {
  fit_spd <- cox_one(dat$time, dat$dev, dat$score)
  score_sd <- sd(dat$score)
  spd <- if (fit_spd$ok) fit_spd$beta * score_sd else NA_real_
  spd_se <- if (fit_spd$ok) fit_spd$se * score_sd else NA_real_

  # Both SMDs compare observed deviation against no observed deviation.
  # X1/X2 are noisy proxies even when score = Z.
  d1 <- dat$dev == 1L
  rate <- mean(dat$dev)
  cov_smd <- max(abs(smd_one(dat$x1[d1], dat$x1[!d1])),
                 abs(smd_one(dat$x2[d1], dat$x2[!d1])))
  score_smd <- smd_one(dat$score[d1], dat$score[!d1])

  # Descriptive marginal exit-time weights for ALL subjects.
  # rESS_percent is rESSexit in the manuscript, not conditional IPCW ESS.
  g_before <- km_before(dat$time, dat$dev)
  w <- 1 / pmax(g_before, 1e-8)
  ress <- sum(w)^2 / sum(w^2)
  k <- max(1L, ceiling(0.01 * nrow(dat)))
  top1 <- sum(sort(w, decreasing = TRUE)[seq_len(k)]) / sum(w)

  # Gap from the conditional DGM coefficient includes noncollapsibility
  # and omitted prognosis; it is not a causal-bias estimate.
  fit_out <- cox_one(dat$time, dat$outcome, dat$treatment)
  naive_loghr <- if (fit_out$ok) fit_out$beta else NA_real_
  coefficient_gap <- if (fit_out$ok) fit_out$beta - unique(dat$beta_true) else NA_real_

  data.table(n = nrow(dat), n_dev = sum(dat$dev),
             deviation_rate = rate, spd = spd, spd_se = spd_se,
             spd_hr = exp(spd), abs_covariate_smd = cov_smd,
             prognostic_score_smd = score_smd,
             rESS = ress, rESS_percent = 100 * ress / nrow(dat),
             top1_weight_share = top1, naive_outcome_loghr = naive_loghr,
             naive_coefficient_gap = coefficient_gap, spd_ok = fit_spd$ok,
             outcome_ok = fit_out$ok)
}

scenario <- data.table(
  scenario = c("A_positive_selection", "A_negative_selection",
               "B_high_spd_treatment_related", "C_low_observed_spd",
               "D_near_zero_baseline_spd"),
  gamma = c(1.0, -1.0, 1.5, 1.5, 0.0),
  alpha = c(0.0, 0.0, 1.5, 1.5, 2.0),
  rho = c(1.0, 1.0, 1.0, 0.3, 1.0),
  interpretation = c(
    "High-risk patients deviate more; no treatment-related deviation",
    "Low-risk patients deviate more; no treatment-related deviation",
    "High-risk and treatment-related deviation",
    "Strong latent selection but noisy observed score",
    "No baseline prognostic sorting; treatment-related deviation"
  )
)

raw_list <- vector("list", nrow(scenario) * B)
z <- 0L
for (s in seq_len(nrow(scenario))) {
  for (b in seq_len(B)) {
    z <- z + 1L
    dat <- simulate_one(N, scenario$gamma[s], scenario$alpha[s],
                        scenario$rho[s], SEED + 100000L * (s - 1L) + b)
    row <- diagnose_one(dat)
    row[, `:=`(scenario = scenario$scenario[s], gamma = scenario$gamma[s],
               alpha = scenario$alpha[s], rho = scenario$rho[s], rep = b,
               interpretation = scenario$interpretation[s])]
    raw_list[[z]] <- row
  }
}
raw <- rbindlist(raw_list, fill = TRUE)

metric_cols <- c("deviation_rate", "spd", "spd_hr", "abs_covariate_smd",
                 "prognostic_score_smd", "rESS_percent", "top1_weight_share",
                 "naive_outcome_loghr", "naive_coefficient_gap")
summary <- raw[, {
  ans <- list(
    n_rep = .N,
    spd_failure_rate = 1 - mean(spd_ok),
    outcome_failure_rate = 1 - mean(outcome_ok)
  )
  for (m in metric_cols) {
    v <- get(m)[is.finite(get(m))]
    ans[[paste0(m, "_mean")]] <- if (length(v)) mean(v) else NA_real_
    ans[[paste0(m, "_sd")]] <- if (length(v) > 1L) sd(v) else NA_real_
    ans[[paste0(m, "_mcse")]] <- if (length(v) > 1L) sd(v) / sqrt(length(v)) else NA_real_
  }
  ans
}, by = .(scenario, gamma, alpha, rho)]

fwrite(raw, file.path(OUT_DIR, "diagnostic_complementarity_raw.csv"))
fwrite(summary, file.path(OUT_DIR, "diagnostic_complementarity_summary.csv"))
cat("Saved diagnostic-complementarity simulation outputs to:", OUT_DIR, "\n")
print(summary[, .(scenario, deviation_rate_mean, spd_mean, spd_hr_mean,
                 abs_covariate_smd_mean, prognostic_score_smd_mean,
                 rESS_percent_mean, top1_weight_share_mean,
                 naive_outcome_loghr_mean, naive_coefficient_gap_mean)])
