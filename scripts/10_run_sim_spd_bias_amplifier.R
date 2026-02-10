# scripts/10_run_sim_spd_bias_amplifier.R
# ------------------------------------------------------------
# Illustrative simulation: SPD as an amplifier of treatment-effect bias.
#
# Narrative:
# In some settings, treatment causes a small increase in deviation risk
# (e.g., mild side effects). When deviation is prognosis-dependent (SPD > 0),
# selective deviation can amplify this imbalance and induce bias in naive
# per-protocol-style analyses.
#
# This script simulates an "amplifier" setting and summarizes how mean bias
# changes across SPD strengths.
# ------------------------------------------------------------

source("R/project_utils.R")

cfg <- init_project(
  n_cores = as.integer(Sys.getenv("N_CORES", "3")),
  seed    = as.integer(Sys.getenv("SEED", "2026")),
  out_dir = Sys.getenv("OUT_DIR", "output")
)
cfg$require_pkgs(c("data.table", "survival", "future.apply", "future"))

future::plan(future::multisession, workers = cfg$n_cores)

# -------------------------
# Knobs
# -------------------------
out_tag <- Sys.getenv("OUT_TAG", "sim_spd_bias_amplifier")
B       <- as.integer(Sys.getenv("B", "500"))
N       <- as.integer(Sys.getenv("N", "2000"))

gamma_vec <- as.numeric(strsplit(Sys.getenv("GAMMA_VEC", "0,0.25,0.5,0.75,1.0,1.25,1.5"), ",", fixed = TRUE)[[1]])

# Data-generating parameters (illustration)
beta_true <- as.numeric(Sys.getenv("BETA_TRUE", "-0.5"))  # true log-HR for treatment (event)
rho_meas  <- as.numeric(Sys.getenv("RHO_MEAS",  "0.7"))   # measurement reliability of prognostic score
alpha_dev <- as.numeric(Sys.getenv("ALPHA_DEV", "0.3"))   # treatment-related deviation log-HR
theta_z   <- as.numeric(Sys.getenv("THETA_Z",   "0.5"))   # prognosis effect on event hazard

out_root <- file.path(cfg$out_dir, out_tag)
raw_dir  <- file.path(out_root, "raw")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(cfg$log_dir, sprintf("spd_bias_progress_%s.log", out_tag))
if (file.exists(log_file)) file.remove(log_file)

cfg$log_line(
  log_file,
  sprintf("[%s] Start amplifier simulation (alpha_dev=%.2f, rho_meas=%.2f)
",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"), alpha_dev, rho_meas)
)

# -------------------------
# DGP: amplifier model
# -------------------------
simulate_amplifier_dataset <- function(N, gamma_true, beta_true, rho_meas, alpha_dev, theta_z, seed) {
  set.seed(as.integer(seed))

  # Latent prognosis and observed score
  z_true <- rnorm(N)
  rho_meas <- max(min(rho_meas, 1), 0)
  z_obs  <- rho_meas * z_true + sqrt(1 - rho_meas^2) * rnorm(N)

  # Random treatment
  A <- rbinom(N, 1, 0.5)

  # Event time (outcome)
  lambda_event <- 0.10
  rate_evt <- lambda_event * exp(theta_z * z_true + beta_true * A)
  t_evt <- rexp(N, rate = rate_evt)

  # Deviation time (depends on prognosis and treatment)
  lambda_dev <- 0.25
  rate_dev <- lambda_dev * exp(gamma_true * z_true + alpha_dev * A)
  t_dev <- rexp(N, rate = rate_dev)

  t_max <- 5
  time <- pmin(t_evt, t_dev, t_max)

  status_evt <- as.integer(t_evt <= t_dev & t_evt <= t_max)
  status_dev <- as.integer(t_dev <= t_evt & t_dev <= t_max)

  data.table::data.table(
    id = seq_len(N),
    time = time,
    status_evt = status_evt,
    status_dev = status_dev,
    z_obs = z_obs,
    A = A
  )
}

run_one_rep <- function(b, gamma_true) {
  dt <- simulate_amplifier_dataset(
    N = N,
    gamma_true = gamma_true,
    beta_true = beta_true,
    rho_meas = rho_meas,
    alpha_dev = alpha_dev,
    theta_z = theta_z,
    seed = cfg$seed + 10000L * b + as.integer(round(100 * gamma_true))
  )

  # Estimated SPD (researcher-observed)
  fit_spd <- survival::coxph(survival::Surv(time, status_dev) ~ z_obs, data = dt)
  spd_est <- unname(stats::coef(fit_spd)[["z_obs"]])

  # Naive event model (adjusts for z_obs but ignores selection induced by deviation)
  fit_evt <- survival::coxph(survival::Surv(time, status_evt) ~ A + z_obs, data = dt)
  beta_est <- unname(stats::coef(fit_evt)[["A"]])

  data.table::data.table(
    rep = b,
    gamma_true = gamma_true,
    spd_est = spd_est,
    beta_est = beta_est,
    bias = beta_est - beta_true
  )
}

cfg$log_line(log_file, "Running replicates...
")

results_list <- list()

for (g in gamma_vec) {
  cat(sprintf("Processing gamma_true = %.2f ...
", g))

  reps <- future.apply::future_lapply(
    seq_len(B),
    function(b) run_one_rep(b, g),
    future.seed = TRUE
  )
  results_list[[as.character(g)]] <- data.table::rbindlist(reps)
}

full_dt <- data.table::rbindlist(results_list)

summ_dt <- full_dt[, .(
  mean_spd = mean(spd_est, na.rm = TRUE),
  mean_bias = mean(bias, na.rm = TRUE),
  se_bias = stats::sd(bias, na.rm = TRUE) / sqrt(.N)
), by = gamma_true]

# -------------------------
# Save
# -------------------------
rep_file <- file.path(raw_dir, sprintf("spd_bias_replicates_%s.csv", out_tag))
sum_file <- file.path(out_root, sprintf("spd_bias_summary_%s.csv", out_tag))

data.table::fwrite(full_dt, rep_file)
data.table::fwrite(summ_dt, sum_file)

cfg$log_line(log_file, sprintf("Done. Saved to %s and %s
", rep_file, sum_file))
cat(sprintf("
Simulation complete.
Saved outputs to:
 - %s
 - %s
", rep_file, sum_file))
