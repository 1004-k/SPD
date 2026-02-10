# scripts/07_run_sim_spd_scalar.R
# ------------------------------------------------------------
# Simulation study for scalar standardized prognostic dependence (SPD).
#
# Focus:
# - SPD = Cox coefficient for deviation per 1 SD higher baseline prognostic score
# - Unit-free invariance to affine score rescaling
# - Operating characteristics of the SPD estimator (bias/coverage vs a pseudo-true reference)
#
# Outputs (under output/<OUT_TAG>/):
# - raw/spd_scalar_replicates_<OUT_TAG>.csv
# - raw/spd_scalar_invariance_<OUT_TAG>.csv
# - raw/spd_scalar_invariance_long_<OUT_TAG>.csv
# - spd_scalar_summary_<OUT_TAG>.csv
# - spd_scalar_reference_<OUT_TAG>.csv
# ------------------------------------------------------------

source("R/project_utils.R")
source("R/spd_scalar.R")

cfg <- init_project(
  n_cores = as.integer(Sys.getenv("N_CORES", "3")),
  seed    = as.integer(Sys.getenv("SEED", "2026")),
  out_dir = Sys.getenv("OUT_DIR", "output")
)
cfg$require_pkgs(c("data.table","survival","future.apply","future"))

# -------------------------
# Data-generating mechanism
# -------------------------
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

  # Latent prognosis
  z_true <- rnorm(N)

  # Observed score with correlation rho_meas to z_true
  rho_meas <- max(min(rho_meas, 1), 0)
  z_obs <- rho_meas * z_true + sqrt(1 - rho_meas^2) * rnorm(N)

  # Deviation time (cause-specific hazard depends on latent prognosis)
  rate_dev <- lambda_dev_base * exp(gamma_true * z_true)

  # Primary event time (competes with deviation; acts as censoring in deviation model)
  rate_evt <- lambda_event_base * exp(theta_z * z_true)

  t_dev <- rexp(N, rate = rate_dev)
  t_evt <- rexp(N, rate = rate_evt)

  time <- pmin(t_dev, t_evt, t_max)
  dev  <- as.integer(t_dev <= t_evt & t_dev <= t_max)

  long_dt <- data.table::data.table(
    id = id,
    tstop = time,
    dev = dev,
    z_obs = z_obs
  )
  list(long_dt = long_dt)
}

# -------------------------
# Knobs (overridable via env vars)
# -------------------------
B      <- as.integer(Sys.getenv("B", "200"))
t_max  <- as.numeric(Sys.getenv("T_MAX", "5"))

N_vec     <- as.integer(strsplit(Sys.getenv("N_VEC", "200,800,2000"), ",", fixed = TRUE)[[1]])
gamma_vec <- as.numeric(strsplit(Sys.getenv("GAMMA_VEC", "0,0.25,0.5,1.0"), ",", fixed = TRUE)[[1]])
rho_vec   <- as.numeric(strsplit(Sys.getenv("RHO_VEC", "1.0,0.7,0.4"), ",", fixed = TRUE)[[1]])

lambda_dev_base   <- as.numeric(Sys.getenv("LAMBDA_DEV", "0.25"))
lambda_event_base <- as.numeric(Sys.getenv("LAMBDA_EVENT", "0.10"))

do_ref <- as.integer(Sys.getenv("DO_REF", "1")) == 1L
N_ref  <- as.integer(Sys.getenv("N_REF", "50000"))
B_ref  <- as.integer(Sys.getenv("B_REF", "3"))

out_tag  <- Sys.getenv("OUT_TAG", sprintf("sim_spd_scalar_B%d", B))
out_root <- file.path(cfg$out_dir, out_tag)
raw_dir  <- file.path(out_root, "raw")
fig_dir  <- file.path(out_root, "figures")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(cfg$log_dir, sprintf("spd_scalar_progress_%s.log", out_tag))
if (file.exists(log_file)) file.remove(log_file)

cfg$log_line(log_file, sprintf("[%s] Start scalar SPD simulation
", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cfg$log_line(log_file, sprintf("B=%d | N_vec=%s | gamma_vec=%s | rho_vec=%s
",
                               B, paste(N_vec, collapse=","), paste(gamma_vec, collapse=","), paste(rho_vec, collapse=",")))
cfg$write_session_info(tag = paste0("sessionInfo_spd_scalar_", out_tag))

# Use multicore/multisession for the replicate loop
future::plan(future::multisession, workers = cfg$n_cores)

# Scenario grid
grid <- data.table::CJ(N = N_vec, gamma_true = gamma_vec, rho_meas = rho_vec, unique = TRUE)
grid[, scen_id := sprintf("N%d_g%.2f_r%.1f", N, gamma_true, rho_meas)]

# Reference (pseudo-true) SPD per (gamma_true, rho_meas) (not by N)
ref_dt <- data.table::CJ(gamma_true = gamma_vec, rho_meas = rho_vec, unique = TRUE)
ref_dt[, beta_ref := NA_real_]
ref_dt[, n_ok := 0L]

if (do_ref) {
  cfg$log_line(log_file, sprintf("[%s] Computing reference SPD with N_ref=%d (B_ref=%d)
",
                                 format(Sys.time(), "%H:%M:%S"), N_ref, B_ref))
  for (j in seq_len(nrow(ref_dt))) {
    g <- ref_dt$gamma_true[j]
    r <- ref_dt$rho_meas[j]
    ref <- compute_reference_spd(
      sim_fun = function(N, seed, ...) simulate_one_dataset_scalar(
        N = N,
        seed = seed,
        t_max = t_max,
        lambda_event_base = lambda_event_base,
        lambda_dev_base   = lambda_dev_base,
        theta_z = 0.50,
        gamma_true = g,
        rho_meas = r
      ),
      N_ref = N_ref,
      B_ref = B_ref,
      seed = 1000 + j * 10
    )
    ref_dt$beta_ref[j] <- ref$beta_ref
    ref_dt$n_ok[j] <- ref$n_ok
    cfg$log_line(log_file, sprintf("  ref: gamma_true=%.2f rho=%.1f -> beta_ref=%.4f (ok=%d)
",
                                   g, r, ref$beta_ref, ref$n_ok))
  }
} else {
  cfg$log_line(log_file, sprintf("[%s] Skipping reference SPD computation (DO_REF=0)
", format(Sys.time(), "%H:%M:%S")))
}

# -------------------------
# Replicate runner
# -------------------------
one_rep <- function(N, gamma_true, rho_meas, b, seed_base = 2026L) {
  dat <- simulate_one_dataset_scalar(
    N = N,
    t_max = t_max,
    lambda_event_base = lambda_event_base,
    lambda_dev_base   = lambda_dev_base,
    theta_z = 0.50,
    gamma_true = gamma_true,
    rho_meas = rho_meas,
    seed = seed_base + 100000L * b + 1000L * N + as.integer(round(100 * gamma_true)) + as.integer(round(10 * rho_meas))
  )

  est <- estimate_spd_scalar(dat$long_dt, score_col = "z_obs", do_invariance = TRUE)

  inv <- est$inv_dt
  inv[, N := N]
  inv[, gamma_true := gamma_true]
  inv[, rho_meas := rho_meas]
  inv[, rep := b]

  list(
    main = data.table::data.table(
      N = N,
      gamma_true = gamma_true,
      rho_meas = rho_meas,
      rep = b,
      n_dev = est$n_dev,
      beta_raw = est$beta_raw,
      se_raw   = est$se_raw,
      ok_raw   = est$ok_raw,
      beta_std = est$beta_std,
      se_std   = est$se_std,
      ok_std   = est$ok_std,
      warn_std = est$warn,
      err_std  = est$err
    ),
    inv = inv
  )
}

rep_list_main <- vector("list", nrow(grid))
rep_list_inv  <- vector("list", nrow(grid))

cfg$log_line(log_file, sprintf("[%s] Running scenarios...
", format(Sys.time(), "%H:%M:%S")))

for (s in seq_len(nrow(grid))) {
  N <- grid$N[s]
  g <- grid$gamma_true[s]
  r <- grid$rho_meas[s]
  scen <- grid$scen_id[s]

  cat(sprintf("[Scenario %d/%d] %s
", s, nrow(grid), scen))
  cfg$log_line(log_file, sprintf("[%s] Scenario %d/%d: %s
", format(Sys.time(), "%H:%M:%S"), s, nrow(grid), scen))

  reps <- future.apply::future_lapply(
    seq_len(B),
    function(b) one_rep(N, g, r, b, seed_base = cfg$seed),
    future.seed = TRUE
  )
  rep_list_main[[s]] <- data.table::rbindlist(lapply(reps, `[[`, "main"), fill = TRUE)
  rep_list_inv[[s]]  <- data.table::rbindlist(lapply(reps, `[[`, "inv"),  fill = TRUE)
}

rep_dt <- data.table::rbindlist(rep_list_main, fill = TRUE)
inv_dt <- data.table::rbindlist(rep_list_inv,  fill = TRUE)

# Merge reference SPD and define target
rep_dt <- merge(rep_dt, ref_dt[, .(gamma_true, rho_meas, beta_ref)], by = c("gamma_true","rho_meas"), all.x = TRUE)

# If reference is unavailable, fall back to the classic attenuation approximation
rep_dt[, beta_target := ifelse(is.finite(beta_ref), beta_ref, gamma_true * rho_meas)]

# Performance summaries
z975 <- stats::qnorm(0.975)
summ <- rep_dt[, {
  ok <- ok_std == TRUE & is.finite(beta_std) & is.finite(se_std)
  B_total <- .N
  B_ok <- sum(ok)
  bt <- beta_target[1]
  list(
    B_total = B_total,
    B_ok = B_ok,
    fail_rate = 1 - (B_ok / B_total),
    n_dev_med = stats::median(n_dev, na.rm = TRUE),
    beta_target = bt,
    bias = if (B_ok > 0) mean(beta_std[ok] - bt) else NA_real_,
    rmse = if (B_ok > 0) sqrt(mean((beta_std[ok] - bt)^2)) else NA_real_,
    mean_se = if (B_ok > 0) mean(se_std[ok]) else NA_real_,
    coverage = if (B_ok > 0) mean((beta_std[ok] - z975*se_std[ok] <= bt) &
                                  (beta_std[ok] + z975*se_std[ok] >= bt)) else NA_real_
  )
}, by = .(N, gamma_true, rho_meas)]

# -------------------------
# Invariance summaries (affine rescaling)
# -------------------------
inv_dt_ok <- inv_dt[ok == TRUE & is.finite(beta_std)]
inv_long_file <- file.path(raw_dir, sprintf("spd_scalar_invariance_long_%s.csv", out_tag))
data.table::fwrite(inv_dt_ok, inv_long_file)

# Wide form for simple differences
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

inv_dt_ok <- unique(inv_dt_ok, by = c("N","gamma_true","rho_meas","rep","scale_id"))
inv_wide <- data.table::dcast(
  inv_dt_ok,
  N + gamma_true + rho_meas + rep ~ scale_id,
  value.var = "beta_std",
  fun.aggregate = safe_mean
)

# If the default scale ids exist, compute differences relative to identity
if (all(c("1","2","3") %in% names(inv_wide))) {
  scale_cols <- c("1","2","3")
  inv_wide[, (scale_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = scale_cols]

  inv_wide[, n_scales_ok := {
    m <- as.matrix(.SD)
    rowSums(is.finite(m))
  }, .SDcols = scale_cols]

  inv_wide[, diff_2_1 := ifelse(is.finite(`2`) & is.finite(`1`), `2` - `1`, NA_real_)]
  inv_wide[, diff_3_1 := ifelse(is.finite(`3`) & is.finite(`1`), `3` - `1`, NA_real_)]
}

# -------------------------
# Save outputs
# -------------------------
rep_file <- file.path(raw_dir, sprintf("spd_scalar_replicates_%s.csv", out_tag))
inv_file <- file.path(raw_dir, sprintf("spd_scalar_invariance_%s.csv", out_tag))
sum_file <- file.path(out_root, sprintf("spd_scalar_summary_%s.csv", out_tag))
ref_file <- file.path(out_root, sprintf("spd_scalar_reference_%s.csv", out_tag))

data.table::fwrite(rep_dt, rep_file)
data.table::fwrite(inv_wide, inv_file)
data.table::fwrite(summ, sum_file)
data.table::fwrite(ref_dt, ref_file)

cfg$log_line(log_file, sprintf("[%s] Done.
", format(Sys.time(), "%H:%M:%S")))
cfg$log_line(log_file, sprintf("Wrote:\n - %s\n - %s\n - %s\n - %s\n - %s\n",
                               rep_file, inv_file, inv_long_file, sum_file, ref_file))

cat("\nSaved outputs to:\n")
cat(sprintf(" - %s\n - %s\n - %s\n - %s\n - %s\n",
            rep_file, inv_file, inv_long_file, sum_file, ref_file))
