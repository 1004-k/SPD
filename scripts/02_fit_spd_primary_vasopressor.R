# scripts/02_fit_spd_primary_vasopressor.R
# ------------------------------------------------------------
# Worked example: estimate scalar SPD for primary vasopressor initiation.
# Model: Surv(time_hr, deviation_event) ~ z
# where z is a 1-SD standardized baseline prognostic score.
# ------------------------------------------------------------

source("R/project_utils.R")

cfg <- init_project(
  n_cores = 1L,
  seed    = as.integer(Sys.getenv("SEED", "2026")),
  out_dir = Sys.getenv("OUT_DIR", "output")
)
cfg$require_pkgs(c("data.table","survival"))

dat <- readRDS(file.path("output", "analysis_dataset_primary_vasopressor.rds"))

# Backward-compatible field names
data.table::setDT(dat)
if (!"time_hr" %in% names(dat) && "t_hr" %in% names(dat)) dat[, time_hr := t_hr]
if (!"deviation_event" %in% names(dat) && "event" %in% names(dat)) dat[, deviation_event := event]

req_cols <- c("time_hr", "deviation_event", "z", "s_hat")
miss <- setdiff(req_cols, names(dat))
if (length(miss) > 0) {
  stop(
    sprintf("Missing required columns in analysis dataset: %s", paste(miss, collapse = ", ")),
    call. = FALSE
  )
}

fit_z <- survival::coxph(survival::Surv(time_hr, deviation_event) ~ z, data = dat)
spd <- unname(stats::coef(fit_z)["z"])
hr_1sd <- exp(spd)
ci <- exp(stats::confint(fit_z)["z", ])

# Optional invariance check: coefficient on raw LP times SD(LP) should match SPD
fit_lp <- survival::coxph(survival::Surv(time_hr, deviation_event) ~ s_hat, data = dat)
coef_lp <- unname(stats::coef(fit_lp)["s_hat"])
spd_from_lp <- coef_lp * stats::sd(dat$s_hat)

res <- data.table::data.table(
  dataset = "mimic-iv-clinical-database-demo-2.2",
  event_definition = "primary_any_vasopressor_initiation",
  n = nrow(dat),
  events = sum(dat$deviation_event),
  spd = spd,
  HR_per_1SD = hr_1sd,
  CI_low = ci[1],
  CI_high = ci[2],
  coef_lp = coef_lp,
  sd_lp = stats::sd(dat$s_hat),
  spd_from_lp = spd_from_lp
)

dir.create("output", showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(res, file.path("output", "results_primary_vasopressor.csv"))
saveRDS(res, file.path("output", "results_primary_vasopressor.rds"))

message("Saved: output/results_primary_vasopressor.csv and .rds")
print(res)
