# scripts/03_fit_spd_sensitivity_norepinephrine.R
# ------------------------------------------------------------
# Worked example sensitivity analysis:
# Deviation event = first norepinephrine (or noradrenaline) initiation during ICU stay.
# ------------------------------------------------------------

source("config/config.R")
source("R/project_utils.R")

cfg <- init_project(
  n_cores = 1L,
  seed    = as.integer(Sys.getenv("SEED", "2026")),
  out_dir = Sys.getenv("OUT_DIR", "output")
)
cfg$require_pkgs(c("data.table","survival"))

dat0 <- readRDS(file.path("output", "analysis_dataset_primary_vasopressor.rds"))

# Load MIMIC tables needed to define norepinephrine-only initiation
inputevents <- read_mimic(DATA_DIR, "icu/inputevents.csv.gz")
d_items     <- read_mimic(DATA_DIR, "icu/d_items.csv.gz")

cohort <- dat0[, .(stay_id, intime, outtime)]
cohort[, intime  := safe_posix(intime)]
cohort[, outtime := safe_posix(outtime)]

inp <- merge(inputevents, cohort, by = "stay_id")
inp <- merge(inp, d_items[, .(itemid, label)], by = "itemid", all.x = TRUE)

inp[, starttime := safe_posix(starttime)]
inp <- inp[starttime >= intime & starttime <= outtime]

norepi <- inp[grepl("norepinephrine|noradrenaline", label, ignore.case = TRUE)]
norepi_t0 <- norepi[, .(t_event_nr = min(starttime)), by = stay_id]

dat <- merge(dat0, norepi_t0, by = "stay_id", all.x = TRUE)

# Time-to-event (hours since ICU admission)
dat[, tstop_hr := as.numeric(difftime(outtime, intime, units = "hours"))]
dat[, deviation_event := !is.na(t_event_nr)]
dat[, time_hr := data.table::fifelse(
  deviation_event,
  as.numeric(difftime(t_event_nr, intime, units = "hours")),
  tstop_hr
)]

dat <- dat[is.finite(time_hr) & time_hr > 0]
dat <- dat[is.finite(tstop_hr) & tstop_hr > 0]

fit_z <- survival::coxph(survival::Surv(time_hr, deviation_event) ~ z, data = dat)
spd <- unname(stats::coef(fit_z)["z"])
hr_1sd <- exp(spd)
ci <- exp(stats::confint(fit_z)["z", ])

# Optional invariance check (raw LP coefficient * SD(LP))
fit_lp <- survival::coxph(survival::Surv(time_hr, deviation_event) ~ s_hat, data = dat)
coef_lp <- unname(stats::coef(fit_lp)["s_hat"])
spd_from_lp <- coef_lp * stats::sd(dat$s_hat)

res <- data.table::data.table(
  dataset = "mimic-iv-clinical-database-demo-2.2",
  event_definition = "sensitivity_norepinephrine_initiation",
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
data.table::fwrite(res, file.path("output", "results_sensitivity_norepinephrine.csv"))
saveRDS(res, file.path("output", "results_sensitivity_norepinephrine.rds"))

message("Saved: output/results_sensitivity_norepinephrine.csv and .rds")
print(res)
