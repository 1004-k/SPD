# scripts/01_build_mimic_primary_vasopressor.R
# ------------------------------------------------------------
# Worked example (MIMIC-IV Clinical Database Demo v2.2)
#
# Builds an analysis dataset for:
# - Time zero: ICU admission (icustays.intime)
# - Deviation event (primary): first vasopressor initiation during ICU stay
# - Baseline prognostic score: mortality logistic model LP using age + baseline labs (within 6h)
# ------------------------------------------------------------

source("config/config.R")
source("R/project_utils.R")

cfg <- init_project(
  n_cores = 1L,
  seed    = as.integer(Sys.getenv("SEED", "2026")),
  out_dir = Sys.getenv("OUT_DIR", "output")
)
cfg$require_pkgs(c("data.table"))

message("Loading MIMIC-IV demo tables...")
patients    <- read_mimic(DATA_DIR, "hosp/patients.csv.gz")
admissions  <- read_mimic(DATA_DIR, "hosp/admissions.csv.gz")
icustays    <- read_mimic(DATA_DIR, "icu/icustays.csv.gz")
labevents   <- read_mimic(DATA_DIR, "hosp/labevents.csv.gz")
d_labitems  <- read_mimic(DATA_DIR, "hosp/d_labitems.csv.gz")
inputevents <- read_mimic(DATA_DIR, "icu/inputevents.csv.gz")
d_items     <- read_mimic(DATA_DIR, "icu/d_items.csv.gz")

# ---- cohort: first ICU stay per patient ----
# (demo illustration; adjust for real analyses)
data.table::setDT(icustays)
data.table::setorder(icustays, subject_id, intime)
cohort <- icustays[, .SD[1], by = subject_id]

cohort <- merge(cohort, admissions[, .(hadm_id, hospital_expire_flag)], by = "hadm_id", all.x = TRUE)
cohort <- merge(cohort, patients[, .(subject_id, anchor_age)], by = "subject_id", all.x = TRUE)
data.table::setnames(cohort, "anchor_age", "age")

cohort[, intime := safe_posix(intime)]
cohort[, outtime := safe_posix(outtime)]

# ---- baseline labs within 6 hours ----
data.table::setDT(d_labitems)
d_labitems[, concept := data.table::fifelse(
  grepl("lactate", label, ignore.case = TRUE), "lactate",
  data.table::fifelse(
    grepl("creatinine", label, ignore.case = TRUE), "creatinine",
    data.table::fifelse(
      grepl("\\bwbc\\b|white blood", label, ignore.case = TRUE), "wbc",
      NA_character_
    )
  )
)]

lab_map <- d_labitems[!is.na(concept), .(itemid, concept, label)]
if (nrow(lab_map) == 0) stop("No lab itemids matched lactate/creatinine/WBC patterns.", call. = FALSE)

labs <- merge(
  labevents,
  cohort[, .(subject_id, hadm_id, stay_id, intime)],
  by = c("subject_id", "hadm_id"),
  allow.cartesian = TRUE
)

labs <- labs[itemid %in% lab_map$itemid]
labs <- merge(labs, lab_map[, .(itemid, concept)], by = "itemid", all.x = TRUE)

labs[, charttime := safe_posix(charttime)]
labs <- labs[charttime >= intime & charttime <= intime + BASELINE_LAB_WINDOW_H * 3600]

data.table::setorder(labs, stay_id, concept, charttime)
labs_first <- labs[, .SD[1], by = .(stay_id, concept)]
labs_wide  <- data.table::dcast(labs_first, stay_id ~ concept, value.var = "valuenum")

cohort <- merge(cohort, labs_wide, by = "stay_id", all.x = TRUE)

# median-impute missing labs (demo illustration)
lab_cols <- intersect(c("lactate", "creatinine", "wbc"), names(cohort))
median_impute(cohort, lab_cols)

# age can be missing for some demo records (rare); impute for stable prediction
if ("age" %in% names(cohort)) {
  median_impute(cohort, "age")
}

# ---- prognostic score (LP) ----
# Simple logistic model for in-hospital mortality: age + 3 labs.
fml <- stats::as.formula("hospital_expire_flag ~ age + lactate + creatinine + wbc")
# Note: The demo cohort can contain missing outcomes (hospital_expire_flag).
# glm() will drop those rows by default, which can cause length mismatches when
# assigning predictions back to the full cohort. We therefore:
#   (1) fit on rows with non-missing outcome
#   (2) predict on the full cohort via newdata=
fit_dt <- cohort[!is.na(hospital_expire_flag)]
if (nrow(fit_dt) < 10) {
  stop("Too few non-missing outcomes to fit the prognostic model.", call. = FALSE)
}

# The warning "fitted probabilities 0 or 1" can occur in small demo samples
# (quasi-separation). We keep the simple glm for transparency.
fit_prog <- stats::glm(
  fml,
  data      = fit_dt,
  family    = stats::binomial(),
  na.action = stats::na.exclude
)

cohort[, s_hat := stats::predict(fit_prog, newdata = cohort, type = "link")]
cohort[, z := as.numeric(scale(s_hat))]

# ---- deviation event: first vasopressor initiation during ICU ----
inp <- merge(inputevents, cohort[, .(stay_id, intime, outtime)], by = "stay_id")
inp <- merge(inp, d_items[, .(itemid, label)], by = "itemid", all.x = TRUE)

inp[, starttime := safe_posix(starttime)]
inp[, intime := safe_posix(intime)]
inp[, outtime := safe_posix(outtime)]

vaso <- inp[
  grepl(
    "norepinephrine|noradrenaline|epinephrine|vasopressin|phenylephrine|dopamine",
    label,
    ignore.case = TRUE
  )
]
vaso <- vaso[starttime >= intime & starttime <= outtime]
vaso_t0 <- vaso[, .(t_event = min(starttime)), by = stay_id]

cohort <- merge(cohort, vaso_t0, by = "stay_id", all.x = TRUE)

# ---- time-to-event dataset ----
cohort[, tstop_hr := as.numeric(difftime(outtime, intime, units = "hours"))]
cohort[, deviation_event := !is.na(t_event)]
cohort[, time_hr := data.table::fifelse(
  deviation_event,
  as.numeric(difftime(t_event, intime, units = "hours")),
  tstop_hr
)]

# Backward-compatible aliases (older scripts may reference these)
cohort[, event := deviation_event]
cohort[, t_hr := time_hr]

# minimal cleaning
cohort <- cohort[is.finite(time_hr) & time_hr > 0]
cohort <- cohort[is.finite(tstop_hr) & tstop_hr > 0]

attr(cohort, "spec") <- list(
  time_zero = "ICU admission (icustays.intime)",
  deviation_primary = "first vasopressor initiation during ICU (inputevents; label-based)",
  baseline_labs_window_h = BASELINE_LAB_WINDOW_H,
  prognostic_model = "glm(binomial): hospital_expire_flag ~ age + lactate + creatinine + wbc",
  standardization = "z = scale(s_hat) (per 1 SD)"
)

dir.create(cfg$out_dir, showWarnings = FALSE)
out_path <- file.path(cfg$out_dir, "analysis_dataset_primary_vasopressor.rds")
saveRDS(cohort, file = out_path)
message(sprintf("Saved: %s", out_path))
