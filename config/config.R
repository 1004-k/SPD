# config/config.R
# -----------------------------
# Local configuration for the MIMIC-IV Demo worked example.
#
# Download the MIMIC-IV Clinical Database Demo (v2.2) ZIP from PhysioNet,
# unzip locally, then set DATA_DIR to the unzipped folder.
#
# Expected structure:
# data/mimic-iv-clinical-database-demo-2.2/
#   hosp/patients.csv.gz
#   hosp/admissions.csv.gz
#   hosp/labevents.csv.gz
#   hosp/d_labitems.csv.gz
#   icu/icustays.csv.gz
#   icu/inputevents.csv.gz
#   icu/d_items.csv.gz

DATA_DIR <- Sys.getenv(
  "MIMIC_DATA_DIR",
  "data/mimic-iv-clinical-database-demo-2.2"
)

# Baseline lab window relative to ICU admission (hours)
BASELINE_LAB_WINDOW_H <- as.numeric(Sys.getenv("BASELINE_LAB_WINDOW_H", "6"))
