# run_all_demo.R
# ------------------------------------------------------------
# Convenience wrapper: run the MIMIC-IV Demo worked example end-to-end.
# ------------------------------------------------------------

source("scripts/01_build_mimic_primary_vasopressor.R")
source("scripts/02_fit_spd_primary_vasopressor.R")
source("scripts/03_fit_spd_sensitivity_norepinephrine.R")
