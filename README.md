# Scalar SPD (standardized prognostic dependence): simulations + worked example

This repository contains **GitHub-ready code** to reproduce the *scalar* SPD results and a small
worked example using **MIMIC-IV Clinical Database Demo (v2.2)**.

**Scalar SPD** is defined as the Cox regression coefficient for deviation per **1 SD** higher
baseline prognostic score. This makes the index **unit-free** and invariant to affine rescaling
of the underlying prognostic score.

## Repository structure

- `R/`
  - `project_utils.R`: logging, package checks, MIMIC table I/O helpers
  - `spd_scalar.R`: scalar SPD estimation + reference computation
- `config/config.R`: local configuration for the MIMIC demo (data path, baseline window)
- `scripts/`
  - `01_build_mimic_primary_vasopressor.R`: build worked-example analysis dataset
  - `02_fit_spd_primary_vasopressor.R`: fit Cox model for primary definition (any vasopressor)
  - `03_fit_spd_sensitivity_norepinephrine.R`: sensitivity (norepinephrine-only)
  - `07_run_sim_spd_scalar.R`: main simulation grid (bias/coverage + invariance)
  - `09_run_sim_spd_scalar_robust_se_worstcase.R`: worst-case robust-SE check
  - `10_run_sim_spd_bias_amplifier.R`: illustrative "amplifier" simulation (SPD vs bias)
- `data/`: **empty by default** (place MIMIC demo here)
- `output/`: created automatically (ignored by git)

## Requirements

R packages used across scripts:
- `data.table`
- `survival`
- `future`
- `future.apply`

Install them in R:
```r
install.packages(c("data.table","survival","future","future.apply"))
```

## Worked example (MIMIC-IV Demo)

### Data (not included)

Download and unzip **MIMIC-IV Clinical Database Demo (v2.2)** from PhysioNet.

Place the unzipped folder at:
```
data/mimic-iv-clinical-database-demo-2.2/
  hosp/...
  icu/...
```

Or set the environment variable:
- `MIMIC_DATA_DIR=/path/to/mimic-iv-clinical-database-demo-2.2`

### Run

From the repository root:
```r
source("scripts/01_build_mimic_primary_vasopressor.R")
source("scripts/02_fit_spd_primary_vasopressor.R")
source("scripts/03_fit_spd_sensitivity_norepinephrine.R")
```

Outputs:
- `output/analysis_dataset_primary_vasopressor.rds`
- `output/results_primary_vasopressor.csv` (main numbers)
- `output/results_sensitivity_norepinephrine.csv` (sensitivity)

## Simulations

### Main scalar SPD simulation grid

```bash
Rscript scripts/07_run_sim_spd_scalar.R
```

Key environment variables (optional):
- `OUT_DIR` (default: `output`)
- `OUT_TAG` (default: `sim_spd_scalar_B{B}`)
- `B` (default: 200)
- `N_VEC` (default: `200,800,2000`)
- `GAMMA_VEC` (default: `0,0.25,0.5,1.0`)
- `RHO_VEC` (default: `1.0,0.7,0.4`)
- `DO_REF` (default: 1) and `N_REF`, `B_REF` for pseudo-true reference computation
- `N_CORES` (default: 3)

### Worst-case robust-SE check

```bash
Rscript scripts/09_run_sim_spd_scalar_robust_se_worstcase.R
```

### Illustrative amplifier simulation

```bash
Rscript scripts/10_run_sim_spd_bias_amplifier.R
```

## Notes

- This code is written to be readable and reproducible rather than maximally optimized.
- The MIMIC demo cohort is small; treat the worked example as an illustration only.
