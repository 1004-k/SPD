# scripts/16_verify_complementarity.R
# ------------------------------------------------------------
# Recomputes Table 3 and eTable S5 from the stored replicate-level results of
# the diagnostic-complementarity simulation (5 scenarios x 500 replicates).
#
# Inputs (default COMP_DIR = results/complementarity):
#   $COMP_DIR/diagnostic_complementarity_raw.csv      one row per replicate
#   $COMP_DIR/diagnostic_complementarity_summary.csv  stored summary (optional)
# Output: $OUT_DIR/complementarity_summary_R.csv
#
# scripts/11 re-runs the simulation itself; this script only summarizes the
# stored replicates.
#
# Usage:  Rscript scripts/16_verify_complementarity.R
# ------------------------------------------------------------

source("R/spd_core.R")   # check_row()

COMP_DIR <- Sys.getenv("COMP_DIR", "results/complementarity")
OUT_DIR  <- Sys.getenv("OUT_DIR", "output/checks")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

raw_path <- file.path(COMP_DIR, "diagnostic_complementarity_raw.csv")
if (!file.exists(raw_path)) {
  stop(sprintf("Could not find %s. Set COMP_DIR to the folder that contains it.", raw_path),
       call. = FALSE)
}
raw <- utils::read.csv(raw_path)

scenarios <- c("A_positive_selection", "A_negative_selection", "B_high_spd_treatment_related",
               "C_low_observed_spd", "D_near_zero_baseline_spd")
metrics <- c("deviation_rate", "spd", "spd_hr", "abs_covariate_smd", "prognostic_score_smd",
             "rESS_percent", "top1_weight_share", "naive_outcome_loghr", "naive_coefficient_gap")

summ <- do.call(rbind, lapply(scenarios, function(s) {
  d <- raw[raw$scenario == s, ]
  out <- data.frame(scenario = s, gamma = d$gamma[1], alpha = d$alpha[1], rho = d$rho[1],
                    n_rep = nrow(d),
                    spd_failure_rate = 1 - mean(as.logical(d$spd_ok)),
                    outcome_failure_rate = 1 - mean(as.logical(d$outcome_ok)),
                    stringsAsFactors = FALSE)
  for (m in metrics) {
    v <- d[[m]][is.finite(d[[m]])]
    out[[paste0(m, "_mean")]] <- mean(v)
    out[[paste0(m, "_sd")]]   <- stats::sd(v)
    out[[paste0(m, "_mcse")]] <- stats::sd(v) / sqrt(length(v))
  }
  out
}))
utils::write.csv(summ, file.path(OUT_DIR, "complementarity_summary_R.csv"), row.names = FALSE)

# Values as printed in Table 3 of the manuscript
REPORTED <- data.frame(
  label     = c("A+: positive", "A-: negative", "B: prognosis + treatment",
                "C: noisy score", "D: treatment only"),
  dev_rate  = c(0.496, 0.513, 0.649, 0.648, 0.721),
  spd       = c(1.002, -1.002, 1.205, 0.242, 0.006),
  spd_mcse  = c(0.0018, 0.0017, 0.0019, 0.0012, 0.0012),
  cov_smd   = c(0.780, 0.956, 0.970, 0.964, 0.124),
  score_smd = c(1.020, -1.280, 1.289, 0.330, -0.146),
  ress      = c(94.2, 92.9, 87.3, 87.3, 79.5),
  gap       = c(0.006, 0.028, -0.311, -0.333, 0.001),
  gap_mcse  = c(0.0065, 0.0051, 0.0097, 0.0104, 0.0091),
  stringsAsFactors = FALSE
)

cat("Table 3: diagnostic-complementarity simulation, recomputed from stored replicates\n")
cat(sprintf("  %-46s %10s %10s %10s  %s\n", "quantity", "reported", "R value", "abs diff", ""))
for (i in seq_along(scenarios)) {
  r <- summ[i, ]; p <- REPORTED[i, ]; lab <- p$label
  check_row(paste(lab, "| deviation rate"), p$dev_rate,  r$deviation_rate_mean,        5e-4)
  check_row(paste(lab, "| mean SPD"),       p$spd,       r$spd_mean,                   5e-4)
  check_row(paste(lab, "| SPD MCSE"),       p$spd_mcse,  r$spd_mcse,                   5e-5)
  check_row(paste(lab, "| max |cov. SMD|"), p$cov_smd,   r$abs_covariate_smd_mean,     5e-4)
  check_row(paste(lab, "| score SMD"),      p$score_smd, r$prognostic_score_smd_mean,  5e-4)
  check_row(paste(lab, "| rESSexit (%)"),   p$ress,      r$rESS_percent_mean,          5e-2)
  check_row(paste(lab, "| coefficient gap"), p$gap,      r$naive_coefficient_gap_mean, 5e-4)
  check_row(paste(lab, "| gap MCSE"),       p$gap_mcse,  r$naive_coefficient_gap_mcse, 5e-5)
  cat("\n")
}

stored_path <- file.path(COMP_DIR, "diagnostic_complementarity_summary.csv")
if (file.exists(stored_path)) {
  st <- utils::read.csv(stored_path)
  st <- st[match(scenarios, st$scenario), ]
  cols <- intersect(setdiff(names(summ), c("scenario", "gamma", "alpha", "rho")), names(st))
  d <- max(abs(as.matrix(st[, cols]) - as.matrix(summ[, cols])))
  cat(sprintf("Largest absolute difference from the stored summary (%d columns): %.1e\n",
              length(cols), d))
}
cat(sprintf("Wrote %s\n", file.path(OUT_DIR, "complementarity_summary_R.csv")))
