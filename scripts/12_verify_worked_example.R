# scripts/12_verify_worked_example.R
# ------------------------------------------------------------
# Recomputes the values reported for the synthetic worked example:
#   Table 4        SPD, HR and bootstrap percentile intervals for the two fitted scores
#   eTable S6      AUC, Brier, calibration, score mean/SD, SPD and intervals (three scores)
#   eTable S7      score-scale sensitivity, all ten conventions
#   contrasts      paired bootstrap differences, gradient boosting minus logistic
#
# The prediction models are not refitted. The script reads the generated
# evaluation cohort and the stored predicted probabilities and recomputes the
# reported values with survival::coxph.
#
# Inputs (default WEX_DIR = results/worked_example):
#   $WEX_DIR/synthetic_evaluation.csv
#   $WEX_DIR/ml_score_bootstrap_raw.csv          (needed for the intervals)
#   $WEX_DIR/ml_score_worked_example_summary.csv (full-precision reported values)
#
# Output: $OUT_DIR/verify_worked_example_R.csv, eTable S7 as score_scale_sensitivity_R.csv
#
# Usage:  Rscript scripts/12_verify_worked_example.R
# ------------------------------------------------------------

source("R/spd_core.R")

WEX_DIR <- Sys.getenv("WEX_DIR", "results/worked_example")
OUT_DIR <- Sys.getenv("OUT_DIR", "output/checks")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

eval_path <- file.path(WEX_DIR, "synthetic_evaluation.csv")
if (!file.exists(eval_path)) {
  stop(sprintf("Could not find %s. Set WEX_DIR to the folder that contains it.", eval_path),
       call. = FALSE)
}
ev <- utils::read.csv(eval_path)

time  <- ev$time
dev   <- ev$deviation
y0    <- ev$reference_outcome
probs <- list(known_risk = ev$oracle_probability,
              linear     = ev$linear_probability,
              boosted    = ev$boosted_probability)
X  <- as.matrix(ev[, paste0("X", 1:5)])
lp <- true_lp(X)
mu <- 0.2 * exp(pmin(pmax(lp, -5), 5))        # known-risk hazard multiplier

cat("Evaluation cohort: n = ", nrow(ev),
    " | deviations = ", sum(dev),
    " | primary outcomes = ", sum(ev$outcome),
    " | reference-treatment events = ", sum(y0), "\n\n", sep = "")

# ============================================================ Table 4 / eTable S6

REPORTED <- data.frame(
  score      = c("Known-DGM reference-risk logit", "Linear logistic", "Gradient boosting"),
  auc        = c(0.778,  0.693,  0.752),
  brier      = c(0.1283, 0.1475, 0.1345),
  cal_int    = c(-0.139, 0.150, -0.230),
  cal_slope  = c(1.030,  1.258,  0.931),
  score_mean = c(-1.390, -1.333, -1.458),
  score_sd   = c(1.353,  0.593,  1.067),
  spd        = c(0.475,  0.349,  0.583),
  boot_low   = c(0.425,  0.281,  0.471),
  boot_high  = c(0.606,  0.415,  0.622),
  stringsAsFactors = FALSE
)

cat("eTable S6 / Table 4 (primary logit scale, probability clip 1e-6)\n")
cat(sprintf("  %-46s %10s %10s %10s  %s\n", "quantity", "reported", "R value", "abs diff", ""))

rows <- list()
for (i in seq_along(probs)) {
  nm <- names(probs)[i]
  p  <- probs[[nm]]
  s  <- sc_logit_clip(p, 1e-6)
  fit <- spd_from_score(time, dev, s, ties = "breslow")
  man <- cox_one_manual(time, dev, standardize(s))
  cal <- calibration(y0, s)

  lab <- REPORTED$score[i]
  check_row(paste(lab, "| AUC"),        REPORTED$auc[i],        auc_roc(y0, p), 5e-4)
  check_row(paste(lab, "| Brier"),      REPORTED$brier[i],      brier(y0, p),   5e-5)
  check_row(paste(lab, "| cal. int."),  REPORTED$cal_int[i],    cal$intercept,  5e-4)
  check_row(paste(lab, "| cal. slope"), REPORTED$cal_slope[i],  cal$slope,      5e-4)
  check_row(paste(lab, "| score mean"), REPORTED$score_mean[i], fit$score_mean, 5e-4)
  check_row(paste(lab, "| score SD"),   REPORTED$score_sd[i],   fit$score_sd,   5e-4)
  check_row(paste(lab, "| SPD"),        REPORTED$spd[i],        fit$spd,        5e-4)
  # coxph versus the independent Newton/step-halving implementation
  check_row(paste(lab, "| coxph vs manual fit"), 0, fit$spd - man$beta, 1e-7)
  cat("\n")

  rows[[nm]] <- data.frame(
    score = lab, auc = auc_roc(y0, p), brier = brier(y0, p),
    calibration_intercept = cal$intercept, calibration_slope = cal$slope,
    score_mean = fit$score_mean, score_sd = fit$score_sd,
    spd = fit$spd, conditional_se = fit$se, hr = fit$hr,
    wald_hr_low = fit$wald_low, wald_hr_high = fit$wald_high,
    stringsAsFactors = FALSE
  )
}
result <- do.call(rbind, rows)

# ============================================================ bootstrap intervals

boot_path <- file.path(WEX_DIR, "ml_score_bootstrap_raw.csv")
if (file.exists(boot_path)) {
  boot <- utils::read.csv(boot_path)
  n_rep <- max(boot$replicate)
  cat(sprintf("Bootstrap percentile intervals from %d stored replicates\n", n_rep))
  # score labels as stored in ml_score_bootstrap_raw.csv
  key <- c(known_risk = "Oracle reference-risk logit",
           linear     = "Linear logistic",
           boosted    = "Gradient boosting")
  for (i in seq_along(key)) {
    v <- boot$SPD[boot$score == key[i]]
    q <- stats::quantile(v, c(0.025, 0.975), names = FALSE, type = 7)
    result$bootstrap_spd_low[i]  <- q[1]
    result$bootstrap_spd_high[i] <- q[2]
    result$bootstrap_hr_low[i]   <- exp(q[1])
    result$bootstrap_hr_high[i]  <- exp(q[2])
    if (n_rep == 2000L) {
      check_row(paste(REPORTED$score[i], "| boot 2.5%"),  REPORTED$boot_low[i],  q[1], 5e-4)
      check_row(paste(REPORTED$score[i], "| boot 97.5%"), REPORTED$boot_high[i], q[2], 5e-4)
    } else {
      cat(sprintf("  %-46s %10s %10.4f, %.4f\n", paste(REPORTED$score[i], "| percentile interval"),
                  "(n/a)", q[1], q[2]))
    }
  }

  # Paired contrasts: gradient boosting minus linear logistic, same replicates.
  cat("\nPaired bootstrap contrasts (boosted minus logistic)\n")
  wide <- function(metric) {
    a <- boot[boot$score == key["boosted"], c("replicate", metric)]
    b <- boot[boot$score == key["linear"],  c("replicate", metric)]
    m <- merge(a, b, by = "replicate", suffixes = c("_boost", "_lin"))
    m[[paste0(metric, "_boost")]] - m[[paste0(metric, "_lin")]]
  }
  contrasts <- do.call(rbind, lapply(c("AUC", "Brier", "SPD"), function(metric) {
    v <- wide(metric)
    q <- stats::quantile(v, c(0.025, 0.975), names = FALSE, type = 7)
    point <- result[3, tolower(metric)] - result[2, tolower(metric)]
    if (metric == "SPD") point <- result$spd[3] - result$spd[2]
    cat(sprintf("  %-46s %10.4f  [%.4f, %.4f]\n", metric, point, q[1], q[2]))
    data.frame(metric = metric, contrast = "boosted_minus_linear",
               point = point, bootstrap_low = q[1], bootstrap_high = q[2],
               stringsAsFactors = FALSE)
  }))
  utils::write.csv(contrasts, file.path(OUT_DIR, "bootstrap_contrasts_R.csv"), row.names = FALSE)
  cat("\n")
} else {
  cat("ml_score_bootstrap_raw.csv not found; skipping interval reproduction.\n\n")
}

# ============================================================ full-precision check

summ_path <- file.path(WEX_DIR, "ml_score_worked_example_summary.csv")
if (file.exists(summ_path)) {
  stored <- utils::read.csv(summ_path)
  stored <- stored[match(c("Oracle reference-risk logit", "Linear logistic", "Gradient boosting"),
                         stored$score), ]
  pairs <- c(AUC = "auc", Brier = "brier",
             calibration_intercept = "calibration_intercept",
             calibration_slope = "calibration_slope",
             score_mean = "score_mean", score_sd = "score_sd",
             SPD = "spd", conditional_SE = "conditional_se")
  if ("bootstrap_spd_low" %in% names(result)) {
    pairs <- c(pairs, bootstrap_SPD_low = "bootstrap_spd_low", bootstrap_SPD_high = "bootstrap_spd_high")
  }
  diffs <- vapply(names(pairs), function(k) max(abs(stored[[k]] - result[[pairs[[k]]]])), numeric(1))
  cat("Largest absolute difference from the stored full-precision values
")
  for (k in names(diffs)) cat(sprintf("  %-46s %10.1e\n", k, diffs[[k]]))
  cat("\n")
}

utils::write.csv(result, file.path(OUT_DIR, "verify_worked_example_R.csv"), row.names = FALSE)

# ============================================================ eTable S7

exact <- list(
  known_risk = sc_logit_exact_known(mu),
  linear     = log(probs$linear)  - log1p(-probs$linear),
  boosted    = log(probs$boosted) - log1p(-probs$boosted)
)

conventions <- list()
for (eps in c(1e-4, 1e-5, 1e-6, 1e-8, 1e-10, 1e-12)) {
  conventions[[sprintf("Probability clip %g", eps)]] <-
    lapply(probs, sc_logit_clip, eps = eps)
}
conventions[["Exact logit, no probability clipping"]] <- exact
conventions[["Logit winsorized at own 0.5/99.5 percentiles"]] <- lapply(exact, sc_winsorize)
conventions[["Untransformed probability"]] <- lapply(probs, sc_probability)
conventions[["Cloglog, probability clip 1e-6"]] <- lapply(probs, sc_cloglog, eps = 1e-6)

S7_REPORTED <- data.frame(
  known_risk = c(0.576, 0.526, 0.475, 0.417, 0.394, 0.368, 0.353, 0.635, 0.617, 0.662),
  linear     = c(0.349, 0.349, 0.349, 0.349, 0.349, 0.349, 0.349, 0.348, 0.375, 0.340),
  boosted    = c(0.583, 0.583, 0.583, 0.583, 0.583, 0.583, 0.583, 0.576, 0.553, 0.565)
)

cat("eTable S7: score-scale sensitivity\n")
cat(sprintf("  %-46s %21s %21s\n", "", "reported", "R value"))
cat(sprintf("  %-46s %7s%7s%7s %7s%7s%7s\n", "convention",
            "known", "linear", "boost", "known", "linear", "boost"))
s7 <- list(); worst <- 0
for (k in seq_along(conventions)) {
  nm <- names(conventions)[k]
  vals <- vapply(c("known_risk", "linear", "boosted"),
                 function(j) spd_from_score(time, dev, conventions[[nm]][[j]], ties = "breslow")$spd,
                 numeric(1))
  worst <- max(worst, max(abs(vals - unlist(S7_REPORTED[k, ]))))
  cat(sprintf("  %-46s %7.3f%7.3f%7.3f %7.3f%7.3f%7.3f\n", nm,
              S7_REPORTED$known_risk[k], S7_REPORTED$linear[k], S7_REPORTED$boosted[k],
              vals[1], vals[2], vals[3]))
  s7[[k]] <- data.frame(convention = nm, known_risk = vals[1],
                        linear = vals[2], boosted = vals[3], stringsAsFactors = FALSE)
}
utils::write.csv(do.call(rbind, s7), file.path(OUT_DIR, "score_scale_sensitivity_R.csv"),
                 row.names = FALSE)
cat(sprintf("\n  largest absolute difference from the reported eTable S7: %.2e\n", worst))
s7_stored <- file.path(Sys.getenv("TIMING_DIR", "results/additional_checks"), "score_scale_sensitivity.csv")
if (file.exists(s7_stored)) {
  st7 <- utils::read.csv(s7_stored)
  r7  <- do.call(rbind, s7)
  if (nrow(st7) == nrow(r7)) {
    cols <- c("known_risk", "linear", "boosted")
    d7 <- max(abs(as.matrix(st7[, cols]) - as.matrix(r7[, cols])))
    cat(sprintf("  largest absolute difference from the stored full-precision values: %.1e\n", d7))
  }
}
cat(sprintf("\nWrote R verification output to %s\n", OUT_DIR))
