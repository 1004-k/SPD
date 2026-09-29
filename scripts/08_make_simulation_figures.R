# scripts/08_make_simulation_figures.R
# ------------------------------------------------------------
# Supplementary figures for the estimator simulations, drawn in base R from
# the output of scripts/07 (baseline setting) and scripts/09.
#
#   eFigure S1  model-based versus robust SE coverage (gamma = 1, rho = 1)
#   eFigure S2  affine-invariance check
#   eFigure S3  bias of the SPD estimate relative to beta_ref
#   eFigure S4  coverage of nominal 95% Wald intervals
#
# Inputs (default SIM_DIR = output):
#   $SIM_DIR/baseline/spd_scalar_summary_baseline.csv
#   $SIM_DIR/baseline/raw/spd_scalar_invariance_baseline.csv
#   $SIM_DIR/sim_spd_scalar_robustSE_worstcase/
#       spd_scalar_robustSE_summary_sim_spd_scalar_robustSE_worstcase.csv
# Output: $FIG_DIR/eFigureS1-S4 (TIFF and PDF)
#
# Usage:  Rscript scripts/08_make_simulation_figures.R
# ------------------------------------------------------------

SIM_DIR    <- Sys.getenv("SIM_DIR", "output")
SIM_TAG    <- Sys.getenv("SIM_TAG", "baseline")
ROBUST_TAG <- Sys.getenv("ROBUST_TAG", "sim_spd_scalar_robustSE_worstcase")
FIG_DIR    <- Sys.getenv("FIG_DIR", "output/figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(SIM_DIR, SIM_TAG, sprintf("spd_scalar_summary_%s.csv", SIM_TAG))
inv_path     <- file.path(SIM_DIR, SIM_TAG, "raw", sprintf("spd_scalar_invariance_%s.csv", SIM_TAG))
robust_path  <- file.path(SIM_DIR, ROBUST_TAG,
                          sprintf("spd_scalar_robustSE_summary_%s.csv", ROBUST_TAG))

save_fig <- function(name, draw, width, height, res = 600) {
  grDevices::tiff(file.path(FIG_DIR, paste0(name, ".tiff")), width = width, height = height,
                  units = "in", res = res, compression = "lzw", type = "cairo")
  draw(); invisible(grDevices::dev.off())
  grDevices::pdf(file.path(FIG_DIR, paste0(name, ".pdf")), width = width, height = height)
  draw(); invisible(grDevices::dev.off())
  cat(sprintf("Wrote %s (.tiff and .pdf)\n", file.path(FIG_DIR, name)))
}

n_style <- data.frame(N = c(200, 800, 2000), lty = c(1, 2, 3), pch = c(19, 17, 15))

# ---------------------------------------------------------------- eFigures S3, S4

panel_plot <- function(summ, column, ylab, ylim, ref_line) {
  op <- graphics::par(mfrow = c(1, 3), mar = c(4.2, 3.2, 2.2, 0.6), oma = c(0, 2.6, 0, 0),
                      mgp = c(2.4, 0.7, 0), las = 1)
  on.exit(graphics::par(op), add = TRUE)
  for (r in sort(unique(summ$rho_meas))) {
    d <- summ[summ$rho_meas == r, ]
    plot(NA, xlim = range(d$gamma_true), ylim = ylim,
         xlab = expression(gamma[true]), ylab = "", main = bquote(rho == .(sprintf("%.1f", r))))
    graphics::abline(h = ref_line, lty = 3, col = "grey60")
    for (k in seq_len(nrow(n_style))) {
      dk <- d[d$N == n_style$N[k], ]
      dk <- dk[order(dk$gamma_true), ]
      graphics::lines(dk$gamma_true, dk[[column]], lty = n_style$lty[k], lwd = 2)
      graphics::points(dk$gamma_true, dk[[column]], pch = n_style$pch[k], cex = 0.9)
    }
    graphics::legend("bottomleft", legend = paste0("N = ", n_style$N), lty = n_style$lty,
                     pch = n_style$pch, lwd = 2, bty = "n", cex = 0.9)
  }
  graphics::mtext(ylab, side = 2, outer = TRUE, line = 0.8, las = 0, cex = 0.95)
}

if (file.exists(summary_path)) {
  summ <- utils::read.csv(summary_path)
  lim_b <- max(abs(summ$bias), na.rm = TRUE) * 1.15
  save_fig("eFigureS3_bias", function() {
    panel_plot(summ, "bias",
               expression("Bias (estimated SPD" - beta[ref] * ")"),
               ylim = c(-lim_b, lim_b), ref_line = 0)
  }, width = 8, height = 5.2)
  save_fig("eFigureS4_coverage", function() {
    panel_plot(summ, "coverage", "Coverage of nominal 95% Wald intervals",
               ylim = range(c(summ$coverage, 0.95), na.rm = TRUE) + c(-0.01, 0.01),
               ref_line = 0.95)
  }, width = 8, height = 5.2)
} else {
  cat(sprintf("%s not found; run scripts/07 first (OUT_TAG=%s). eFigures S3 and S4 skipped.\n",
              summary_path, SIM_TAG))
}

# ---------------------------------------------------------------- eFigure S2

if (file.exists(inv_path)) {
  inv <- utils::read.csv(inv_path)
  d <- c(inv$diff_2_1, inv$diff_3_1)
  d <- d[is.finite(d)]
  save_fig("eFigureS2_affine_invariance", function() {
    op <- graphics::par(mar = c(4.4, 5.8, 1.0, 1.0), mgp = c(4.0, 0.8, 0), las = 1)
    on.exit(graphics::par(op), add = TRUE)
    graphics::hist(d, breaks = 50, main = "", col = "white", ylab = "Frequency",
                   xlab = expression("Difference in SPD estimate (affine" - "identity)"))
    graphics::abline(v = 0, lty = 3)
  }, width = 7.2, height = 4.8)
  cat(sprintf("  affine-invariance differences: n = %d, max |difference| = %.1e\n",
              length(d), max(abs(d))))
} else {
  cat(sprintf("%s not found; run scripts/07 first. eFigure S2 skipped.\n", inv_path))
}

# ---------------------------------------------------------------- eFigure S1

if (file.exists(robust_path)) {
  rb <- utils::read.csv(robust_path)
  rb <- rb[order(rb$setting, rb$N), ]
  series <- list(
    list(setting = "baseline", col = "coverage_model",  lty = 1, pch = 19, lab = "Baseline, model-based SE"),
    list(setting = "baseline", col = "coverage_robust", lty = 2, pch = 1,  lab = "Baseline, robust SE"),
    list(setting = "fastdev",  col = "coverage_model",  lty = 1, pch = 17, lab = "Faster deviation, model-based SE"),
    list(setting = "fastdev",  col = "coverage_robust", lty = 2, pch = 2,  lab = "Faster deviation, robust SE")
  )
  yr <- range(c(rb$coverage_model, rb$coverage_robust, 0.95), na.rm = TRUE) + c(-0.01, 0.01)
  save_fig("eFigureS1_robust_se_coverage", function() {
    op <- graphics::par(mar = c(4.4, 4.6, 1.0, 1.0), mgp = c(2.8, 0.8, 0), las = 1)
    on.exit(graphics::par(op), add = TRUE)
    plot(NA, xlim = range(rb$N), ylim = yr, xaxt = "n",
         xlab = "Sample size (N)", ylab = "Coverage (nominal 95%)")
    graphics::axis(1, at = sort(unique(rb$N)))
    graphics::abline(h = 0.95, lty = 3, col = "grey60")
    for (s in series) {
      d <- rb[rb$setting == s$setting, ]
      graphics::lines(d$N, d[[s$col]], lty = s$lty, lwd = 2)
      graphics::points(d$N, d[[s$col]], pch = s$pch)
    }
    graphics::legend("bottomright", legend = vapply(series, `[[`, "", "lab"),
                     lty = vapply(series, `[[`, 0, "lty"), pch = vapply(series, `[[`, 0, "pch"),
                     lwd = 2, bty = "n", cex = 0.85)
  }, width = 7.2, height = 4.8)
} else {
  cat(sprintf("%s not found; run scripts/09 first. eFigure S1 skipped.\n", robust_path))
}
