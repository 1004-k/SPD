# scripts/14_make_figures.R
# ------------------------------------------------------------
# Main-manuscript figures in base R.
#
#   Figure 1  Selection schematic (A -> D <- Z -> Y, A -> Y) and reporting
#             workflow. No data; drawn from the layout below.
#   Figure 2  Diagnostic-complementarity simulation. Panel A: mean SPD with
#             +/- 1.96 MCSE bars. Panel B: mean SPD against the outcome
#             coefficient gap.
#
# Input for Figure 2: $COMP_DIR/diagnostic_complementarity_summary.csv
#   default COMP_DIR = results/complementarity (the reported results);
#   set COMP_DIR to the output folder of scripts/11 to plot a re-run.
# Output: $FIG_DIR/Figure1_R.tiff, $FIG_DIR/Figure2_R.tiff (plus PDFs)
#
# Usage:  Rscript scripts/14_make_figures.R
# ------------------------------------------------------------

COMP_DIR <- Sys.getenv("COMP_DIR", "results/complementarity")
FIG_DIR  <- Sys.getenv("FIG_DIR", "output/figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(name, draw, width, height, res = 400) {
  grDevices::tiff(file.path(FIG_DIR, paste0(name, ".tiff")), width = width, height = height,
                  units = "in", res = res, compression = "lzw", type = "cairo")
  draw(); grDevices::dev.off()
  grDevices::pdf(file.path(FIG_DIR, paste0(name, ".pdf")), width = width, height = height)
  draw(); grDevices::dev.off()
  cat(sprintf("Wrote %s.tiff and %s.pdf\n", file.path(FIG_DIR, name), file.path(FIG_DIR, name)))
}

box <- function(x, y, w, h, label, fill, cex = 1) {
  graphics::rect(x, y, x + w, y + h, col = fill, border = "#555555", lwd = 1.2)
  graphics::text(x + w / 2, y + h / 2, label, cex = cex)
}

arrow <- function(x0, y0, x1, y1, col = "#555555") {
  graphics::arrows(x0, y0, x1, y1, length = 0.10, lwd = 1.8, col = col)
}

# ---------------------------------------------------------------- Figure 1

draw_figure1 <- function() {
  op <- graphics::par(mfrow = c(2, 1), mar = c(0.4, 0.6, 0.6, 0.6), xaxs = "i", yaxs = "i")
  on.exit(graphics::par(op), add = TRUE)

  # Panel A: selection schematic
  plot(NA, xlim = c(0, 10), ylim = c(0, 4), axes = FALSE, xlab = "", ylab = "")
  graphics::text(0.1, 3.75, "A. A possible selection mechanism", adj = c(0, 1), font = 2, cex = 1.15)
  box(1.00, 2.00, 1.75, 0.78, "Treatment A",            "#e8f1fb")
  box(4.20, 2.00, 1.75, 0.78, "Deviation D",            "#fff0dc")
  box(7.40, 2.00, 1.75, 0.78, "Prognosis Z\nscore S(X)", "#eaf4eb")
  box(4.20, 0.25, 1.75, 0.78, "Outcome Y",              "#eee7f5")
  arrow(2.80, 2.39, 4.15, 2.39)                    # A -> D
  arrow(7.35, 2.39, 6.00, 2.39, col = "#a14b31")   # Z -> D
  arrow(1.90, 1.94, 4.15, 0.72)                    # A -> Y
  arrow(8.25, 1.94, 6.00, 0.72)                    # Z -> Y
  graphics::rect(4.07, 1.87, 6.08, 2.91, border = "#a14b31", lwd = 1.1, lty = 2)
  graphics::text(5.08, 3.16, "Conditioning on adherence can open A -> D <- Z -> Y",
                 cex = 0.95)
  graphics::text(5.00, 0.08, adj = c(0.5, 0), cex = 0.9,
                 "SPD describes the score-deviation association; it is not an estimate of a causal arrow.")

  # Panel B: reporting workflow
  plot(NA, xlim = c(0, 10), ylim = c(0, 3), axes = FALSE, xlab = "", ylab = "")
  graphics::text(0.1, 2.80, "B. Reporting and follow-up", adj = c(0, 1), font = 2, cex = 1.15)
  steps <- list(
    c("Define",  "Reference prognosis,\nscore scale and\ndeviation rules"),
    c("Fit",     "Standardize score;\ncause-specific\ndeviation model"),
    c("Report",  "Signed SPD, HR,\nuncertainty and\nscore quality"),
    c("Inspect", "Balance, overlap,\nconditional weights,\npost-baseline drivers"))
  xs <- c(0.05, 2.60, 5.15, 7.70)
  for (i in seq_along(steps)) {
    graphics::rect(xs[i], 0.65, xs[i] + 2.2, 2.35, col = "#f1f3f5", border = "#666666", lwd = 1)
    graphics::text(xs[i] + 1.1, 2.02, steps[[i]][1], font = 2, cex = 1.0)
    graphics::text(xs[i] + 1.1, 1.30, steps[[i]][2], cex = 0.88)
    if (i < length(steps)) arrow(xs[i] + 2.25, 1.5, xs[i] + 2.52, 1.5)
  }
  graphics::text(5.0, 0.22, adj = c(0.5, 0), cex = 0.9,
                 "No universal threshold; keep the planned causal design and adjustment.")
}

save_fig("Figure1_R", draw_figure1, width = 10.5, height = 6.2)

# ---------------------------------------------------------------- Figure 2

summ_path <- file.path(COMP_DIR, "diagnostic_complementarity_summary.csv")
have_summary <- file.exists(summ_path)
if (!have_summary) {
  cat(sprintf("\n%s not found; Figure 1 was written, Figure 2 skipped.\n", summ_path))
  cat("Run scripts/11_run_sim_spd_diagnostic_complementarity.R first, then set COMP_DIR.\n")
}
d <- if (have_summary) utils::read.csv(summ_path) else NULL

# Scenario order as in Table 3.
if (have_summary) {
ord <- c("A_positive_selection", "A_negative_selection", "B_high_spd_treatment_related",
         "C_low_observed_spd", "D_near_zero_baseline_spd")
d <- d[match(ord, d$scenario), ]
labels  <- c("A+\npositive", "A-\nnegative", "B\nprognosis +\ntreatment",
             "C\nnoisy score", "D\ntreatment\nonly")
short   <- c("A+", "A-", "B", "C", "D")
cols    <- c("#2d6a9f", "#b55d44", "#8d4f9e", "#d28b2c", "#4c956c")
}

draw_figure2 <- function() {
  op <- graphics::par(mfrow = c(1, 2), mar = c(5.0, 4.6, 4.0, 1.0),
                      mgp = c(2.7, 0.8, 0), las = 1, oma = c(2.4, 0, 2.2, 0))
  on.exit(graphics::par(op), add = TRUE)
  x <- seq_along(ord)

  plot(NA, xlim = c(0.6, 5.4), ylim = c(-1.3, 1.5), xaxt = "n",
       xlab = "", ylab = "Mean SPD (log-HR per 1-SD)",
       main = "A. Prognostic sorting and score reliability",
       cex.main = 1.05, font.main = 1)
  graphics::grid(nx = NA, ny = NULL, col = "grey88", lty = 1)
  graphics::abline(h = 0, col = "#666666", lwd = 0.9)
  graphics::axis(1, at = x, labels = labels, tick = FALSE, cex.axis = 0.82,
                 padj = 0.6, las = 1)
  graphics::arrows(x, d$spd_mean - 1.96 * d$spd_mcse, x, d$spd_mean + 1.96 * d$spd_mcse,
                   angle = 90, code = 3, length = 0.05, col = "#222222")
  graphics::points(x, d$spd_mean, pch = 19, cex = 1.5, col = cols)

  plot(d$spd_mean, d$naive_coefficient_gap_mean, xlim = c(-1.25, 1.5), ylim = c(-0.40, 0.08),
       pch = 19, cex = 1.5, col = cols, xlab = "Mean SPD",
       ylab = "Mean outcome coefficient gap",
       main = "B. Outcome coefficient gaps are not causal bias",
       cex.main = 1.05, font.main = 1)
  graphics::grid(col = "grey88", lty = 1)
  graphics::abline(h = 0, v = 0, col = "#666666", lwd = 0.9)
  graphics::text(d$spd_mean, d$naive_coefficient_gap_mean, labels = short,
                 pos = 4, offset = 0.6, cex = 0.95)

  graphics::mtext("Diagnostic-complementarity simulation: 2,000 subjects x 500 replicates per scenario",
                  outer = TRUE, side = 3, line = 0.3, font = 2, cex = 1.05)
  graphics::mtext(paste("Error bars in A: +/-1.96 MCSE. Gap in B: treatment-only Cox coefficient",
                        "minus conditional DGM coefficient (-0.5).\n",
                        "Selection, omitted prognosis and Cox noncollapsibility all contribute to that gap."),
                  outer = TRUE, side = 1, line = 0.6, cex = 0.82)
}

if (have_summary) save_fig("Figure2_R", draw_figure2, width = 10.8, height = 4.8)
