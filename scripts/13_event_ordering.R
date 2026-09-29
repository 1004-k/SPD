# scripts/13_event_ordering.R
# ------------------------------------------------------------
# Controlled event-ordering illustration (eTable S8, eFigure S5), in R.
#
# Two modes:
#
#   MODE = "verify"  (default)
#     Reads the stored constructed cohort and recomputes the reported values.
#     This reproduces eTable S8 exactly, because the cohort is fixed input.
#
#   MODE = "generate"
#     Builds a new constructed cohort with R's random number generator. The
#     design is the same, but the draws differ from those of the stored
#     cohort, so the numbers will not equal the reported ones.
#
# Inputs (verify mode):  $TIMING_DIR/timing_illustration_cohort.csv
# Outputs:               $OUT_DIR/timing_illustration_summary_R.csv
#                        $FIG_DIR/eFigureS5_timing_R.tiff (and .pdf)
#
# Usage:  Rscript scripts/13_event_ordering.R
#         MODE=generate SEED=20260919 Rscript scripts/13_event_ordering.R
# ------------------------------------------------------------

source("R/spd_core.R")

MODE       <- Sys.getenv("MODE", "verify")
TIMING_DIR <- Sys.getenv("TIMING_DIR", "results/additional_checks")
OUT_DIR    <- Sys.getenv("OUT_DIR", "output/checks")
FIG_DIR    <- Sys.getenv("FIG_DIR", "output/figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (MODE == "verify") {
  path <- file.path(TIMING_DIR, "timing_illustration_cohort.csv")
  if (!file.exists(path)) {
    stop(sprintf("Could not find %s. Set TIMING_DIR, or run with MODE=generate.", path),
         call. = FALSE)
  }
  co    <- utils::read.csv(path)
  score <- co$score
  dev   <- co$deviation
  early <- co$time_earlier
  later <- co$time_later
  cat("Mode: verify (reading the stored constructed cohort)\n\n")
} else {
  seed <- as.integer(Sys.getenv("SEED", "20260919"))
  n    <- as.integer(Sys.getenv("N", "2000"))
  set.seed(seed)
  score <- stats::rnorm(n)
  dev   <- stats::rbinom(n, 1, 0.5)
  ids   <- which(dev == 1L)
  # Imperfect relation between score and departure order.
  rank_order <- order(score[ids] + stats::rnorm(length(ids)))
  event_times <- sort(stats::runif(length(ids), 0.1, 4.9))
  early <- rep(5, n); later <- rep(5, n)
  early[ids[rank_order]] <- rev(event_times)   # largest noisy rank leaves first
  later[ids[rank_order]] <- event_times        # largest noisy rank leaves last
  stopifnot(identical(sort(early), sort(later)))
  cat("Mode: generate (new R random draws; values differ from the reported table)\n\n")
}

# The multiset of observed times must be identical across the two constructions.
stopifnot(isTRUE(all.equal(sort(early), sort(later))))

endpoint_smd <- smd(score[dev == 1L], score[dev == 0L])

cases <- list("Higher-score deviators earlier" = early,
              "Higher-score deviators later"   = later,
              "Earlier case, all times doubled" = 2 * early)

rows <- lapply(names(cases), function(lab) {
  tt  <- cases[[lab]]
  fit <- spd_from_score(tt, dev, score, ties = "breslow")
  man <- cox_one_manual(tt, dev, standardize(score))
  stopifnot(abs(fit$spd - man$beta) < 1e-7)
  data.frame(case = lab, n = length(score), deviations = sum(dev),
             deviation_proportion = mean(dev), endpoint_score_SMD = endpoint_smd,
             SPD = fit$spd, HR = fit$hr, stringsAsFactors = FALSE)
})
res <- do.call(rbind, rows)

# A common change of time units cannot alter a partial likelihood that depends
# only on the ordering of observed times.
stopifnot(abs(res$SPD[1] - res$SPD[3]) < 1e-12)

cat("eTable S8: identical final summaries with different score-event ordering\n")
print(res, row.names = FALSE, digits = 6)

if (MODE == "verify") {
  cat("\n  comparison with the reported values\n")
  cat(sprintf("  %-46s %10s %10s %10s  %s\n", "quantity", "reported", "R value", "abs diff", ""))
  check_row("deviations",              1006,     res$deviations[1],           0)
  check_row("deviation proportion",    0.503,    res$deviation_proportion[1], 5e-4)
  check_row("final score SMD",        -0.058,    res$endpoint_score_SMD[1],   5e-4)
  check_row("SPD, higher-score earlier", 0.104,  res$SPD[1],                  5e-4)
  check_row("SPD, higher-score later",  -0.188,  res$SPD[2],                  5e-4)
  check_row("SPD, doubled times",        0.104,  res$SPD[3],                  5e-4)
  stored_path <- file.path(TIMING_DIR, "timing_illustration_summary.csv")
  if (file.exists(stored_path)) {
    st <- utils::read.csv(stored_path)
    st <- st[match(res$case, st$case), ]
    d <- max(abs(c(st$endpoint_score_SMD - res$endpoint_score_SMD, st$SPD - res$SPD, st$HR - res$HR)))
    cat(sprintf("  largest absolute difference from the stored full-precision values: %.1e\n", d))
  }
}

utils::write.csv(res, file.path(OUT_DIR, "timing_illustration_summary_R.csv"), row.names = FALSE)

# ---------------------------------------------------------------- eFigure S5

draw_figure <- function() {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.4, 3.0, 0.8),
                      mgp = c(2.5, 0.8, 0), las = 1)
  on.exit(graphics::par(op), add = TRUE)
  grid_t <- seq(0, 5, length.out = 201)
  high   <- score >= stats::median(score)
  panels <- list(list(t = early, spd = res$SPD[1], title = "A. Higher-score deviators earlier"),
                 list(t = later, spd = res$SPD[2], title = "B. Higher-score deviators later"))
  for (k in seq_along(panels)) {
    pk <- panels[[k]]
    plot(NA, xlim = c(0, 5), ylim = c(0, 0.65), xlab = "Follow-up time",
         ylab = if (k == 1) "Cumulative proportion with deviation" else "",
         main = pk$title, cex.main = 1.0, font.main = 1)
    graphics::grid(col = "grey88", lty = 1)
    for (g in list(list(m = high,  col = "#24658c", lab = "Score at/above median"),
                   list(m = !high, col = "#bd5c36", lab = "Score below median"))) {
      y <- vapply(grid_t, function(u) mean(dev[g$m] == 1L & pk$t[g$m] <= u), numeric(1))
      graphics::lines(grid_t, y, type = "s", col = g$col, lwd = 2.2)
    }
    graphics::text(0.15, 0.62, adj = c(0, 1), cex = 0.92,
                   labels = sprintf("SPD = %.3f\nFinal score SMD = %.3f",
                                    pk$spd, endpoint_smd))
    if (k == 2) {
      graphics::legend("bottomright", bty = "n", cex = 0.88, lwd = 2.2,
                       col = c("#24658c", "#bd5c36"),
                       legend = c("Score at/above median", "Score below median"))
    }
  }
  graphics::mtext("Same deviators and event-time distribution; different score-event ordering",
                  outer = TRUE, line = -1.4, cex = 1.05)
}

tiff_path <- file.path(FIG_DIR, "eFigureS5_timing_R.tiff")
grDevices::tiff(tiff_path, width = 10, height = 4, units = "in", res = 600,
                compression = "lzw", type = "cairo")
draw_figure(); invisible(grDevices::dev.off())

grDevices::pdf(file.path(FIG_DIR, "eFigureS5_timing_R.pdf"), width = 10, height = 4)
draw_figure(); invisible(grDevices::dev.off())

cat(sprintf("\nWrote %s and the PDF companion\n", tiff_path))
cat(sprintf("Wrote %s\n", file.path(OUT_DIR, "timing_illustration_summary_R.csv")))
