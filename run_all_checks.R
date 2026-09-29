# run_all_checks.R
# ------------------------------------------------------------
# Recomputes the reported results that come from stored data (a few seconds):
#   scripts/12  Table 4, eTables S6 and S7, paired bootstrap contrasts
#   scripts/13  eTable S8 and eFigure S5
#   scripts/16  Table 3 and eTable S5
#   scripts/14  Figures 1 and 2
# Set NATIVE=1 to also run scripts/15 (R-only re-derivation of the worked
# example; its numbers differ from the reported ones).
#
# Run from the repository root:
#   Rscript run_all_checks.R
#
# Results are written to output/checks and figures to output/figures.
# ------------------------------------------------------------

stopifnot(file.exists("R/spd_core.R"))

run_step <- function(path, label) {
  cat("\n", strrep("-", 78), "\n", label, "\n", strrep("-", 78), "\n", sep = "")
  ok <- tryCatch({ source(path, echo = FALSE, local = new.env()); TRUE },
                 error = function(e) { cat("FAILED: ", conditionMessage(e), "\n", sep = ""); FALSE })
  invisible(ok)
}

results <- c(
  "12 worked example" = run_step("scripts/12_verify_worked_example.R",
                                 "12  Table 4, eTables S6 and S7, bootstrap intervals"),
  "13 event ordering" = run_step("scripts/13_event_ordering.R",
                                 "13  eTable S8 and eFigure S5"),
  "16 complementarity" = run_step("scripts/16_verify_complementarity.R",
                                  "16  Table 3 and eTable S5"),
  "14 figures"        = run_step("scripts/14_make_figures.R",
                                 "14  Figures 1 and 2")
)

if (identical(Sys.getenv("NATIVE"), "1")) {
  results["15 R-only worked example"] <-
    run_step("scripts/15_worked_example_r_native.R",
             "15  R-only re-derivation of the worked example")
}

cat("\n", strrep("=", 78), "\n", sep = "")
for (nm in names(results)) cat(sprintf("  %-34s %s\n", nm, if (results[[nm]]) "completed" else "FAILED"))
cat(strrep("=", 78), "\n")
if (!all(results)) quit(save = "no", status = 1L)
