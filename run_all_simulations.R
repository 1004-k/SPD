# run_all_simulations.R
# ------------------------------------------------------------
# Re-runs the simulation studies. This takes time: about 40 minutes on two
# cores for everything below.
#   scripts/07  estimator performance: baseline (Table 2, eTable S1),
#               faster deviation (eTable S2), more primary outcomes (eTable S3)
#   scripts/09  model-based versus robust SE, worst case (eTable S4)
#   scripts/08  eFigures S1-S4 from the output of 07 and 09
#   scripts/11  diagnostic-complementarity simulation (Table 3, eTable S5)
#
# Run from the repository root:
#   Rscript run_all_simulations.R
# The number of parallel workers for 07 and 09 is set with N_CORES (default 3).
# Output is written to output/.
# ------------------------------------------------------------

rscript <- file.path(R.home("bin"), "Rscript")

run <- function(script, env = character()) {
  cat("\n", strrep("-", 78), "\n", script, if (length(env)) paste0("  [", paste(names(env), env, sep = "=", collapse = " "), "]"),
      "\n", strrep("-", 78), "\n", sep = "")
  old <- if (length(env)) Sys.getenv(names(env), unset = NA, names = TRUE) else character()
  if (length(env)) do.call(Sys.setenv, as.list(env))
  on.exit({
    for (nm in names(old)) {
      if (is.na(old[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, stats::setNames(list(old[[nm]]), nm))
    }
  })
  status <- system2(rscript, shQuote(script))
  if (!identical(as.integer(status), 0L)) stop(sprintf("%s failed (exit status %s)", script, status), call. = FALSE)
  invisible(TRUE)
}

sim07 <- c(B = "500", N_REF = "50000", B_REF = "3", SEED = "2026")

run("scripts/07_run_sim_spd_scalar.R", c(sim07, OUT_TAG = "baseline"))
run("scripts/07_run_sim_spd_scalar.R", c(sim07, LAMBDA_DEV = "0.50", OUT_TAG = "fastdev"))
run("scripts/07_run_sim_spd_scalar.R", c(sim07, LAMBDA_EVENT = "0.20", OUT_TAG = "fastevent"))
run("scripts/09_run_sim_spd_scalar_robust_se_worstcase.R", c(B = "500", SEED = "2026"))
run("scripts/08_make_simulation_figures.R")
run("scripts/11_run_sim_spd_diagnostic_complementarity.R", c(N = "2000", B = "500", SEED = "20260918"))

cat("\nAll simulations finished. Output is in output/.\n")
