# scripts/15_worked_example_r_native.R
# ------------------------------------------------------------
# R-only re-derivation of the synthetic worked example.
#
# The script regenerates the development and evaluation cohorts in R and refits
# the prediction models, so it does not reproduce Table 4, eTable S6 or
# eTable S7:
#   1. R's random number generator gives different cohorts from the stored ones.
#   2. The flexible score here is an unpenalized logistic regression with natural
#      splines and the X1:X2 interaction. The gradient-boosting model used for
#      the reported analysis has no exact R counterpart.
# scripts/12 recomputes the reported values from the stored cohort and
# predictions. This script is for studying the example further in R.
#
# Output: $OUT_DIR/worked_example_r_native.csv
#         $OUT_DIR/worked_example_r_native_scale.csv
#
# Usage:  Rscript scripts/15_worked_example_r_native.R
#         BOOT=2000 SEED=20260919 Rscript scripts/15_worked_example_r_native.R
# ------------------------------------------------------------

source("R/spd_core.R")
suppressPackageStartupMessages(library(splines))

SEED    <- as.integer(Sys.getenv("SEED", "20260919"))
N_DEV   <- as.integer(Sys.getenv("N_DEVELOPMENT", "6000"))
N_EVAL  <- as.integer(Sys.getenv("N_EVALUATION", "2000"))
BOOT    <- as.integer(Sys.getenv("BOOT", "500"))
OUT_DIR <- Sys.getenv("OUT_DIR", "output/r_native")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(strrep("=", 78), "\n")
cat("R-only re-derivation of the worked example. These are not the reported\n")
cat("Table 4 / eTable S6 values; see the header of this file.\n")
cat(strrep("=", 78), "\n\n")

HORIZON <- 5

draw_cohort <- function(n) {
  X <- matrix(stats::rnorm(n * 5), ncol = 5,
              dimnames = list(NULL, paste0("X", 1:5)))
  list(X = X, lp = true_lp(X))
}

reference_event <- function(lp) {
  as.integer(stats::rexp(length(lp), rate = 0.04 * exp(pmin(pmax(lp, -5), 5))) <= HORIZON)
}

fit_learners <- function(X, y) {
  df <- as.data.frame(X); df$y <- y
  linear <- stats::glm(y ~ X1 + X2 + X3 + X4 + X5, data = df, family = stats::binomial())
  # Flexible alternative available in base R: smooth terms for each predictor
  # plus the one pairwise interaction present in the generating mechanism.
  flexible <- stats::glm(y ~ ns(X1, 4) + ns(X2, 4) + ns(X3, 4) + ns(X4, 4) + ns(X5, 4) + X1:X2,
                         data = df, family = stats::binomial())
  list(linear = linear, flexible = flexible)
}

predict_prob <- function(model, X) {
  stats::predict(model, newdata = as.data.frame(X), type = "response")
}

set.seed(SEED)
dev_c  <- draw_cohort(N_DEV)
y_dev  <- reference_event(dev_c$lp)
models <- fit_learners(dev_c$X, y_dev)

eval_c <- draw_cohort(N_EVAL)
lp     <- eval_c$lp
y0     <- reference_event(lp)
a      <- stats::rbinom(N_EVAL, 1, 0.5)
t_dev  <- stats::rexp(N_EVAL, rate = 0.15 * exp(pmin(pmax(0.80 * lp + 0.60 * a, -5), 5)))
t_out  <- stats::rexp(N_EVAL, rate = 0.04 * exp(pmin(pmax(lp - 0.50 * a, -5), 5)))
time   <- pmin(t_dev, t_out, HORIZON)
dev    <- as.integer(t_dev <= t_out & t_dev <= HORIZON)
outcome <- as.integer(t_out < t_dev & t_out <= HORIZON)

mu     <- 0.2 * exp(pmin(pmax(lp, -5), 5))
known_risk <- -expm1(-mu)                          # 1 - exp(-mu), accurate for small mu
probs  <- list("Known reference risk" = known_risk,
               "Main-effects logistic" = predict_prob(models$linear, eval_c$X),
               "Splines + X1:X2"       = predict_prob(models$flexible, eval_c$X))

cat(sprintf("Development events %d/%d | evaluation: deviations %d, outcomes %d, reference events %d\n\n",
            sum(y_dev), N_DEV, sum(dev), sum(outcome), sum(y0)))

rows <- lapply(names(probs), function(nm) {
  p   <- probs[[nm]]
  s   <- sc_logit_clip(p, 1e-6)
  fit <- spd_from_score(time, dev, s, ties = "breslow")
  cal <- calibration(y0, s)
  data.frame(score = nm, auc = auc_roc(y0, p), brier = brier(y0, p),
             calibration_intercept = cal$intercept, calibration_slope = cal$slope,
             score_mean = fit$score_mean, score_sd = fit$score_sd,
             spd = fit$spd, conditional_se = fit$se, hr = fit$hr,
             stringsAsFactors = FALSE)
})
result <- do.call(rbind, rows)

# ---------------------------------------------------------------- bootstrap

if (BOOT > 0L) {
  cat(sprintf("Two-sample bootstrap with refitting: %d replicates\n", BOOT))
  draws <- array(NA_real_, dim = c(BOOT, length(probs), 3),
                 dimnames = list(NULL, names(probs), c("spd", "auc", "brier")))
  pb <- utils::txtProgressBar(min = 0, max = BOOT, style = 3)
  for (b in seq_len(BOOT)) {
    set.seed(SEED + 1000000L + b)
    di <- sample.int(N_DEV,  N_DEV,  replace = TRUE)
    ei <- sample.int(N_EVAL, N_EVAL, replace = TRUE)
    mb <- tryCatch(suppressWarnings(fit_learners(dev_c$X[di, , drop = FALSE], y_dev[di])),
                   error = function(e) NULL)
    if (is.null(mb)) next
    pb_list <- list(known_risk[ei],
                    suppressWarnings(predict_prob(mb$linear,   eval_c$X[ei, , drop = FALSE])),
                    suppressWarnings(predict_prob(mb$flexible, eval_c$X[ei, , drop = FALSE])))
    for (j in seq_along(pb_list)) {
      pj <- pb_list[[j]]
      fj <- tryCatch(spd_from_score(time[ei], dev[ei], sc_logit_clip(pj, 1e-6), "breslow"),
                     error = function(e) NULL)
      if (!is.null(fj)) draws[b, j, "spd"] <- fj$spd
      draws[b, j, "auc"]   <- auc_roc(y0[ei], pj)
      draws[b, j, "brier"] <- brier(y0[ei], pj)
    }
    utils::setTxtProgressBar(pb, b)
  }
  close(pb)
  q <- function(v) stats::quantile(v[is.finite(v)], c(0.025, 0.975), names = FALSE, type = 7)
  for (j in seq_along(probs)) {
    ci <- q(draws[, j, "spd"])
    result$bootstrap_spd_low[j]  <- ci[1]
    result$bootstrap_spd_high[j] <- ci[2]
    result$bootstrap_hr_low[j]   <- exp(ci[1])
    result$bootstrap_hr_high[j]  <- exp(ci[2])
  }
  contrasts <- do.call(rbind, lapply(c("spd", "auc", "brier"), function(m) {
    v  <- draws[, 3, m] - draws[, 2, m]
    ci <- q(v)
    data.frame(metric = m, contrast = "flexible_minus_linear",
               point = result[[if (m == "spd") "spd" else m]][3] -
                       result[[if (m == "spd") "spd" else m]][2],
               bootstrap_low = ci[1], bootstrap_high = ci[2], stringsAsFactors = FALSE)
  }))
  utils::write.csv(contrasts, file.path(OUT_DIR, "worked_example_r_native_contrasts.csv"),
                   row.names = FALSE)
  cat("\nPaired contrasts, flexible minus main-effects:\n")
  print(contrasts, row.names = FALSE, digits = 4)
}

cat("\nNative R worked example\n")
print(result, row.names = FALSE, digits = 4)
utils::write.csv(result, file.path(OUT_DIR, "worked_example_r_native.csv"), row.names = FALSE)

# ---------------------------------------------------------------- score scale

exact <- list(sc_logit_exact_known(mu),
              log(probs[[2]]) - log1p(-probs[[2]]),
              log(probs[[3]]) - log1p(-probs[[3]]))
conv <- list()
for (eps in c(1e-4, 1e-6, 1e-10)) {
  conv[[sprintf("Probability clip %g", eps)]] <- lapply(probs, sc_logit_clip, eps = eps)
}
conv[["Exact logit, no probability clipping"]] <- exact
conv[["Logit winsorized at own 0.5/99.5 percentiles"]] <- lapply(exact, sc_winsorize)
conv[["Cloglog, probability clip 1e-6"]] <- lapply(probs, sc_cloglog, eps = 1e-6)

scale_tab <- do.call(rbind, lapply(names(conv), function(nm) {
  v <- vapply(conv[[nm]], function(s) spd_from_score(time, dev, s, "breslow")$spd, numeric(1))
  data.frame(convention = nm, known_risk = v[1], linear = v[2], flexible = v[3],
             stringsAsFactors = FALSE)
}))
cat("\nScore-scale sensitivity in this R cohort\n")
print(scale_tab, row.names = FALSE, digits = 4)
utils::write.csv(scale_tab, file.path(OUT_DIR, "worked_example_r_native_scale.csv"),
                 row.names = FALSE)

cat(sprintf("\nWrote native R output to %s\n", OUT_DIR))
cat("These are not the reported values; scripts/12 recomputes those.\n")
