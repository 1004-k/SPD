# R/spd_core.R
# ------------------------------------------------------------
# Core helpers for scalar standardized prognostic dependence (SPD).
#
# Dependencies: base R and survival.
#
# Provided:
#   spd_from_score()      SPD for a score vector (standardize -> cause-specific Cox)
#   cox_one_manual()      Newton fit with step-halving (Breslow partial likelihood)
#   smd()                 standardized mean difference, pooled-variance denominator
#   km_before()           Kaplan-Meier survival immediately before each observed time
#   ress_exit()           descriptive exit-time weight concentration (rESSexit, %)
#   auc_roc(), brier()    prediction metrics without extra packages
#   calibration()         intercept and slope of an unstandardized logit score
#   score conventions:    sc_logit_clip(), sc_logit_exact_known(), sc_winsorize(),
#                         sc_probability(), sc_cloglog()
#
# Convention notes
#   * SPD is the Cox coefficient for a score standardized in the analysis sample
#     (sample SD, denominator n - 1), so it is the log hazard ratio per 1 SD.
#   * Breslow tie handling is the default here, as in the worked example and
#     the diagnostic-complementarity simulation. The estimator simulation in
#     scripts/07 uses Efron; pass ties = "efron" to match it.
# ------------------------------------------------------------

suppressPackageStartupMessages(library(survival))

standardize <- function(x) {
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) stop("standardize(): non-finite or zero SD", call. = FALSE)
  (x - mean(x)) / s
}

# ---------------------------------------------------------------- SPD

spd_from_score <- function(time, event, score, ties = c("breslow", "efron")) {
  ties <- match.arg(ties)
  z <- standardize(score)
  fit <- survival::coxph(survival::Surv(time, event) ~ z, ties = ties)
  beta <- as.numeric(stats::coef(fit)[1])
  se   <- sqrt(as.numeric(stats::vcov(fit)[1, 1]))
  list(spd = beta, se = se, hr = exp(beta),
       wald_low = exp(beta - 1.96 * se), wald_high = exp(beta + 1.96 * se),
       score_mean = mean(score), score_sd = stats::sd(score))
}

# Independent one-covariate fit: Newton steps with step-halving whenever a
# proposed update lowers the Breslow partial log likelihood. Used to check
# coxph() rather than to replace it.
cox_one_manual <- function(time, event, x, tol = 1e-10, max_iter = 100L) {
  o  <- order(time, method = "radix")
  t  <- time[o]; d <- event[o]; xx <- x[o]
  ev <- which(d == 1L)
  if (length(ev) < 2L) stop("cox_one_manual(): fewer than two events", call. = FALSE)
  # First index whose time is >= this event time, i.e. the start of the Breslow
  # risk set. left.open = TRUE counts strictly smaller times, so +1 lands on the
  # first tied observation and every tie shares one risk set.
  risk_start <- findInterval(t[ev], t, left.open = TRUE) + 1L

  score_info <- function(beta) {
    eta <- beta * xx
    r   <- exp(eta - max(eta))
    s0  <- rev(cumsum(rev(r)))[risk_start]
    s1  <- rev(cumsum(rev(r * xx)))[risk_start]
    s2  <- rev(cumsum(rev(r * xx * xx)))[risk_start]
    mu  <- s1 / s0
    list(u    = sum(xx[ev] - mu),
         info = sum(s2 / s0 - mu^2),
         ll   = sum(eta[ev] - log(s0) - max(eta)))
  }

  beta <- 0
  for (i in seq_len(max_iter)) {
    st <- score_info(beta)
    if (!is.finite(st$info) || st$info <= 0) stop("cox_one_manual(): non-positive information", call. = FALSE)
    step <- st$u / st$info
    frac <- 1
    while (frac > 1e-8 && score_info(beta + frac * step)$ll < st$ll - 1e-9) frac <- frac / 2
    beta <- beta + frac * step
    if (abs(step) < tol) {
      return(list(beta = beta, se = 1 / sqrt(score_info(beta)$info), iterations = i))
    }
  }
  stop("cox_one_manual(): did not converge", call. = FALSE)
}

# ---------------------------------------------------------------- diagnostics

smd <- function(a, b) {
  den <- sqrt((stats::var(a) + stats::var(b)) / 2)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  (mean(a) - mean(b)) / den
}

km_before <- function(time, event) {
  o  <- order(time, method = "radix")
  tt <- time[o]; ee <- event[o]
  n  <- length(tt)
  ans <- numeric(n); surv <- 1; i <- 1L
  while (i <= n) {
    j <- i
    while (j < n && tt[j + 1L] == tt[i]) j <- j + 1L
    ans[i:j] <- surv                                  # value strictly before this time
    dd <- sum(ee[i:j] == 1L, na.rm = TRUE)
    at_risk <- n - i + 1L
    if (dd > 0 && at_risk > 0) surv <- surv * max(0, 1 - dd / at_risk)
    i <- j + 1L
  }
  out <- numeric(n); out[o] <- ans; out
}

# Descriptive cross-subject exit-time weight concentration, as a percentage of N.
# Not conditional IPCW, not a risk-set-specific effective sample size, and not
# evidence about positivity.
ress_exit <- function(time, event, floor = 1e-8) {
  w <- 1 / pmax(km_before(time, event), floor)
  100 * (sum(w)^2 / sum(w^2)) / length(w)
}

# ---------------------------------------------------------------- metrics

# Mann-Whitney form of the ROC area; ties contribute one half.
auc_roc <- function(y, p) {
  y <- as.integer(y)
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

brier <- function(y, p) mean((p - as.integer(y))^2)

# Intercept and slope from a logistic regression of the reference outcome on the
# UNSTANDARDIZED logit score. Ideal values are 0 and 1.
calibration <- function(y, logit_score) {
  fit <- stats::glm(y ~ logit_score, family = stats::binomial())
  cf <- as.numeric(stats::coef(fit))
  list(intercept = cf[1], slope = cf[2])
}

# ---------------------------------------------------------------- score scales

sc_logit_clip <- function(p, eps = 1e-6) {
  q <- pmin(pmax(p, eps), 1 - eps)
  log(q) - log1p(-q)
}

# Exact logit of the known reference risk p0 = 1 - exp(-mu), evaluated as
# mu + log(1 - exp(-mu)) to avoid cancellation as the risk approaches one.
sc_logit_exact_known <- function(mu) mu + log(-expm1(-mu))

sc_winsorize <- function(s, probs = c(0.005, 0.995)) {
  q <- stats::quantile(s, probs = probs, names = FALSE, type = 7)
  pmin(pmax(s, q[1]), q[2])
}

sc_probability <- function(p) p

sc_cloglog <- function(p, eps = 1e-6) {
  q <- pmin(pmax(p, eps), 1 - eps)
  log(-log1p(-q))
}

# Latent predictor of the worked-example data-generating mechanism.
true_lp <- function(X) {
  (0.55 * X[, 1] + 0.45 * X[, 2] + 0.65 * (X[, 3]^2 - 1) +
     0.60 * X[, 1] * X[, 2] + 0.45 * sin(X[, 4]) + 0.25 * X[, 5]) / 1.55
}

# ---------------------------------------------------------------- reporting

check_row <- function(label, reported, obtained, tol) {
  diff <- abs(reported - obtained)
  cat(sprintf("  %-46s %10.4f %10.4f %10.2e  %s\n",
              label, reported, obtained, diff,
              if (is.finite(diff) && diff <= tol) "ok" else "CHECK"))
  invisible(is.finite(diff) && diff <= tol)
}
