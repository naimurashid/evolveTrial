# Tests for the exact posterior probabilities of the phase 2 engine (step P2-1, src/posterior_exact.cpp). Fast checks
# only; the full gate is dev/perf/phase2/g1_check.R.

ns <- function(name) getFromNamespace(name, "evolveTrial")

ref_two <- function(sh, sc, c = log(2)) {
  # Pr(sc1 G1 + sc2 G2 > c) = Q1(c / sc1) + int_0^{c / sc1} f1(g) Q2((c - sc1 g) / sc2) dg, split at quantiles of both
  # terms; for a shape below 1 each piece is integrated in v = g^a1, where the integrand is smooth.
  x <- c / sc[1]
  h <- function(g) pgamma((c - sc[1] * g) / sc[2], sh[2], lower.tail = FALSE)
  qs <- qgamma(c(1e-6, 1e-3, 0.05, 0.25, 0.5, 0.75, 0.95, 0.999, 1 - 1e-6), sh[1])
  qj <- (c - sc[2] * qgamma(c(1e-6, 1e-3, 0.05, 0.5, 0.95, 0.999, 1 - 1e-6), sh[2])) / sc[1]
  br <- sort(unique(c(0, qs[qs < x], qj[qj > 0 & qj < x], x)))
  tot <- pgamma(x, sh[1], lower.tail = FALSE)
  for (p in seq_len(length(br) - 1L)) {
    lo <- br[p]; hi <- br[p + 1L]
    if (sh[1] < 1) {
      a1 <- sh[1]
      f <- function(v) { g <- v^(1 / a1); exp(-g - lgamma(a1 + 1)) * h(g) }
      tot <- tot + stats::integrate(f, lo^a1, hi^a1, rel.tol = 1e-11, abs.tol = 1e-15, subdivisions = 5000L,
                                  stop.on.error = FALSE)$value
    } else {
      f <- function(g) dgamma(g, sh[1]) * h(g)
      tot <- tot + stats::integrate(f, lo, hi, rel.tol = 1e-11, abs.tol = 1e-15, subdivisions = 5000L,
                                  stop.on.error = FALSE)$value
    }
  }
  tot
}

test_that("one positive weight reduces to pgamma", {
  pm <- ns(".exact_pr_median")
  E <- c(3L, 1L, 0L); PT <- c(20, 8, 2); a <- rep(0.3, 3); b <- rep(0.3, 3); lens <- c(3, 3, 3)
  thr <- 2.2
  got <- pm(E, PT, a, b, lens, thr)
  expect_equal(got$method, 1L)
  expect_identical(got$below, pgamma(log(2), 3.3, scale = thr / 20.3, lower.tail = FALSE))
  expect_identical(got$above, 1 - got$below)
})

test_that("two positive weights agree with integrate (weak priors, design-like inputs)", {
  pm <- ns(".exact_pr_median")
  for (E in list(c(0L, 0L, 0L), c(1L, 0L, 0L), c(2L, 1L, 0L), c(7L, 3L, 1L))) {
    a <- rep(0.2167, 3); b <- a; PT <- c(25, 12, 3) + 5 * E; lens <- c(3, 3, 3)
    for (thr in c(3.8, 4.527, 5.999)) {
      got <- pm(E, PT, a, b, lens, thr)
      sh <- a[1:2] + E[1:2]; sc <- c(3, thr - 3) / (b[1:2] + PT[1:2])
      expect_lt(abs(got$below - ref_two(sh, sc)), 1e-9)
      expect_true(got$converged)
    }
  }
})

test_that("more weights agree with the phase 1 sampler (Monte Carlo, 4 standard errors)", {
  pm <- ns(".exact_pr_median")
  set.seed(11)
  E <- c(4L, 0L, 2L, 1L, 0L); PT <- c(30, 20, 11, 6, 2); a <- rep(0.4, 5); b <- a; lens <- rep(2, 5)
  thr <- 8.7
  got <- pm(E, PT, a, b, lens, thr)
  lam <- evolveTrial::draw_posterior_hazard_samples_cpp(5L, E, PT, a, b, 200000L)
  mc <- mean(evolveTrial::calculate_median_survival_matrix_cpp(lam, lens) < thr)
  expect_lt(abs(mc - got$below) / sqrt(got$below * (1 - got$below) / 2e5), 4)
})

test_that("edge cases follow the phase 1 tally semantics", {
  pm <- ns(".exact_pr_median")
  E <- c(2L, 1L); PT <- c(10, 5); a <- c(0.3, 0.3); b <- a; lens <- c(3, 3)
  # thr <= 0: no hazard below thr, med > thr surely
  for (thr in c(0, -1, -Inf)) {
    g <- pm(E, PT, a, b, lens, thr)
    expect_identical(c(g$below, g$above), c(0, 1))
  }
  # NaN threshold: every comparison is FALSE, both tallies are 0
  g <- pm(E, PT, a, b, lens, NaN)
  expect_identical(c(g$below, g$above), c(0, 0))
  # +Inf threshold: last hazard positive almost surely, so the median is finite
  g <- pm(E, PT, a, b, lens, Inf)
  expect_identical(c(g$below, g$above), c(1, 0))
  # zero shape in the last interval (alpha 0, no events): lambda_last = 0, cumulative hazard capped at the grid
  g <- pm(c(2L, 0L), PT, c(0.3, 0), b, lens, Inf)
  expect_identical(g$below, pgamma(log(2), 2.3, scale = 3 / 10.3, lower.tail = FALSE))
  # all positive-weight shapes zero: S = 0
  g <- pm(c(0L, 0L), PT, c(0, 0), b, lens, 4)
  expect_identical(c(g$below, g$above), c(0, 1))
  # non-positive posterior rate is replaced by 1e-6, as draw_posterior_hazard_samples_cpp does
  g <- pm(E, c(-0.3, 5), a, b, lens, 2)
  expect_identical(g$status, 1L)
  expect_identical(g$below, pgamma(log(2), 2.3, scale = 2 / 1e-6, lower.tail = FALSE))
  # invalid input: NaN both
  g <- pm(E, c(NaN, 5), a, b, lens, 2)
  expect_identical(g$status, 2L)
  expect_true(is.nan(g$below) && is.nan(g$above))
  g <- pm(c(-1L, 0L), PT, c(0.3, 0.3), b, lens, 2)
  expect_true(is.nan(g$below))
  # vectorized hooks agree with the list hook
  thr <- c(1, 3, 4.5, 30)
  expect_identical(ns(".exact_pr_median_below")(E, PT, a, b, lens, thr),
                   vapply(thr, function(t) pm(E, PT, a, b, lens, t)$below, 0))
  expect_identical(ns(".exact_pr_median_above")(E, PT, a, b, lens, thr),
                   vapply(thr, function(t) pm(E, PT, a, b, lens, t)$above, 0))
})

test_that("equal scales merge into one gamma; extreme scale ratios stay accurate", {
  gst <- ns(".exact_gamma_sum_tail")
  g <- gst(c(0.4, 1.1), c(0.5, 0.5), log(2))
  expect_equal(g$method, 1L)
  expect_identical(g$tail, pgamma(log(2), 1.5, scale = 0.5, lower.tail = FALSE))
  for (sc in list(c(0.12, 0.00093), c(5, 1e-5), c(1e-3, 2e-3), c(40, 0.01))) {
    sh <- c(0.318, 1.07)
    expect_lt(abs(gst(sh, sc, log(2))$tail - ref_two(sh, sc)), 1e-9)
  }
  # large shapes (300 events)
  expect_lt(abs(gst(c(300.2, 120.2), c(0.0015, 0.0021), log(2))$tail -
                  ref_two(c(300.2, 120.2), c(0.0015, 0.0021))), 1e-9)
})

test_that("between-arm probabilities are pnorm on the Laplace output, with rnorm edge semantics", {
  plr <- ns(".exact_pr_loghr")
  K <- 9L
  EC <- c(5L, 4L, 3L, 2L, 1L, 1L, 0L, 0L, 0L); ET <- c(3L, 3L, 2L, 1L, 1L, 0L, 0L, 0L, 0L)
  PTC <- c(60, 45, 30, 20, 12, 8, 4, 2, 1); PTT <- PTC * 1.1; a <- rep(0.48, K); lens <- rep(24 / K, K)
  m <- c(log(0.75), 0, log(1.2))
  got <- plr(EC, PTC, ET, PTT, a, a, lens, 0, 1, 1e-6, 50L, m)
  bp <- evolveTrial::ph_beta_mode_var_cpp(EC, PTC, ET, PTT, a, a, 0, 1, 1e-6, 50L)
  expect_identical(c(got$mean, got$sd), c(bp$mean, bp$sd))
  expect_identical(got$pr_ge, pnorm(m, bp$mean, bp$sd, lower.tail = FALSE))
  expect_identical(got$pr_lt, pnorm(m, bp$mean, bp$sd))
  set.seed(5)
  o <- ns(".chunk_ph_sample")(EC, PTC, ET, PTT, a, a, lens, 0, 1, 1e-6, 50L, 100000L)
  for (i in seq_along(m)) {
    p <- got$pr_ge[i]
    expect_lt(abs(mean(o$logHR >= m[i]) - p) / sqrt(p * (1 - p) / 1e5), 4)
  }
  lp <- ns(".exact_loghr_probs")
  expect_true(all(is.nan(lp(NaN, 1, 0))))
  expect_true(all(is.nan(lp(0, NaN, 0))))
  expect_true(all(is.nan(lp(0, -1, 0))))
  expect_identical(unname(lp(0.2, 0, c(0.1, 0.2, 0.3))), cbind(c(1, 1, 0), c(0, 0, 1)))
  expect_identical(unname(lp(Inf, 1, 5)[1, ]), c(1, 0))
})
