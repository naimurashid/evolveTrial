# Unit tests for the R-semantics floating-point helpers in src/r_semantics.h (engine C++ port, step 2b).
# Each helper must reproduce the matching R function bit for bit (identical()), including R's long double
# accumulation in sum() and mean(). Inputs are seeded.

rs <- function(name) getFromNamespace(name, "evolveTrial")

# identical() over a list of cases, reporting how many cases differ and the first one
expect_all_identical <- function(got, want, what) {
  ok <- vapply(seq_along(want), function(i) identical(got[[i]], want[[i]]), logical(1))
  nbad <- sum(!ok)
  info <- if (nbad > 0L) {
    i <- which(!ok)[1L]
    sprintf("%s: %d of %d cases differ; first is case %d: got %s, want %s", what, nbad, length(ok), i,
            format(got[[i]], digits = 17), format(want[[i]], digits = 17))
  } else sprintf("%s: %d cases", what, length(ok))
  testthat::expect(nbad == 0L, info)
  invisible(length(ok))
}

# ---- seeded input generator ------------------------------------------------------------------------------------------
gen_inputs <- function(n_cases, seed, max_len = 300L, min_len = 0L) {
  set.seed(seed)
  fixed_lens <- c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L)
  out <- vector("list", n_cases)
  for (i in seq_len(n_cases)) {
    len <- if (i <= length(fixed_lens) * 8L) fixed_lens[(i - 1L) %% length(fixed_lens) + 1L] else
      sample(min_len:max_len, 1L)
    fam <- (i - 1L) %% 9L
    out[[i]] <- switch(fam + 1L,
      rexp(len, rate = 1 / runif(1, 2, 12)),                                        # exponential
      rweibull(len, shape = 1.3, scale = runif(1, 3.8, 9) / log(2)^(1 / 1.3)),      # Weibull like the engine's event times
      1e8 + runif(len, 0, 1e-3),                                                     # 1e8 scale with small offsets
      round(runif(len, 0, 30), 1),                                                   # decimal fractions (multiples of 0.1)
      sample(c(0.1, 0.2, 0.3, 2.5, 7), len, replace = TRUE),                         # heavy ties
      sort(rweibull(len, 1.3, 6)),                                                   # sorted
      sort(rweibull(len, 1.3, 6), decreasing = TRUE),                                # reverse sorted
      rnorm(len, mean = 0, sd = 10^runif(1, -3, 6)),                                 # mixed sign, wide scale
      { v <- rep(1e-16, len); if (len > 0L) v[sample(len, 1L)] <- runif(1, 0.5, 2); v }   # one large value among tiny ones
    )
  }
  out
}

test_that("long double accumulation: a double accumulator would give a different answer", {
  x <- c(1, rep(1e-16, 10))
  expect_false(identical(sum(x), Reduce(`+`, x)))        # precondition: the example separates the two
  expect_identical(rs("rsem_sum")(x), sum(x))
  expect_identical(rs("rsem_mean")(x), mean(x))
  y <- c(1e8, rep(1e-8, 1000), -1e8)
  expect_identical(rs("rsem_sum")(y), sum(y))
  expect_identical(rs("rsem_mean")(y), mean(y))
  if (.Machine$sizeof.longdouble == 16L) expect_identical(rs("rsem_ld_digits")(), 64L)  # x87 extended on x86-64
})

test_that("rsem_sum equals sum() on 2,000 seeded inputs plus a 1e5 vector and NA cases", {
  cases <- gen_inputs(2000L, 11L)
  cases[[length(cases) + 1L]] <- rweibull(1e5, 1.3, 7)
  cases[[length(cases) + 1L]] <- 1e8 + runif(1e5, 0, 1e-3)
  expect_all_identical(lapply(cases, rs("rsem_sum")), lapply(cases, sum), "sum")
  na_cases <- list(c(NA, 1, 2), c(1, NaN, 3), c(Inf, -Inf), c(Inf, 1), c(-Inf, 5), 1.7e308 + c(0, 0.1e308), c(1.7e308, 1.7e308), NA_real_)
  expect_all_identical(lapply(na_cases, rs("rsem_sum")), lapply(na_cases, sum), "sum special values")
})

test_that("rsem_mean equals mean() on 2,000 seeded inputs plus a 1e5 vector and special values", {
  cases <- gen_inputs(2000L, 12L)
  cases[[length(cases) + 1L]] <- rweibull(1e5, 1.3, 7)
  cases[[length(cases) + 1L]] <- 1e8 + runif(1e5, 0, 1e-3)
  expect_all_identical(lapply(cases, rs("rsem_mean")), lapply(cases, mean), "mean")
  sp <- list(c(NA, 1, 2), c(1, NaN, 3), c(Inf, -Inf), c(Inf, 1), c(1.7e308, 1.7e308), c(-1.7e308, -1.7e308, 1e300),
             numeric(0), 5, c(1.7e308, 1.7e308, -1.7e308))
  expect_all_identical(lapply(sp, rs("rsem_mean")), lapply(sp, mean), "mean special values")
})

test_that("rsem_median equals stats::median on 2,000 seeded inputs", {
  cases <- gen_inputs(2000L, 13L)
  cases[[length(cases) + 1L]] <- rweibull(1e5, 1.3, 7)
  cases[[length(cases) + 1L]] <- rweibull(1e5 + 1L, 1.3, 7)
  expect_all_identical(lapply(cases, rs("rsem_median")), lapply(cases, stats::median), "median")
  sp <- list(numeric(0), c(NA, 1, 2), c(1, NaN), 4, c(2, 1), c(1e8 + 1e-3, 1e8), c(0.1, 0.2), c(0.1, 0.7, 0.3, 0.2))
  expect_all_identical(lapply(sp, rs("rsem_median")), lapply(sp, stats::median), "median special values")
})

test_that("rsem_quantile7 equals stats::quantile(type = 7) on 2,000 seeded inputs", {
  cases <- gen_inputs(2000L, 14L, min_len = 1L)
  set.seed(15L)
  ps <- c(0, 1, 0.5, 0.25, 0.75, 1 / 3, 2 / 3, 0.9, 1 / 7, 0.123456789)
  pcase <- lapply(seq_along(cases), function(i) if (i %% 2L == 0L) sample(ps, 1L) else runif(1L))
  got <- Map(function(x, p) rs("rsem_quantile7")(x, p), cases, pcase)
  want <- Map(function(x, p) as.numeric(stats::quantile(x, probs = p, type = 7, names = FALSE)), cases, pcase)
  expect_all_identical(got, want, "quantile7")
  # the whole probability vector the engine uses, K = 3..16, on engine-like event times
  set.seed(16L)
  got2 <- list(); want2 <- list()
  for (i in 1:300) {
    K <- sample(3:16, 1L); x <- pmin(rweibull(sample(50:300, 1L), 1.3, runif(1, 5, 14)), 120)
    pr <- seq(0, 1, length.out = K + 1)
    got2[[i]] <- vapply(pr, function(p) rs("rsem_quantile7")(x, p), numeric(1))
    want2[[i]] <- as.numeric(stats::quantile(x, probs = pr, na.rm = TRUE, type = 7))
  }
  expect_all_identical(got2, want2, "quantile7 over seq probs")
  expect_identical(rs("rsem_quantile7")(numeric(0), 0.5), NA_real_)
})

test_that("rsem_find_interval equals findInterval for every flag combination", {
  set.seed(17L)
  n_cases <- 2000L
  flags <- expand.grid(rc = c(FALSE, TRUE), ai = c(FALSE, TRUE), lo = c(FALSE, TRUE))
  got <- vector("list", n_cases); want <- vector("list", n_cases)
  for (i in seq_len(n_cases)) {
    nb <- sample(c(0L, 1L, 2L, 3L, 5L, 9L, 13L, 17L), 1L)
    br <- sort(if (i %% 4L == 0L) sample(c(0, 3, 6, 6, 9, 12), nb, replace = TRUE) else round(runif(nb, 0, 24), sample(0:3, 1L)))
    x <- c(br, br + 1e-12, br - 1e-12, -1, 0, 24, 25, 1e9, runif(5, -2, 30), if (i %% 10L == 0L) NaN)
    f <- flags[(i - 1L) %% nrow(flags) + 1L, ]
    got[[i]] <- rs("rsem_find_interval")(x, br, f$rc, f$ai, f$lo)
    want[[i]] <- findInterval(x, br, rightmost.closed = f$rc, all.inside = f$ai, left.open = f$lo)
  }
  expect_all_identical(got, want, "findInterval")
  # the flag values the engine uses (state_management.R:27-31): left.open = FALSE, rightmost.closed = FALSE, all.inside default
  br <- c(0, 3, 6, 9, 12, 15, 18, 21, 24); x <- c(0, 3, 24, 30, 11.999999999999998, 1e-300)
  expect_identical(rs("rsem_find_interval")(x, br, FALSE, FALSE, FALSE), findInterval(x, br, left.open = FALSE, rightmost.closed = FALSE))
})

test_that("rsem_mean_logical equals mean() of a logical vector for all (count, n), n up to 2,000", {
  mism <- 0L; first <- NULL; total <- 0
  for (n in 1:2000) {
    v <- logical(n)
    want <- numeric(n + 1L)
    want[1L] <- mean(v)
    for (cn in seq_len(n)) { v[cn] <- TRUE; want[cn + 1L] <- mean(v) }
    got <- rs("rsem_mean_logical")(0:n, rep(n, n + 1L))
    bad <- which(!(got == want | (is.nan(got) & is.nan(want))))   # exact double equality
    total <- total + n + 1
    if (length(bad)) { mism <- mism + length(bad); if (is.null(first)) first <- c(bad[1L] - 1L, n) }
  }
  expect(mism == 0L, sprintf("mean_logical: %d mismatches in %.0f cases; first (count, n) = %s", mism, total,
                             paste(first, collapse = ", ")))
  expect_identical(total, 2003000)
  expect_identical(rs("rsem_mean_logical")(0L, 0L), mean(logical(0)))
})

test_that("rsem_mean_logical also equals mean() where plain count / n would not (n of a few thousand, as the posterior draws)", {
  # For n up to 2,000 count / n and mean() never differ; they first differ at n = 3,511 (count 2,613). The engine takes
  # mean() over the posterior draws (hundreds to 10,000), so cover those sizes with every count.
  ns <- c(800L, 3511L, 5000L, 5998L, 7107L, 8846L, 10000L)
  n_plain_differs <- 0L
  for (n in ns) {
    v <- logical(n); want <- numeric(n + 1L); want[1L] <- mean(v)
    for (cn in seq_len(n)) { v[cn] <- TRUE; want[cn + 1L] <- mean(v) }
    got <- rs("rsem_mean_logical")(0:n, rep(n, n + 1L))
    expect(all(got == want), sprintf("mean_logical n = %d: %d mismatches of %d", n, sum(got != want), n + 1L))
    n_plain_differs <- n_plain_differs + sum(want != (0:n) / n)
  }
  expect_gt(n_plain_differs, 0L)   # the cases really separate long-double division from plain division
})

test_that("rsem_seq_probs equals seq(0, 1, length.out = K + 1) for K = 1..40", {
  got <- lapply(1:40, function(K) rs("rsem_seq_probs")(K))
  expect_all_identical(got, lapply(1:40, function(K) seq(0, 1, length.out = K + 1)), "seq_probs (K + 1 double)")
  expect_all_identical(got, lapply(1:40, function(K) seq(0, 1, length.out = K + 1L)), "seq_probs (K + 1L integer)")
  # and the observation recorded in the inventory: i * (1 / K), not i / K
  expect_identical(rs("rsem_seq_probs")(7L)[4:5], (3:4) * (1 / 7))
})

test_that("rsem_rebalance_cutpoints equals rebalance_cutpoints_from_events", {
  ref <- getFromNamespace("rebalance_cutpoints_from_events", "evolveTrial")
  set.seed(18L)
  n_cases <- 2000L
  got <- vector("list", n_cases); want <- vector("list", n_cases)
  for (i in seq_len(n_cases)) {
    n <- sample(50:300, 1L); K <- sample(3:16, 1L)
    med <- runif(1, 3.8, 9)
    x <- rweibull(n, 1.3, med / log(2)^(1 / 1.3))
    if (i %% 5L == 0L) x[sample(n, max(1L, n %/% 10L))] <- 130 + runif(max(1L, n %/% 10L), 0, 50)  # beyond the cap
    if (i %% 7L == 0L) x <- round(x, 1)                                                              # ties
    if (i %% 11L == 0L) x[sample(n, 3L)] <- 0                                                         # values at 0
    if (i %% 13L == 0L) x[sample(n, 4L)] <- 120                                                       # values at the cap
    mfu <- if (i %% 3L == 0L) 24 else 120
    got[[i]] <- rs("rsem_rebalance_cutpoints")(x, mfu, K)
    want[[i]] <- ref(x, mfu, K)
  }
  expect_all_identical(got, want, "rebalance_cutpoints (random)")
  # degenerate inputs: too few events, all events equal (the 1e-6 bump), events all at the cap, NaN removed by na.rm
  sp <- list(list(c(1, 2), 120, 8L), list(rep(5, 60), 120, 8L), list(rep(120, 60), 120, 8L), list(rep(0, 60), 120, 8L),
             list(c(rep(3, 30), rep(4, 30)), 120, 12L), list(c(runif(60, 0, 10), NA), 120, 6L), list(runif(7), 120, 8L),
             list(runif(8), 120, 8L), list(c(0.5, 1e-7, 1e-7 + 1e-9, 2, 3, 4, 5), 24, 3L))
  expect_all_identical(lapply(sp, function(a) rs("rsem_rebalance_cutpoints")(a[[1]], a[[2]], a[[3]])),
                       lapply(sp, function(a) ref(a[[1]], a[[2]], a[[3]])), "rebalance_cutpoints (degenerate)")
})
