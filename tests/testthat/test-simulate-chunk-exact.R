# Tests for the phase 2 chunk engine (engine speed-up phase 2, step P2-2; src/simulate_chunk_exact.cpp,
# R/simulate_chunk_exact.R): the counter-based generator, the replay identity that ties the copied trial logic to the
# phase 1 engine, split invariance, the reducer, and the driver option. The chunk tests use the step 2a captures
# dev/perf/evidence/args/args_R01.rds and args_R23.rds and skip without them.

ns <- function(name) getFromNamespace(name, "evolveTrial")
cap_path <- function(D) testthat::test_path("..", "..", "dev", "perf", "evidence", "args", paste0("args_", D, ".rds"))

build_frame_x <- function(d, scenario, R, return_variance = TRUE) {
  rsp <- ns("run_simulation_pure")
  args_i <- utils::modifyList(d$base_args, d$scens[[scenario]], keep.null = TRUE)
  args_i$return_percentiles <- FALSE
  args_i$percentile_probs <- c(0, 0.25, 0.5, 0.75, 0.9, 1.0)
  args_i$return_variance <- return_variance
  args_i <- args_i[intersect(names(args_i), names(formals(rsp)))]
  args_i$num_simulations <- as.integer(R)
  old <- options(evolveTrial.chunk_probe = function(env) env, evolveTrial.looks_stop_at_final = TRUE)
  on.exit(options(old), add = TRUE)
  env <- do.call(rsp, args_i)
  stopifnot(is.environment(env), isTRUE(env$looks_stop_at_final))
  env
}

# ---- generator -------------------------------------------------------------------------------------------------------
test_that("xoshiro256++ matches the reference outputs for state (1, 2, 3, 4)", {
  # independent implementation (Python integers) of the published algorithm; the first output is rotl(5, 23) + 1
  out <- ns(".exact_rng_raw")(c(0, 1, 0, 2, 0, 3, 0, 4), 5L)
  expect_identical(out, c(0, 41943041, 0, 58720359, 835584, 58720359, 836097, 2571370674, 2148704769, 2596549581))
})

test_that("per-trial uniforms match the documented seeding and lie in (0, 1)", {
  u <- ns(".exact_rng_uniforms")(c(5L, 7L), 3L, 3L)
  # reference from an independent implementation of the seeding in the header of src/simulate_chunk_exact.cpp
  expect_identical(u, c(0.7864174013990898, 0.612948370221046, 0.7474315999287281))
  big <- ns(".exact_rng_uniforms")(c(123L, 456L), 1L, 100000L)
  expect_true(all(big > 0 & big < 1))
  expect_gt(stats::ks.test(big, "punif")$p.value, 1e-3)   # deterministic input, so not flaky
  expect_false(identical(u, ns(".exact_rng_uniforms")(c(5L, 7L), 4L, 3L)))
  expect_false(identical(u, ns(".exact_rng_uniforms")(c(5L, 8L), 3L, 3L)))
  expect_identical(u, ns(".exact_rng_uniforms")(c(5L, 7L), 3L, 3L))
  # first uniforms of consecutive trials are not correlated
  first <- vapply(1:4000, function(i) ns(".exact_rng_uniforms")(c(9L, 9L), i, 1L), numeric(1))
  expect_lt(abs(stats::cor(first[-1], first[-4000])), 0.06)
  expect_gt(stats::ks.test(first, "punif")$p.value, 1e-3)
})

# ---- the trial logic is phase 1's ---------------------------------------------------------------------------------------
test_that("the replay instantiation reproduces the phase 1 engine exactly (R01, R23, both scenarios)", {
  for (D in c("R01", "R23")) {
    skip_if_not(file.exists(cap_path(D)), paste("capture", D, "not available"))
    d <- readRDS(cap_path(D))
    for (s in seq_along(d$scens)) {
      env <- build_frame_x(d, s, 120)
      set.seed(4242); a <- ns(".simulate_chunk_cpp")(1:120, env$cfg); ra <- .Random.seed
      set.seed(4242); b <- ns(".simulate_chunk_exact_replay")(1:120, env$cfg); rb <- .Random.seed
      expect_identical(b, a)
      expect_identical(rb, ra)
    }
  }
})

# ---- exact engine --------------------------------------------------------------------------------------------------------
test_that("the exact engine depends only on the seed and the trial index, and draws nothing from R", {
  for (D in c("R01", "R23")) {
    skip_if_not(file.exists(cap_path(D)), paste("capture", D, "not available"))
    d <- readRDS(cap_path(D))
    env <- build_frame_x(d, 2, 150)
    ex <- ns(".simulate_chunk_exact"); red <- ns(".exact_reduce")
    set.seed(1); r0 <- .Random.seed
    whole <- ex(1:150, env$cfg, c(11L, 22L))
    expect_identical(.Random.seed, r0)
    parts <- list(ex(1:7, env$cfg, c(11L, 22L)), ex(8:100, env$cfg, c(11L, 22L)), ex(101:150, env$cfg, c(11L, 22L)))
    a <- red(list(whole), env$cfg, 150L)
    expect_identical(red(parts, env$cfg, 150L), a)
    expect_identical(red(rev(parts), env$cfg, 150L), a)          # order of blocks does not matter
    one <- ex(c(77L, 3L), env$cfg, c(11L, 22L))                   # any trial on its own
    expect_identical(one$final_n[1, ], whole$final_n[77, ])
    expect_identical(one$final_time[2], whole$final_time[3])
    expect_false(identical(red(list(ex(1:150, env$cfg, c(11L, 23L))), env$cfg, 150L), a))
    expect_named(a, names(ns(".simulate_chunk_cpp")(1:5, env$cfg)))
    expect_identical(a$sum_final_n, colSums(whole$final_n))         # integer-valued, so exact
    expect_identical(a$final_n_raw[[1]], unname(whole$final_n[, 1]))
    expect_error(red(parts[-2], env$cfg, 150L), "missing")
    expect_error(red(c(parts, parts[1]), env$cfg, 150L), "twice")
  }
})

test_that("the exact engine refuses configurations outside the supported set", {
  skip_if_not(file.exists(cap_path("R01")), "capture R01 not available")
  d <- readRDS(cap_path("R01"))
  env <- build_frame_x(d, 1, 20)
  bad <- env$cfg; bad$random_dropout <- TRUE
  expect_error(ns(".simulate_chunk_exact")(1:5, bad, c(1L, 2L)), "outside the supported set")
  bad <- env$cfg; bad$supported <- FALSE
  expect_error(ns(".simulate_chunk_exact")(1:5, bad, c(1L, 2L)), "outside the supported set")
  expect_error(ns(".simulate_chunk_exact")(1:5, env$cfg, 1L), "two non-NA integers")
})

# ---- driver -----------------------------------------------------------------------------------------------------------
test_that("run_scenarios with engine 'exact': one sample.int draw per scenario, identical over workers and blocks", {
  skip_if_not(file.exists(cap_path("R23")), "capture R23 not available")
  skip_on_os("windows")
  d <- readRDS(cap_path("R23"))
  old <- options(evolveTrial.engine = "exact", evolveTrial.looks_stop_at_final = TRUE, evolveTrial.exact_block_size = NULL)
  on.exit(options(old), add = TRUE)
  go <- function(w, par, bs = NULL, seed = 5) {
    options(evolveTrial.exact_block_size = bs)
    ba <- d$base_args; ba$num_simulations <- 200L; ba$num_workers <- w; ba$parallel_replicates <- par
    tab <- run_scenarios(ba, d$scens, seed = seed, return_variance = TRUE)
    list(tab = tab, rng = .Random.seed)
  }
  a <- go(1, FALSE)
  set.seed(5); for (s in seq_along(d$scens)) sample.int(.Machine$integer.max, 2L, replace = TRUE)
  expect_identical(a$rng, .Random.seed)
  expect_identical(go(2, TRUE, 9), a)
  expect_identical(go(3, TRUE, 64), a)
  expect_false(identical(go(1, FALSE, seed = 6)$tab, a$tab))
  options(evolveTrial.engine = "bogus")
  expect_error(go(1, FALSE), "must be \"cpp\" or \"exact\"")
})

test_that("use_cpp_engine = FALSE selects the R body whatever the engine option", {
  skip_if_not(file.exists(cap_path("R01")), "capture R01 not available")
  d <- readRDS(cap_path("R01"))
  old <- options(evolveTrial.use_cpp_engine = FALSE, evolveTrial.looks_stop_at_final = TRUE)
  on.exit(options(old), add = TRUE)
  go <- function(eng) {
    options(evolveTrial.engine = eng)
    ba <- d$base_args; ba$num_simulations <- 40L; ba$parallel_replicates <- FALSE
    list(run_scenarios(ba, d$scens, seed = 3, return_variance = TRUE), .Random.seed)
  }
  expect_identical(go("exact"), go("cpp"))
})

# ---- step P2-4: the engine's slicer against the 3bb0694 re-slicer ------------------------------------------------------
test_that("the P2-4 slicer reproduces the reference re-slicer exactly on random registries and looks", {
  cmp <- ns(".exact_slice_compare")
  set.seed(20261104)
  n_case <- 0L
  n_bad <- 0L
  for (it in seq_len(3000)) {
    n <- sample(c(0:3, sample.int(300, 1)), 1)
    K <- sample(1:8, 1)
    max_fu <- sample(c(24, 36, 60, Inf), 1)
    enroll <- sort(cumsum(stats::rexp(n, 1)))
    te <- stats::rweibull(n, shape = runif(1, 0.5, 2), scale = runif(1, 3, 40))
    rc <- if (runif(1) < 0.7) rep(Inf, n) else stats::rexp(n, 1 / 30)
    if (n > 3 && runif(1) < 0.2) te[sample.int(n, 2)] <- te[1]                        # tied times
    cal <- if (n > 0) runif(1, -1, max(enroll) + 80) else runif(1, 0, 50)
    if (n > 0 && runif(1) < 0.1) cal <- enroll[sample.int(n, 1)]                     # look at an enrollment time
    cutmax <- if (is.finite(max_fu)) max_fu else 60
    cut <- switch(sample(4, 1),
                  seq(0, cutmax, length.out = K + 1),                                  # base cutpoints
                  c(0, sort(runif(K - 1, 0, cutmax)), cutmax),                          # rebalanced-like
                  c(0, sort(sample(c(1, 2, 3), K - 1, replace = TRUE)), cutmax),        # repeated cutpoints
                  c(0, rev(sort(runif(K - 1, 0, cutmax))), cutmax))                     # not sorted (reference path)
    if (n > 0 && runif(1) < 0.1) te[1] <- cut[min(2, K + 1)]                           # event exactly at a cutpoint
    if (n > 0 && runif(1) < 0.02) te[sample.int(n, 1)] <- NaN                          # NaN time (NA median path)
    out <- cmp(enroll, te, rc, cal, max_fu, cut, K)
    if (!identical(out$new, out$ref)) n_bad <- n_bad + 1L
    n_case <- n_case + 1L
  }
  expect_identical(n_case, 3000L)
  expect_identical(n_bad, 0L)
})
