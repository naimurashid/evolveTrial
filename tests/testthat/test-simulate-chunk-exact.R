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

# P3-1 (Codex review P3): the replay instantiation runs the exact file's own copy of the trial loop with R's stream and the
# Monte Carlo posteriors, so comparing it with the R body itself checks the exact engine's amd8 chronology (j = 0 look,
# cutoff, last-slot rule, phantom vs-ref arrival, post-enrolment looks), dropout and Inf trial time directly. Frames from
# tests/testthat/fixtures (always available), looks_stop_at_final passed explicitly.
build_frame_p31 <- function(D, scenario, R, lsaf, mod = list()) {
  d <- readRDS(testthat::test_path("fixtures", paste0("args_", D, ".rds")))
  rsp <- ns("run_simulation_pure")
  args_i <- utils::modifyList(utils::modifyList(d$base_args, mod, keep.null = TRUE), d$scens[[scenario]], keep.null = TRUE)
  args_i$return_variance <- TRUE
  args_i$looks_stop_at_final <- lsaf
  args_i <- args_i[intersect(names(args_i), names(formals(rsp)))]
  args_i$num_simulations <- as.integer(R)
  old <- options(evolveTrial.chunk_probe = function(env) env)
  on.exit(options(old), add = TRUE)
  do.call(rsp, args_i)
}

test_that("P3-1: the replay instantiation equals the R body on looks_stop_at_final FALSE frames (exact file's loop)", {
  half23 <- 0.5 * (52 + 51) / 2
  frames <- list(
    list(D = "R01", lsaf = FALSE, mod = list(), R = 60L),
    list(D = "R01", lsaf = FALSE, mod = list(max_trial_time = 10), R = 60L),
    list(D = "R01", lsaf = FALSE, mod = list(min_follow_up_at_final = 0), R = 60L),
    list(D = "R01", lsaf = FALSE, mod = list(max_trial_time = Inf, random_dropout = TRUE, censor_max_time_sim = 30,
                                             rebalance_after_events = NULL), R = 60L),
    list(D = "R01", lsaf = TRUE, mod = list(max_trial_time = Inf, futility_threshold_hc_prob = NULL,
                                            median_pfs_success_threshold_arms = NULL), R = 40L),
    list(D = "R23", lsaf = FALSE, mod = list(), R = 15L),
    list(D = "R23", lsaf = FALSE, mod = list(max_trial_time = half23), R = 15L),
    list(D = "R23", lsaf = FALSE, mod = list(min_follow_up_at_final = 0, random_dropout = TRUE), R = 15L)
  )
  for (f in frames) for (s in 1:2) {
    env <- build_frame_p31(f$D, s, f$R, f$lsaf, f$mod)
    info <- paste(f$D, s, f$lsaf, paste(names(f$mod), unlist(f$mod), collapse = " "))
    expect_true(isTRUE(env$cfg$supported), info = info)
    a <- env$simulate_chunk_r(seq_len(f$R), seed = 5150L); ra <- .Random.seed
    set.seed(5150L); b <- ns(".simulate_chunk_exact_replay")(seq_len(f$R), env$cfg); rb <- .Random.seed
    expect_identical(b, a, info = info)
    expect_identical(rb, ra, info = info)
  }
})

test_that("P3-1 (Codex P2): an arrival time that overflows to Inf errors in the R body, cpp and the exact engine", {
  d <- readRDS(testthat::test_path("fixtures", "args_R01.rds"))
  ba <- utils::modifyList(d$base_args, list(max_trial_time = Inf, overall_accrual_rate = 1e-308,
                                            min_follow_up_at_final = 0, max_total_patients_per_arm = c(Experimental = 1L),
                                            rebalance_after_events = NULL, futility_threshold_hc_prob = NULL,
                                            random_dropout = FALSE), keep.null = TRUE)
  # one arrival per trial; it overflows with probability Pr(Exp(1) > 1.798) = 0.17, so 60 trials overflow in every engine
  # (the exact engine has its own streams and errors on its own overflowing trial)
  ba$num_simulations <- 60L; ba$parallel_replicates <- FALSE
  run <- function(engine, use_cpp) {
    old <- options(evolveTrial.engine = engine, evolveTrial.use_cpp_engine = use_cpp, evolveTrial.looks_stop_at_final = FALSE)
    on.exit(options(old), add = TRUE)
    run_scenarios(ba, d$scens[1], seed = 1, return_variance = TRUE)
  }
  expect_error(run("cpp", FALSE), "missing value where TRUE/FALSE needed")   # the R body
  expect_error(run("cpp", TRUE), "missing value where TRUE/FALSE needed")
  expect_error(run("exact", TRUE), "missing value where TRUE/FALSE needed")
  expect_error(run(NULL, TRUE), "missing value where TRUE/FALSE needed")     # the default (exact)
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
  bad <- env$cfg; bad$random_dropout <- NA   # P3-1: TRUE is supported; a non-logical value is still refused
  expect_error(ns(".simulate_chunk_exact")(1:5, bad, c(1L, 2L)), "must be TRUE or FALSE")
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

# ---- review M1: the ExactPost wiring against the phase 1 decision rule -------------------------------------------------------
# .exact_trace runs the engine's trials with every exact probability recorded as (inputs, threshold, value). For each call
# site of src/simulate_chunk_exact.cpp the test checks (a) the threshold is the one phase 1 uses at that site
# (src/simulate_chunk.cpp:375 interim HC futility_thr_interim, below; :762 final HC success_thr, above, then futility_thr,
# below; :514 interim PH log_margin_interim, ge; :741-742 final PH log_eff_final, lt, then log_margin_final, ge), (b) the
# value equals the exact functions called directly on the recorded inputs, and (c) it agrees with the phase 1 Monte Carlo
# posterior at 10^5 draws within 4 standard errors. A swap of thresholds or of tails fails (a) or (b); the mutations tried
# are listed in dev/perf/phase2/REVIEW_FIXES.md. Fixed seeds, so the outcome is deterministic.
mc_close <- function(mc, p, N) {
  se <- max(sqrt(p * (1 - p) / N), 1 / N)
  abs(mc - p) / se
}

test_that("M1 HC: thresholds, tails and values at every exact call site (interim, final success, final futility)", {
  nrec <- c(`1` = 0L, `2` = 0L, `3` = 0L)
  for (spec in list(list("R03", 2), list("R01", 2), list("R04", 2))) {
    D <- spec[[1]]
    skip_if_not(file.exists(cap_path(D)), paste("capture", D, "not available"))
    d <- readRDS(cap_path(D))
    env <- build_frame_x(d, spec[[2]], 150)
    cfg <- env$cfg
    # the thresholds that identify the call sites must differ, or a swap could not be seen
    expect_false(any(cfg$success_thr %in% c(cfg$futility_thr, cfg$futility_thr_interim)))
    tr <- ns(".exact_trace")(1:150, cfg, c(31L, 7L))
    expect_identical(tr$result, ns(".simulate_chunk_exact")(1:150, cfg, c(31L, 7L)))   # tracing changes nothing
    P <- tr$post
    expect_true(all(P$kind %in% 1:3))
    expect_true(all(P$thr[P$kind == 1] %in% cfg$futility_thr_interim))      # interim: futility_thr_interim
    expect_true(all(P$thr[P$kind == 2] %in% cfg$success_thr))               # final efficacy: success_thr
    expect_true(all(P$thr[P$kind == 3] %in% cfg$futility_thr))              # final futility: futility_thr
    expect_true(all(P$n[P$kind == 1] == cfg$draws_interim))
    expect_true(all(P$n[P$kind %in% 2:3] == cfg$draws_final))
    # a final futility evaluation follows a final efficacy evaluation on the same inputs (phase 1 order)
    i3 <- which(P$kind == 3)
    expect_true(all(P$kind[i3 - 1L] == 2L))
    expect_identical(P$ev1[i3], P$ev1[i3 - 1L])
    expect_identical(P$pt1[i3], P$pt1[i3 - 1L])
    below <- ns(".exact_pr_median_below"); above <- ns(".exact_pr_median_above")
    direct <- vapply(seq_along(P$kind), function(i) {
      f <- if (P$kind[i] == 2L) above else below
      f(P$ev1[[i]], P$pt1[[i]], cfg$prior_alpha, cfg$prior_beta, P$lens[[i]], P$thr[i])
    }, numeric(1))
    expect_identical(P$prob, direct)
    for (kd in 1:3) nrec[as.character(kd)] <- nrec[as.character(kd)] + sum(P$kind == kd)
    # Monte Carlo: up to 6 records per kind, preferring probabilities away from 0 and 1
    set.seed(9000 + match(D, c("R01", "R03", "R04")))
    for (kd in 1:3) {
      ix <- which(P$kind == kd)
      if (!length(ix)) next
      ix <- ix[order(abs(P$prob[ix] - 0.5))][seq_len(min(6L, length(ix)))]
      for (i in ix) {
        lam <- evolveTrial::draw_posterior_hazard_samples_cpp(cfg$K, P$ev1[[i]], P$pt1[[i]], cfg$prior_alpha,
                                                              cfg$prior_beta, 100000L)
        med <- evolveTrial::calculate_median_survival_matrix_cpp(lam, P$lens[[i]])
        mc <- if (kd == 2L) mean(med > P$thr[i]) else mean(med < P$thr[i])
        expect_lt(mc_close(mc, P$prob[i], 1e5), 4)
      }
    }
  }
  expect_true(all(nrec > 0))   # all three HC call sites were exercised
})

test_that("M1 PH: thresholds, tails and values at every exact call site (interim, final efficacy, final futility)", {
  n_int <- 0L; n_fin <- 0L
  for (spec in list(list("R23", 2), list("R26", 2), list("R25", 1))) {
    D <- spec[[1]]
    skip_if_not(file.exists(cap_path(D)), paste("capture", D, "not available"))
    d <- readRDS(cap_path(D))
    env <- build_frame_x(d, spec[[2]], 150)
    cfg <- env$cfg
    # the interim/final margins may coincide in a design, but the efficacy threshold must differ from the margins
    expect_false(cfg$log_eff_final %in% c(cfg$log_margin_interim, cfg$log_margin_final))
    tr <- ns(".exact_trace")(1:150, cfg, c(31L, 7L))
    expect_identical(tr$result, ns(".simulate_chunk_exact")(1:150, cfg, c(31L, 7L)))
    P <- tr$post
    expect_true(all(P$kind %in% 4:5))
    k <- P$kind
    is_fin_eff <- which(k == 5L)                       # ph_lt is called only for the final efficacy probability
    expect_true(all(P$thr[is_fin_eff] == cfg$log_eff_final))
    # final order: ph_lt(log_eff_final) then ph_ge(log_margin_final) on the same inputs
    fin_fut <- is_fin_eff + 1L
    expect_true(all(k[fin_fut] == 4L))
    expect_true(all(P$thr[fin_fut] == cfg$log_margin_final))
    expect_true(all(P$bmean[fin_fut] == P$bmean[is_fin_eff]) && all(P$bsd[fin_fut] == P$bsd[is_fin_eff]))
    interim <- setdiff(which(k == 4L), fin_fut)
    expect_true(all(P$thr[interim] == cfg$log_margin_interim))      # interim futility: log_margin_interim, ge
    expect_true(all(P$n[interim] == cfg$draws_vsref_interim))
    expect_true(all(P$n[is_fin_eff] == cfg$draws_final))
    n_int <- n_int + length(interim); n_fin <- n_fin + length(is_fin_eff)
    # value: pnorm on the Laplace output computed here from the recorded inputs
    direct <- vapply(seq_along(k), function(i) {
      bp <- evolveTrial::ph_beta_mode_var_cpp(P$ev1[[i]], P$pt1[[i]], P$ev2[[i]], P$pt2[[i]], cfg$prior_alpha,
                                               cfg$prior_beta, cfg$ph_mu, cfg$ph_sigma, cfg$ph_tol, cfg$ph_max_iter)
      c(bp$mean, bp$sd, if (k[i] == 4L) pnorm(P$thr[i], bp$mean, bp$sd, lower.tail = FALSE) else pnorm(P$thr[i], bp$mean, bp$sd))
    }, numeric(3))
    expect_identical(rbind(P$bmean, P$bsd, P$prob), direct)
    # Monte Carlo: the phase 1 sampler, up to 6 records per kind
    set.seed(9100 + match(D, c("R23", "R25", "R26")))
    for (kd in 4:5) {
      ix <- which(k == kd)
      ix <- ix[order(abs(P$prob[ix] - 0.5))][seq_len(min(6L, length(ix)))]
      for (i in ix) {
        o <- ns(".chunk_ph_sample")(P$ev1[[i]], P$pt1[[i]], P$ev2[[i]], P$pt2[[i]], cfg$prior_alpha, cfg$prior_beta,
                                    P$lens[[i]], cfg$ph_mu, cfg$ph_sigma, cfg$ph_tol, cfg$ph_max_iter, 100000L)
        mc <- if (kd == 4L) mean(o$logHR >= P$thr[i]) else mean(o$logHR < P$thr[i])
        expect_lt(mc_close(mc, P$prob[i], 1e5), 4)
      }
    }
  }
  expect_gt(n_int, 0L); expect_gt(n_fin, 0L)
})

# ---- review M2: the patient-level samplers of CounterRng ---------------------------------------------------------------------
# The trace returns every draw the trial loop made: arrival gap, arm assignment (position among the eligible arms and their
# normalized probabilities) and event time with the shape and scale passed. The gap must be Exp(rate = accrual_rate), the
# event time Weibull(shape, scale) as passed (and the pair must be one of the design's arms), and the arm frequencies must
# follow the normalized probabilities. Kolmogorov-Smirnov and chi-square tests with fixed seeds, so deterministic; the pass
# threshold is p > 1e-3 (a correct sampler fails an unlucky seed with probability 1e-3 per test, and these seeds pass; a
# wrong parameterization gives p below 1e-100 at these sample sizes, see REVIEW_FIXES.md).
test_that("M2: arrival gaps, event times and arm assignments follow pexp, pweibull and the randomization probabilities", {
  skip_if_not(file.exists(cap_path("R23")), "capture R23 not available")
  d <- readRDS(cap_path("R23"))
  env <- build_frame_x(d, 2, 400)
  cfg <- env$cfg
  cfg$rand_probs <- c(0.25, 0.75)                      # unequal, so an arm mix-up or a flat draw is visible
  expect_gt(abs(cfg$shape[1] - cfg$shape[2]) + abs(cfg$scale[1] - cfg$scale[2]), 0)
  expect_false(isTRUE(all.equal(cfg$accrual_rate, 1)))   # rate used as scale would be invisible at rate 1
  tr <- ns(".exact_trace")(1:400, cfg, c(2026L, 1104L))
  D <- tr$draws
  expect_gt(length(D$gap), 5000L)
  # arrival gaps
  expect_true(all(D$gap_scale == 1 / cfg$accrual_rate))
  ks_gap <- stats::ks.test(stats::pexp(D$gap, rate = cfg$accrual_rate), "punif")
  expect_gt(ks_gap$p.value, 1e-3)
  # event times
  pairs <- paste(cfg$shape, cfg$scale)
  expect_true(all(paste(D$ev_shape, D$ev_scale) %in% pairs))
  u_ev <- stats::pweibull(D$ev_t, shape = D$ev_shape, scale = D$ev_scale)
  ks_ev <- stats::ks.test(u_ev, "punif")
  expect_gt(ks_ev$p.value, 1e-3)
  for (a in 1:2) {                                      # and per arm
    sel <- D$ev_shape == cfg$shape[a] & D$ev_scale == cfg$scale[a]
    if (sum(sel) > 200 && !identical(pairs[1], pairs[2])) expect_gt(stats::ks.test(u_ev[sel], "punif")$p.value, 1e-3)
  }
  # arm assignment: records with both arms eligible, chosen position against the normalized probabilities
  n2 <- D$arm_n == 2L
  expect_gt(sum(n2), 3000L)
  off <- c(0L, cumsum(D$arm_n))[seq_along(D$arm_n)]
  p1 <- D$arm_probs[off[n2] + 1L]; p2 <- D$arm_probs[off[n2] + 2L]
  expect_true(all(abs(p1 - 0.25) < 1e-15) && all(abs(p2 - 0.75) < 1e-15))
  cnt <- tabulate(D$arm_pick[n2], 2L)
  ch <- stats::chisq.test(cnt, p = c(0.25, 0.75))
  expect_gt(ch$p.value, 1e-3)
  # one arm draw and one event time per enrolled patient; the last gap of a trial may end the accrual loop
  expect_identical(length(D$arm_pick), length(D$ev_t))
  expect_gte(length(D$gap), length(D$ev_t))
})

# ---- review L1: quadrature non-convergence stops the engine ----------------------------------------------------------------
test_that("L1: an exact probability whose quadrature did not converge stops the run with an informative error", {
  skip_if_not(file.exists(cap_path("R01")), "capture R01 not available")
  d <- readRDS(cap_path("R01"))
  env <- build_frame_x(d, 2, 30)
  ok <- ns(".simulate_chunk_exact")(1:30, env$cfg, c(5L, 6L))
  force <- ns(".exact_test_force_nonconverged")
  expect_false(force(TRUE))
  on.exit(force(FALSE), add = TRUE)
  expect_error(ns(".simulate_chunk_exact")(1:30, env$cfg, c(5L, 6L)), "quadrature did not converge")
  expect_error(ns(".exact_trace")(1:30, env$cfg, c(5L, 6L)), "quadrature did not converge")
  expect_true(force(FALSE))
  expect_identical(ns(".simulate_chunk_exact")(1:30, env$cfg, c(5L, 6L)), ok)   # the seam is off again
})

# ---- review L2: one notice when the exact engine is requested on a path that does not use it ------------------------------
test_that("L2: engine 'exact' on a path without the exact engine gives one message and draws nothing", {
  skip_if_not(file.exists(cap_path("R01")), "capture R01 not available")
  d <- readRDS(cap_path("R01"))
  note <- ns(".exact_note_unused"); en <- ns(".exact_note_env")
  old <- options(evolveTrial.looks_stop_at_final = TRUE, evolveTrial.engine = "cpp", evolveTrial.use_cpp_engine = TRUE)
  on.exit({ options(old); en$noted <- NULL }, add = TRUE)
  en$noted <- NULL
  run <- function() {
    ba <- d$base_args; ba$num_simulations <- 20L; ba$parallel_replicates <- FALSE
    run_scenarios(ba, d$scens, seed = 3, return_variance = TRUE)
  }
  # explicit "cpp" engine (no longer the default since 0.0.0.9002): silent, and the notice stays armed
  expect_no_message(run())
  expect_null(en$noted)
  # exact engine in use: silent
  options(evolveTrial.engine = "exact")
  expect_no_message(run())
  expect_null(en$noted)
  # exact requested, R body selected: one message, then silence; results as the "cpp" engine with the R body
  options(evolveTrial.use_cpp_engine = FALSE)
  expect_message(a <- run(), "evolveTrial.engine.*exact.*use_cpp_engine is FALSE")
  expect_no_message(b <- run())
  options(evolveTrial.engine = "cpp")
  expect_identical(run(), a)
  expect_identical(a, b)
  # the notice draws nothing from R's generator
  en$noted <- NULL
  options(evolveTrial.engine = "exact")
  set.seed(77); r0 <- .Random.seed
  expect_message(note("a test path"), "a test path")
  expect_identical(.Random.seed, r0)
  expect_no_message(note("a test path"))
  # other simulators: the same notice, once
  en$noted <- NULL
  expect_message(try(run_simulation_binary(num_simulations = 5L, arm_names = "A", true_response_prob = 0.3, n1_per_arm = 10,
                                           n_total_per_arm = 20, p0 = 0.2), silent = TRUE), "run_simulation_binary")
  expect_true(en$noted)
  # an invalid option value is not reported here (run_simulation_pure rejects it) and does not error
  en$noted <- NULL
  options(evolveTrial.engine = "bogus")
  expect_no_message(note("x"))
})
