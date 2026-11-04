# Regression tests for the two chunk engines against the original R functions (engine speed-up, default switch
# 2026-11-04; report dev/perf/phase2/DEFAULT_SWITCH_REPORT.md).
#
#   1. phase 1 ("cpp") against the R body, end to end through run_scenarios(): identical() tables and RNG state;
#   2. exact posterior tail probabilities against Monte Carlo from the ORIGINAL R samplers (draw_posterior_hazard_samples_r
#      with calculate_median_survival_piecewise, ph_beta_mode_var_r with rnorm, sample_vs_ref_medians_ph_r);
#   3. exact engine against the R body on operating characteristics, within Monte Carlo tolerances;
#   4. golden values of both engines at fixed seeds (tests/testthat/fixtures/engine_golden.rds);
#   5. default engine selection and the silent fallback.
#
# Designs are the step 2a captures R01 (single-arm, K = 8), R03 (single-arm, K = 9) and R23 (between-arm PH, K = 15),
# copied to tests/testthat/fixtures/ so the file runs without dev/. Each capture has a null scenario (1) and an
# alternative scenario (2). The larger versions of 3 run only when EVOLVETRIAL_LONG_TESTS is non-empty (the R body costs
# about 13 s per 1,000 trials on R01 and 75 s per 1,000 on R23).
#
# Already covered elsewhere and not repeated here (chunk level, through the chunk_probe seam):
#   test-simulate-chunk-cpp.R        C++ chunk identical to the R body on R01, return_variance TRUE/FALSE, percentiles,
#                                    sequential path, and the dispatcher's fallback on each unsupported item;
#   test-simulate-chunk-cpp-vsref.R  the same on R23, reference arm listed second, and vs-ref fallback items;
#   test-simulate-chunk-exact.R      exact engine split invariance, replay identity with phase 1, explicit-option refusal;
#   test-posterior-exact.R           exact posteriors against pgamma, integrate and the C++ samplers.

ns <- function(name) getFromNamespace(name, "evolveTrial")
fx <- function(D) readRDS(testthat::test_path("fixtures", paste0("args_", D, ".rds")))
golden_path <- function() testthat::test_path("fixtures", "engine_golden.rds")
long_tests <- function() nzchar(Sys.getenv("EVOLVETRIAL_LONG_TESTS"))

# Run all scenarios of a capture sequentially. engine = NULL leaves options(evolveTrial.engine) unset (the default).
# Returns the table as a plain data.frame and R's RNG state after the call.
run_design <- function(d, n, seed, engine = NULL, use_cpp = TRUE, rv = TRUE, looks_stop = TRUE, base_mod = list()) {
  old <- options(evolveTrial.engine = engine, evolveTrial.use_cpp_engine = use_cpp,
                 evolveTrial.looks_stop_at_final = looks_stop, evolveTrial.exact_block_size = NULL)
  on.exit(options(old), add = TRUE)
  ba <- utils::modifyList(d$base_args, base_mod, keep.null = TRUE)
  ba$num_simulations <- as.integer(n); ba$parallel_replicates <- FALSE
  tab <- run_scenarios(ba, d$scens, seed = seed, return_variance = rv)
  list(tab = as.data.frame(tab), rng = .Random.seed)
}

# ---- 1. phase 1 against the R body, end to end -----------------------------------------------------------------------
test_that("engine 'cpp' is identical to the R body through run_scenarios (R01, R23; null and alternative; rv T/F)", {
  for (spec in list(list(D = "R01", n = 40L), list(D = "R23", n = 10L))) {
    d <- fx(spec$D)
    for (rv in c(TRUE, FALSE)) {
      a <- run_design(d, spec$n, seed = 31, engine = "cpp", rv = rv)
      b <- run_design(d, spec$n, seed = 31, engine = "cpp", use_cpp = FALSE, rv = rv)
      info <- paste(spec$D, "rv", rv)
      expect_identical(a, b, info = info)
      expect_identical(sort(unique(a$tab$scenario)), 1:2, info = info)
      expect_identical(any(grepl("^Var_", names(a$tab))), rv, info = info)
    }
  }
})

# ---- 2. exact posteriors against the original R samplers -------------------------------------------------------------
# z = (MC - p) / sqrt(p (1 - p) / n_draws); |z| < 4. The floor 1 / n_draws in the variance keeps p near 0 or 1 from
# giving a zero denominator (then the MC tally must sit within about 4 draws of n p). Fixed seeds, so deterministic.
mc_z <- function(mc, p, n) (mc - p) / sqrt(max(p * (1 - p), 1 / n) / n)

test_that("single-arm exact Pr(median < thr) matches draw_posterior_hazard_samples_r + calculate_median_survival_piecewise", {
  pm <- ns(".exact_pr_median")
  states <- list(
    # single interval, weak prior
    list(E = 5L, PT = 40, a = 0.3, b = 0.3, lens = 24, thr = c(4, 6, 9)),
    # single interval, strong prior
    list(E = 2L, PT = 10, a = 8, b = 40, lens = 24, thr = c(3, 4.5)),
    # R01-like PWE (K = 8, 3-month intervals, prior 0.217), early interim with empty late intervals
    list(E = c(3L, 2L, 1L, 0L, 0L, 0L, 0L, 0L), PT = c(28, 15, 6, 1, 0, 0, 0, 0), a = rep(0.2167, 8), b = rep(0.2167, 8),
         lens = rep(3, 8), thr = c(3.8, 4.527, 7)),
    # R01-like PWE, later look with most intervals informed
    list(E = c(6L, 5L, 4L, 3L, 2L, 1L, 1L, 0L), PT = c(80, 62, 45, 32, 20, 12, 6, 2), a = rep(0.2167, 8),
         b = rep(0.2167, 8), lens = rep(3, 8), thr = c(3.8, 6, 9)),
    # PWE, no events at all (posterior median mostly beyond the grid: the last-hazard extension)
    list(E = rep(0L, 4), PT = c(12, 6, 2, 0), a = rep(0.5, 4), b = rep(0.5, 4), lens = rep(6, 4), thr = c(10, 24, 40)),
    # PWE, large prior and unequal interval lengths
    list(E = c(4L, 2L, 1L), PT = c(30, 18, 5), a = c(5, 5, 5), b = c(25, 25, 25), lens = c(2, 4, 18), thr = c(3, 4, 6))
  )
  n <- 40000L
  set.seed(20261104)
  for (i in seq_along(states)) {
    s <- states[[i]]
    K <- length(s$E)
    lam <- ns("draw_posterior_hazard_samples_r")(K, s$E, s$PT, s$a, s$b, num_samples = n)
    med <- apply(lam, 1L, ns("calculate_median_survival_piecewise"), interval_lengths = s$lens)
    for (thr in s$thr) {
      got <- pm(s$E, s$PT, s$a, s$b, s$lens, thr)
      expect_true(isTRUE(got$converged), info = paste(i, thr))
      expect_lt(abs(got$below + got$above - 1), 1e-12)
      expect_lt(abs(mc_z(mean(med < thr), got$below, n)), 4, label = sprintf("state %d thr %g z", i, thr))
    }
  }
})

test_that("between-arm exact Pr(logHR >= m) and Pr(logHR < m) match ph_beta_mode_var_r + rnorm", {
  plr <- ns(".exact_pr_loghr")
  K <- 15L; lens <- rep(1.6, K)
  shape <- function(tot, decay) as.integer(round(tot * decay^(0:(K - 1)) / sum(decay^(0:(K - 1)))))
  states <- list(
    # R23-like early look, few events, prior sd 1
    list(EC = shape(8, 0.7), ET = shape(6, 0.7), PTC = 40 * 0.8^(0:(K - 1)), PTT = 42 * 0.8^(0:(K - 1)),
         a = 0.305, mu = 0, sigma = 1),
    # later look, more events
    list(EC = shape(30, 0.85), ET = shape(22, 0.85), PTC = 120 * 0.85^(0:(K - 1)), PTT = 135 * 0.85^(0:(K - 1)),
         a = 0.305, mu = 0, sigma = 1),
    # tight log-HR prior
    list(EC = shape(12, 0.8), ET = shape(9, 0.8), PTC = 60 * 0.8^(0:(K - 1)), PTT = 66 * 0.8^(0:(K - 1)),
         a = 0.305, mu = 0, sigma = 0.2),
    # large hazard prior, prior mean away from zero
    list(EC = shape(10, 0.8), ET = shape(14, 0.8), PTC = 50 * 0.8^(0:(K - 1)), PTT = 48 * 0.8^(0:(K - 1)),
         a = 5, mu = log(0.8), sigma = 1)
  )
  m <- c(log(0.6), log(1 - 0.0338), 0, log(1.0338), log(1.5))
  n <- 100000L
  set.seed(4711)
  for (i in seq_along(states)) {
    s <- states[[i]]
    a <- rep(s$a, K)
    got <- plr(s$EC, s$PTC, s$ET, s$PTT, a, a, lens, s$mu, s$sigma, 1e-6, 50L, m)
    bp <- ns("ph_beta_mode_var_r")(s$EC, s$PTC, s$ET, s$PTT, a, a, s$mu, s$sigma)
    expect_equal(c(got$mean, got$sd), c(bp$mean, bp$sd), tolerance = 1e-6, info = paste("state", i))
    draws <- rnorm(n, bp$mean, bp$sd)
    for (j in seq_along(m)) {
      expect_lt(abs(got$pr_ge[j] + got$pr_lt[j] - 1), 1e-12)
      expect_lt(abs(mc_z(mean(draws >= m[j]), got$pr_ge[j], n)), 4, label = sprintf("state %d m %d z", i, j))
    }
  }
})

test_that("between-arm exact log-HR tails match the full original sampler sample_vs_ref_medians_ph_r", {
  plr <- ns(".exact_pr_loghr")
  d <- fx("R23")
  args <- d$base_args
  K <- length(args$prior_alpha_params_model)
  sl <- function(E, PT) list(metrics = list(events_per_interval = E, person_time_per_interval = PT))
  EC <- c(4L, 3L, 3L, 2L, 1L, 1L, 1L, rep(0L, K - 7L)); ET <- c(3L, 2L, 2L, 1L, 1L, 1L, rep(0L, K - 6L))
  PTC <- 45 * 0.82^(0:(K - 1)); PTT <- 47 * 0.82^(0:(K - 1))
  m <- c(log(0.7), log(1 - args$compare_arms_hr_margin), 0, log(1.4))
  got <- plr(EC, PTC, ET, PTT, args$prior_alpha_params_model, args$prior_beta_params_model,
             diff(args$interval_cutpoints_sim), args$ph_loghr_prior_mean, args$ph_loghr_prior_sd, 1e-6, 50L, m)
  n <- 20000L
  set.seed(99)
  o <- ns("sample_vs_ref_medians_ph_r")(sl(EC, PTC), sl(ET, PTT), args, n)
  for (j in seq_along(m)) {
    expect_lt(abs(mc_z(mean(o$logHR >= m[j]), got$pr_ge[j], n)), 4, label = sprintf("m %d z", j))
  }
})

# ---- 3. exact engine against the R body on operating characteristics -------------------------------------------------
# For each arm row and each metric X in {Type_I_Error_or_Power, PET_Futility, Exp_N}, the two runs are independent
# (different random streams), so with V = the returned Var_X (the Monte Carlo variance of the estimate, s^2 / n):
#   z = (X_exact - X_R) / sqrt(V_exact + V_R),   required |z| < 4.
# Known direction: the R body uses finite posterior draws, which stop slightly more often for futility, so the one-sided
# checks are  power_exact >= power_R - 4 SE  and  PET_Futility_exact <= PET_Futility_R + 4 SE  (implied by |z| < 4, kept
# explicit so a sign flip in a future engine is named in the failure). When both variances are zero (a metric that is
# 0 or 1 in both runs, e.g. the reference arm's power) the values must be equal.
oc_compare <- function(D, n_exact, n_r, seed_exact, seed_r) {
  d <- fx(D)
  ex <- run_design(d, n_exact, seed = seed_exact, engine = "exact")$tab
  rb <- run_design(d, n_r, seed = seed_r, engine = "cpp", use_cpp = FALSE)$tab
  expect_identical(ex[, c("Arm_Name", "scenario")], rb[, c("Arm_Name", "scenario")])
  for (X in c("Type_I_Error_or_Power", "PET_Futility", "Exp_N")) {
    se <- sqrt(ex[[paste0("Var_", X)]] + rb[[paste0("Var_", X)]])
    diff <- ex[[X]] - rb[[X]]
    for (r in seq_len(nrow(ex))) {
      info <- sprintf("%s scenario %d %s %s: exact %.4f R %.4f se %.4f", D, ex$scenario[r], ex$Arm_Name[r], X,
                      ex[[X]][r], rb[[X]][r], se[r])
      if (se[r] == 0) {
        expect_identical(ex[[X]][r], rb[[X]][r], info = info)
      } else {
        expect_lt(abs(diff[r] / se[r]), 4, label = info)
      }
    }
  }
  se_pow <- sqrt(ex$Var_Type_I_Error_or_Power + rb$Var_Type_I_Error_or_Power)
  expect_true(all(ex$Type_I_Error_or_Power >= rb$Type_I_Error_or_Power - 4 * se_pow))
  se_fut <- sqrt(ex$Var_PET_Futility + rb$Var_PET_Futility)
  expect_true(all(ex$PET_Futility <= rb$PET_Futility + 4 * se_fut))
  invisible(list(exact = ex, r = rb))
}

test_that("exact engine agrees with the R body on power, type I error and expected N (routine size)", {
  oc_compare("R01", n_exact = 4000L, n_r = 200L, seed_exact = 11, seed_r = 12)
  oc_compare("R23", n_exact = 2000L, n_r = 60L, seed_exact = 21, seed_r = 22)
})

test_that("exact engine agrees with the R body on operating characteristics (long: EVOLVETRIAL_LONG_TESTS)", {
  skip_if_not(long_tests(), "set EVOLVETRIAL_LONG_TESTS to run")
  oc_compare("R01", n_exact = 40000L, n_r = 4000L, seed_exact = 111, seed_r = 112)
  oc_compare("R03", n_exact = 40000L, n_r = 4000L, seed_exact = 131, seed_r = 132)
  oc_compare("R23", n_exact = 20000L, n_r = 1000L, seed_exact = 121, seed_r = 122)
})

test_that("engine 'cpp' and the R body stay identical at larger size and on the parallel path (long)", {
  skip_if_not(long_tests(), "set EVOLVETRIAL_LONG_TESTS to run")
  skip_on_os("windows")
  for (spec in list(list(D = "R01", n = 600L), list(D = "R23", n = 200L))) {
    d <- fx(spec$D)
    expect_identical(run_design(d, spec$n, seed = 41, engine = "cpp"),
                     run_design(d, spec$n, seed = 41, engine = "cpp", use_cpp = FALSE), info = spec$D)
    par <- function(use_cpp) {
      old <- options(evolveTrial.engine = "cpp", evolveTrial.use_cpp_engine = use_cpp, evolveTrial.looks_stop_at_final = TRUE)
      on.exit(options(old), add = TRUE)
      ba <- d$base_args; ba$num_simulations <- spec$n; ba$parallel_replicates <- TRUE; ba$num_workers <- 2L
      as.data.frame(run_scenarios(ba, d$scens, seed = 42, return_variance = TRUE))
    }
    expect_identical(par(TRUE), par(FALSE), info = paste(spec$D, "parallel"))
  }
})

# ---- 4. golden values --------------------------------------------------------------------------------------------------
# Full run_scenarios() tables and the RNG state after the call, at fixed seeds, for both engines. Any change to either
# engine's output fails here. Tables are compared at tolerance 1e-10 (a changed trial outcome moves a proportion by at
# least 1 / n, far above it; the tolerance only absorbs last-bit libm differences across platforms); RNG states must be
# identical. To regenerate after an intended change: EVOLVETRIAL_UPDATE_GOLDEN=1, run this file, and commit the RDS.
golden_specs <- list(
  exact_R01 = list(D = "R01", engine = "exact", n = 500L, seed = 101, rv = TRUE),
  exact_R03 = list(D = "R03", engine = "exact", n = 500L, seed = 103, rv = FALSE),
  exact_R23 = list(D = "R23", engine = "exact", n = 300L, seed = 123, rv = TRUE),
  cpp_R01 = list(D = "R01", engine = "cpp", n = 200L, seed = 201, rv = TRUE),
  cpp_R03 = list(D = "R03", engine = "cpp", n = 200L, seed = 203, rv = FALSE),
  cpp_R23 = list(D = "R23", engine = "cpp", n = 100L, seed = 223, rv = TRUE),
  # P3-1 frames (looks_stop_at_final FALSE unless stated, the amd8 chronology); the six entries above are unchanged
  exact_R01_p31 = list(D = "R01", engine = "exact", n = 500L, seed = 311, rv = TRUE, lsaf = FALSE,
                       mod = list(max_trial_time = Inf, random_dropout = TRUE, censor_max_time_sim = 30,
                                  rebalance_after_events = NULL, min_follow_up_at_final = 0)),
  exact_R03_p31 = list(D = "R03", engine = "exact", n = 500L, seed = 313, rv = TRUE, lsaf = FALSE,
                       mod = list(median_pfs_success_threshold_arms = NULL, futility_threshold_hc_prob = NULL,
                                  max_trial_time = 12)),
  exact_R23_p31 = list(D = "R23", engine = "exact", n = 300L, seed = 323, rv = TRUE, lsaf = FALSE,
                       mod = list(random_dropout = TRUE, censor_max_time_sim = 30, rebalance_after_events = NULL,
                                  max_trial_time = 26)),
  cpp_R01_p31 = list(D = "R01", engine = "cpp", n = 200L, seed = 411, rv = TRUE, lsaf = FALSE,
                     mod = list(random_dropout = TRUE, censor_max_time_sim = 30, max_trial_time = 10)),
  cpp_R03_p31 = list(D = "R03", engine = "cpp", n = 200L, seed = 413, rv = TRUE, lsaf = TRUE,
                     mod = list(max_trial_time = Inf, median_pfs_futility_threshold_arms = NULL, futility_threshold_hc_prob = NULL)),
  cpp_R23_p31 = list(D = "R23", engine = "cpp", n = 100L, seed = 423, rv = TRUE, lsaf = FALSE,
                     mod = list(min_follow_up_at_final = 0, max_trial_time = Inf))
)
golden_run <- function(g) run_design(fx(g$D), g$n, seed = g$seed, engine = g$engine, rv = g$rv,
                                     looks_stop = if (is.null(g$lsaf)) TRUE else g$lsaf,
                                     base_mod = if (is.null(g$mod)) list() else g$mod)

if (nzchar(Sys.getenv("EVOLVETRIAL_UPDATE_GOLDEN"))) {
  saveRDS(lapply(golden_specs, golden_run), golden_path(), version = 2)
}

test_that("golden values: exact and cpp engines reproduce the stored outputs", {
  skip_if_not(file.exists(golden_path()), "golden fixture not available")
  gold <- readRDS(golden_path())
  expect_identical(names(gold), names(golden_specs))
  for (k in names(golden_specs)) {
    now <- golden_run(golden_specs[[k]])
    expect_equal(now$tab, gold[[k]]$tab, tolerance = 1e-10, info = k)
    expect_identical(now$rng, gold[[k]]$rng, info = k)
  }
  # the two engines are different engines: their stored outputs differ
  expect_false(isTRUE(all.equal(gold$exact_R01$tab, gold$cpp_R01$tab)))
})

# ---- 5. default selection ------------------------------------------------------------------------------------------------
test_that("default (option unset) runs the exact engine on supported paths, silently", {
  d <- fx("R01")
  en <- ns(".exact_note_env"); en$noted <- NULL
  on.exit(en$noted <- NULL, add = TRUE)
  expect_null(getOption("evolveTrial.engine"))
  expect_no_message(def <- run_design(d, 60L, seed = 7))
  expect_identical(def, run_design(d, 60L, seed = 7, engine = "exact"))
  expect_false(identical(def$tab, run_design(d, 60L, seed = 7, engine = "cpp")$tab))
  # one sample.int(.Machine$integer.max, 2) draw per scenario, as the exact engine documents
  set.seed(7); for (s in seq_along(d$scens)) sample.int(.Machine$integer.max, 2L, replace = TRUE)
  expect_identical(def$rng, .Random.seed)
  expect_null(en$noted)
  d23 <- fx("R23")
  expect_identical(run_design(d23, 30L, seed = 8), run_design(d23, 30L, seed = 8, engine = "exact"))
})

test_that("default falls back silently to the phase 1 dispatcher on unsupported paths, as engine 'cpp' does", {
  en <- ns(".exact_note_env"); en$noted <- NULL
  on.exit(en$noted <- NULL, add = TRUE)
  d <- fx("R01"); d23 <- fx("R23")
  # (a) interim efficacy stopping (item 9, single-arm) and person-time milestones (item 5, between-arm). Before P3-1 this
  #     case used looks_stop_at_final FALSE, which both engines now run.
  eff <- list(efficacy_stopping_rule_hc = TRUE, efficacy_threshold_hc_prob = 0.99)
  expect_no_message(a <- run_design(d, 15L, seed = 9, looks_stop = FALSE, base_mod = eff))
  expect_identical(a, run_design(d, 15L, seed = 9, engine = "cpp", looks_stop = FALSE, base_mod = eff))
  ptm <- list(person_time_milestones = c(0.3, 0.6))
  expect_no_message(a <- run_design(d23, 6L, seed = 9, looks_stop = FALSE, base_mod = ptm))
  expect_identical(a, run_design(d23, 6L, seed = 9, engine = "cpp", looks_stop = FALSE, base_mod = ptm))
  # (b) a between-arm frame without the HR margin (vs-ref item V6)
  nomargin <- list(compare_arms_hr_margin = NULL)
  expect_no_message(a <- run_design(d23, 6L, seed = 10, base_mod = nomargin))
  expect_identical(a, run_design(d23, 6L, seed = 10, engine = "cpp", base_mod = nomargin))
  # (c) use_cpp_engine FALSE: the R body, no notice under the default
  expect_no_message(a <- run_design(d, 15L, seed = 11, use_cpp = FALSE))
  expect_identical(a, run_design(d, 15L, seed = 11, engine = "cpp", use_cpp = FALSE))
  # (d) EVOLVETRIAL_USE_CPP FALSE (item 11): the R body
  old_env <- Sys.getenv("EVOLVETRIAL_USE_CPP", unset = NA)
  Sys.setenv(EVOLVETRIAL_USE_CPP = "FALSE")
  on.exit(if (is.na(old_env)) Sys.unsetenv("EVOLVETRIAL_USE_CPP") else Sys.setenv(EVOLVETRIAL_USE_CPP = old_env), add = TRUE)
  expect_no_message(a <- run_design(d, 15L, seed = 12))
  expect_identical(a, run_design(d, 15L, seed = 12, engine = "cpp"))
  if (is.na(old_env)) Sys.unsetenv("EVOLVETRIAL_USE_CPP") else Sys.setenv(EVOLVETRIAL_USE_CPP = old_env)
  expect_null(en$noted)
  # an explicit "exact" still refuses the unsupported frame (the distinction the default relies on)
  expect_error(run_design(d, 30L, seed = 9, engine = "exact", looks_stop = FALSE, base_mod = eff), "outside the supported set")
})

test_that("the L2 notice fires only for an explicit option, never for the default", {
  en <- ns(".exact_note_env"); note <- ns(".exact_note_unused")
  en$noted <- NULL
  on.exit({ en$noted <- NULL; options(evolveTrial.engine = NULL) }, add = TRUE)
  options(evolveTrial.engine = NULL)
  expect_no_message(note("hybrid path"))
  expect_no_message(try(run_simulation_binary(num_simulations = 5L, arm_names = "A", true_response_prob = 0.3,
                                              n1_per_arm = 10, n_total_per_arm = 20, p0 = 0.2), silent = TRUE))
  expect_null(en$noted)
  options(evolveTrial.engine = "exact")
  expect_message(note("hybrid path"), "hybrid path")
  expect_true(en$noted)
})

# ---- P3-0, revised in P3-1: amd8 chronology frames with a binding cutoff -----------------------------------------------
# P3-0 kept every looks_stop_at_final = FALSE frame on the R body. P3-1 ported the amd8 chronology to both engines, so the
# same frames now run phase 1 (identical to the R body) and, under the default, the exact engine.
test_that("P6-off frame with a binding max_trial_time: cpp equals the R body, the default runs the exact engine", {
  en <- ns(".exact_note_env"); en$noted <- NULL
  on.exit(en$noted <- NULL, add = TRUE)
  # R01: 32 patients at 2 per month (about 16 months of accrual); R23 at its own rate, cut well before accrual ends
  for (spec in list(list(D = "R01", n = 20L, cut = 10), list(D = "R23", n = 6L, cut = NULL))) {
    d <- fx(spec$D)
    cut <- if (is.null(spec$cut)) {
      0.5 * sum(d$base_args$max_total_patients_per_arm) / d$base_args$overall_accrual_rate
    } else spec$cut
    mod <- list(max_trial_time = cut)
    a <- run_design(d, spec$n, seed = 13, engine = "cpp", looks_stop = FALSE, base_mod = mod)
    b <- run_design(d, spec$n, seed = 13, use_cpp = FALSE, looks_stop = FALSE, base_mod = mod)
    expect_identical(a, b, info = spec$D)
    # the cutoff binds: expected enrolment per arm stays below the arm's maximum
    expect_true(all(a$tab$Exp_N < d$base_args$max_total_patients_per_arm[a$tab$Arm_Name]), info = spec$D)
    expect_no_message(def <- run_design(d, spec$n, seed = 13, looks_stop = FALSE, base_mod = mod))
    expect_identical(def, run_design(d, spec$n, seed = 13, engine = "exact", looks_stop = FALSE, base_mod = mod), info = spec$D)
    expect_true(all(def$tab$Exp_N < d$base_args$max_total_patients_per_arm[def$tab$Arm_Name]), info = spec$D)
  }
  expect_null(en$noted)
})

# ---- P3-1: the first slice of the general loop ----------------------------------------------------------------------------
# Items 4 (amd8 chronology), 6 (dropout on), 7 (max_trial_time Inf), 8 (no rebalancing), 10 (median-threshold fallbacks)
# and 13c (no interim futility rule) run on both engines. Gate B (dev/perf/phase3/P3-1_REPORT.md) is the full check with
# branch counts; these frames are its routine-size subset.
p31_frames <- function(D) {
  if (D == "R23") {
    half <- 0.5 * sum(fx("R23")$base_args$max_total_patients_per_arm) / fx("R23")$base_args$overall_accrual_rate
    return(list(
      off = list(lsaf = FALSE, mod = list()),
      off_cut = list(lsaf = FALSE, mod = list(max_trial_time = half)),
      off_F0 = list(lsaf = FALSE, mod = list(min_follow_up_at_final = 0)),
      off_inf = list(lsaf = FALSE, mod = list(max_trial_time = Inf)),
      on_inf = list(lsaf = TRUE, mod = list(max_trial_time = Inf)),
      off_drop = list(lsaf = FALSE, mod = list(random_dropout = TRUE, censor_max_time_sim = 30)),
      on_drop = list(lsaf = TRUE, mod = list(random_dropout = TRUE, censor_max_time_sim = 30)),
      off_norebal = list(lsaf = FALSE, mod = list(rebalance_after_events = NULL))))
  }
  list(
    off = list(lsaf = FALSE, mod = list()),
    off_cut = list(lsaf = FALSE, mod = list(max_trial_time = 10)),
    off_F0 = list(lsaf = FALSE, mod = list(min_follow_up_at_final = 0)),
    off_inf = list(lsaf = FALSE, mod = list(max_trial_time = Inf)),
    on_inf = list(lsaf = TRUE, mod = list(max_trial_time = Inf)),
    off_drop = list(lsaf = FALSE, mod = list(random_dropout = TRUE, censor_max_time_sim = 30)),
    on_drop = list(lsaf = TRUE, mod = list(random_dropout = TRUE, censor_max_time_sim = 30)),
    off_norebal = list(lsaf = FALSE, mod = list(rebalance_after_events = NULL)),
    off_succ_null = list(lsaf = FALSE, mod = list(median_pfs_success_threshold_arms = NULL)),
    off_fut_null = list(lsaf = FALSE, mod = list(median_pfs_futility_threshold_arms = NULL)),
    off_chain = list(lsaf = FALSE, mod = list(median_pfs_success_threshold_arms = NULL, median_pfs_futility_threshold_arms = NULL,
                                              futility_median_arms = NULL)),
    off_futprob_null = list(lsaf = FALSE, mod = list(futility_threshold_hc_prob = NULL)),
    on_futprob_null = list(lsaf = TRUE, mod = list(futility_threshold_hc_prob = NULL)))
}
# every option of run_simulation_pure left at its default except the required arguments (return_variance aside)
p31_unset <- function(d, keep = character(0)) {
  req <- c("arm_names", "reference_arm_name", "compare_arms_option", "weibull_shape_true_arms", "weibull_median_true_arms",
           "interval_cutpoints_sim", "max_follow_up_sim", "censor_max_time_sim", "prior_alpha_params_model",
           "prior_beta_params_model", "num_posterior_draws", "cohort_size_per_arm", "max_total_patients_per_arm",
           "overall_accrual_rate", "randomization_probs", keep)
  d$base_args <- d$base_args[intersect(names(d$base_args), req)]
  d$scens <- lapply(d$scens, function(s) s[intersect(names(s), req)])
  d
}
run_unset <- function(d, n, seed, engine = NULL, use_cpp = TRUE) {
  old <- options(evolveTrial.engine = engine, evolveTrial.use_cpp_engine = use_cpp, evolveTrial.looks_stop_at_final = NULL,
                 evolveTrial.exact_block_size = NULL)
  on.exit(options(old), add = TRUE)
  ba <- d$base_args; ba$num_simulations <- as.integer(n); ba$parallel_replicates <- FALSE
  tab <- run_scenarios(ba, d$scens, seed = seed, return_variance = TRUE)
  list(tab = as.data.frame(tab), rng = .Random.seed)
}
p31_check <- function(D, n, seed) {
  d <- fx(D)
  for (nm in names(p31_frames(D))) {
    f <- p31_frames(D)[[nm]]
    info <- paste(D, nm)
    a <- run_design(d, n, seed = seed, engine = "cpp", looks_stop = f$lsaf, base_mod = f$mod)
    b <- run_design(d, n, seed = seed, use_cpp = FALSE, looks_stop = f$lsaf, base_mod = f$mod)
    expect_identical(a, b, info = info)
    ex <- run_design(d, n, seed = seed, looks_stop = f$lsaf, base_mod = f$mod)   # default: the exact engine
    expect_identical(ex, run_design(d, n, seed = seed, engine = "exact", looks_stop = f$lsaf, base_mod = f$mod), info = info)
    expect_false(anyNA(a$tab), info = info)
    expect_false(anyNA(ex$tab), info = info)
  }
}

test_that("P3-1 frames: cpp is identical to the R body, the default runs the exact engine, no NA reaches the table", {
  p31_check("R01", 15L, 3101)
  p31_check("R23", 5L, 3123)
})

test_that("P3-1 frames, long tier: larger sizes and R03 (EVOLVETRIAL_LONG_TESTS)", {
  skip_if_not(long_tests(), "set EVOLVETRIAL_LONG_TESTS to run")
  p31_check("R01", 150L, 3201)
  p31_check("R03", 150L, 3203)
  p31_check("R23", 40L, 3223)
  for (D in c("R01", "R23")) {
    du <- p31_unset(fx(D), keep = if (D == "R23") c("use_ph_model_vs_ref", "compare_arms_hr_margin",
                                                     "efficacy_threshold_vs_ref_prob", "futility_threshold_vs_ref_prob"))
    expect_identical(run_unset(du, 60L, 3211, engine = "cpp"), run_unset(du, 60L, 3211, use_cpp = FALSE), info = D)
  }
})

test_that("P3-1: a call with every option unset runs the exact engine; cpp equals the R body", {
  en <- ns(".exact_note_env"); en$noted <- NULL
  on.exit(en$noted <- NULL, add = TRUE)
  d <- p31_unset(fx("R01"))
  rsp <- ns("run_simulation_pure")
  # the frame is admitted: looks_stop_at_final FALSE, dropout on, max_trial_time Inf, no rebalancing, NULL thresholds
  old <- options(evolveTrial.chunk_probe = function(env) env, evolveTrial.looks_stop_at_final = NULL)
  args_i <- utils::modifyList(d$base_args, d$scens[[1]]); args_i$num_simulations <- 5L; args_i$return_variance <- TRUE
  env <- do.call(rsp, args_i[intersect(names(args_i), c(names(formals(rsp)), "num_simulations"))])
  options(old)
  expect_identical(ns(".chunk_cpp_violations")(env), character(0))
  expect_true(isTRUE(env$cfg$supported))
  expect_false(env$cfg$looks_stop_at_final)
  expect_true(env$cfg$random_dropout)
  expect_identical(env$cfg$max_trial_time, Inf)
  expect_false(env$cfg$rebalance)
  expect_identical(env$cfg$rebalance_threshold, NA_real_)   # explicit flag; the threshold is never read
  expect_false(env$cfg$futility_interim)
  expect_identical(env$cfg$success_thr, Inf)
  expect_identical(env$cfg$futility_thr, -Inf)
  expect_identical(run_unset(d, 20L, 3301, engine = "cpp"), run_unset(d, 20L, 3301, use_cpp = FALSE))
  expect_no_message(def <- run_unset(d, 200L, 3302))
  expect_identical(def, run_unset(d, 200L, 3302, engine = "exact"))
  expect_false(anyNA(def$tab))
  # no thresholds: every trial ends inconclusive at the final analysis and none stops early
  expect_true(all(def$tab$Pr_Final_Inconclusive == 1))
  expect_null(en$noted)
})

test_that("P3-1: min_follow_up_at_final < 0 trips the amd8 invariant with the R body's message in every engine", {
  for (D in c("R01", "R23")) {
    d <- fx(D)
    mod <- list(min_follow_up_at_final = -1)
    msg <- "Chronology violation: enrollment at"
    for (eng in list(list(engine = "cpp", use_cpp = FALSE), list(engine = "cpp", use_cpp = TRUE),
                     list(engine = "exact", use_cpp = TRUE))) {
      expect_error(run_design(d, 4L, seed = 3401, engine = eng$engine, use_cpp = eng$use_cpp, looks_stop = FALSE,
                              base_mod = mod), msg, info = paste(D, eng$engine, eng$use_cpp))
    }
    # the same message text, number for number, from the R body and from phase 1
    m_r <- tryCatch(run_design(d, 4L, seed = 3401, use_cpp = FALSE, looks_stop = FALSE, base_mod = mod),
                    error = conditionMessage)
    m_c <- tryCatch(run_design(d, 4L, seed = 3401, engine = "cpp", looks_stop = FALSE, base_mod = mod),
                    error = conditionMessage)
    expect_identical(m_c, m_r, info = D)
  }
})

test_that("P3-1: one frame of each still-excluded item falls back to the R body (2, 5, 9, V2, V5, V6, V9)", {
  en <- ns(".exact_note_env"); en$noted <- NULL
  on.exit(en$noted <- NULL, add = TRUE)
  d <- fx("R01"); d23 <- fx("R23")
  hc <- list(
    "5" = list(person_time_milestones = c(0.4)),
    "9" = list(efficacy_stopping_rule_hc = TRUE, efficacy_threshold_hc_prob = 0.99)
  )
  vs <- list(
    "V2" = list(use_ph_model_vs_ref = FALSE),
    "V5" = list(efficacy_stopping_rule_vs_ref = TRUE),
    "V6" = list(compare_arms_hr_margin = NULL),
    "V9" = list(min_events_per_arm = c(Doublet = 11, Triplet = 11))
  )
  chk <- function(dd, mods, n) for (k in names(mods)) {
    expect_no_message(a <- run_design(dd, n, seed = 3501, looks_stop = FALSE, base_mod = mods[[k]]))
    expect_identical(a, run_design(dd, n, seed = 3501, use_cpp = FALSE, looks_stop = FALSE, base_mod = mods[[k]]), info = k)
    expect_error(run_design(dd, n, seed = 3501, engine = "exact", looks_stop = FALSE, base_mod = mods[[k]]),
                 "outside the supported set", info = k)
  }
  chk(d, hc, 6L)
  chk(d23, vs, 3L)
  # item 2: two arms, each against a historical control (compare_arms_option FALSE)
  two <- list(compare_arms_option = FALSE)
  expect_identical(run_design(d23, 3L, seed = 3502, looks_stop = FALSE, base_mod = two),
                   run_design(d23, 3L, seed = 3502, use_cpp = FALSE, looks_stop = FALSE, base_mod = two))
  # items 3 and V4 (unnamed caps) make the R body itself error in make_state() end to end, so they are checked at the
  # frame level (test-simulate-chunk-cpp.R, test-simulate-chunk-cpp-vsref.R), as are V3 and item 12 (progress bar)
  expect_null(en$noted)
})

# ---- Codex review (P2): invalid priors and overflowed posterior scales ------------------------------------------------------
# R01 with zero enrolment caps: the final analysis sees the prior only. (a) Negative prior shapes are outside the supported
# set (item 13g): the default runs the R body, an explicit exact engine refuses. (b) A prior rate of 1e-310 makes w / rate
# overflow: the exact engine reaches the R body's decision (final futility) instead of a silent inconclusive. (c) The same
# with prior shape 0.01 is not resolved by the huge-scale rule: the exact engine stops with an error.
test_that("Codex P2: invalid priors fall back to the R body; overflowed scales decide as the R body; unresolved tails stop", {
  en <- ns(".exact_note_env"); en$noted <- NULL
  on.exit(en$noted <- NULL, add = TRUE)
  d <- fx("R01"); K <- length(d$base_args$prior_alpha_params_model)
  zero <- d$base_args$max_total_patients_per_arm; zero[] <- 0
  dec <- c("Pr_Final_Efficacy", "Pr_Final_Futility", "Pr_Final_Inconclusive")
  # (a)
  for (m in list(list(prior_alpha_params_model = rep(-1, K)), list(prior_alpha_params_model = rep(0, K)),
                 list(prior_beta_params_model = c(0, rep(1, K - 1))))) {
    mod <- c(list(max_total_patients_per_arm = zero), m)
    # a zero prior rate makes the R body warn (posterior rate replaced by 1e-6) in every engine that runs it
    expect_no_message(a <- suppressWarnings(run_design(d, 10L, seed = 3, base_mod = mod)))
    expect_identical(a, suppressWarnings(run_design(d, 10L, seed = 3, use_cpp = FALSE, base_mod = mod)))
    expect_identical(a, suppressWarnings(run_design(d, 10L, seed = 3, engine = "cpp", base_mod = mod)))
    expect_error(run_design(d, 10L, seed = 3, engine = "exact", base_mod = mod), "13g well-formed design inputs")
  }
  neg <- run_design(d, 10L, seed = 3, use_cpp = FALSE,
                    base_mod = list(max_total_patients_per_arm = zero, prior_alpha_params_model = rep(-1, K)))
  expect_true(all(neg$tab$Pr_Final_Efficacy == 1))   # the R body's decision in the Codex reproduction
  # (b)
  mod <- list(max_total_patients_per_arm = zero, prior_alpha_params_model = rep(1, K),
              prior_beta_params_model = c(1e-310, rep(1, K - 1)))
  ex <- run_design(d, 10L, seed = 3, base_mod = mod)
  rb <- run_design(d, 10L, seed = 3, use_cpp = FALSE, base_mod = mod)
  expect_identical(ex$tab[, dec], rb$tab[, dec])
  expect_true(all(ex$tab$Pr_Final_Futility == 1))
  # (c)
  mod$prior_alpha_params_model <- rep(0.01, K)
  expect_error(run_design(d, 10L, seed = 3, base_mod = mod), "exact posterior probability is not finite")
})
