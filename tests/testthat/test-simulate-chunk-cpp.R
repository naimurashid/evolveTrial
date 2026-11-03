# Tests for the C++ chunk engine (engine C++ port, step 3): the sample() port, the supported-set predicate and the
# dispatcher, and chunk identity with the R body on the R01 capture under each return-value setting.
# The chunk tests need the step 2a capture dev/perf/evidence/args/args_R01.rds and skip without it.

ns <- function(name) getFromNamespace(name, "evolveTrial")
rng_now <- function() get(".Random.seed", envir = globalenv())

r01_path <- function() testthat::test_path("..", "..", "dev", "perf", "evidence", "args", "args_R01.rds")

# The run_simulation_pure frame for one scenario of a capture, built as run_scenarios' run_one builds args_i
# (simulation_driver.R, run_scenarios) and returned by the chunk_probe seam before any chunk runs.
build_frame <- function(d, scenario, R, return_variance = TRUE, return_percentiles = FALSE) {
  rsp <- ns("run_simulation_pure")
  args_i <- utils::modifyList(d$base_args, d$scens[[scenario]], keep.null = TRUE)
  args_i$return_percentiles <- return_percentiles
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

# A modified copy of a frame: the three chunk closures are re-pointed at the copy and cfg is repacked from it.
copy_frame <- function(env, modify = function(e) NULL) {
  e2 <- list2env(as.list(env, all.names = TRUE), parent = parent.env(env))
  modify(e2)
  for (nm in c("simulate_chunk", "simulate_chunk_cpp", "simulate_chunk_r")) {
    f <- get(nm, envir = e2, inherits = FALSE)
    environment(f) <- e2
    assign(nm, f, envir = e2)
  }
  e2$cfg <- ns(".pack_chunk_cfg")(e2)
  e2
}

run_both <- function(env, idx, seed, fn_names = c("simulate_chunk_r", "simulate_chunk_cpp")) {
  out <- lapply(fn_names, function(nm) {
    res <- get(nm, envir = env)(idx, seed = seed)
    list(res = res, rng = rng_now())
  })
  names(out) <- fn_names
  out
}

# ---- sample() port (plan P2) ----------------------------------------------------------------------------------------
test_that("the sample() port matches base::sample on 1, 2 and 3 arms, equal and unequal probabilities", {
  smp <- ns(".chunk_sample_one")
  cases <- list(1, 0.37, c(0.5, 0.5), c(0.3, 0.7), c(0.7, 0.3), c(1, 3), c(0, 1), rep(1 / 3, 3), c(0.2, 0.5, 0.3),
                c(0.5, 0.25, 0.25), c(2, 1, 1), c(0.5, 0, 0.5), c(0.1, 0.1, 0.8))
  n_bad <- 0L; n_all <- 0L
  for (p in cases) {
    arms <- c("A", "B", "C")[seq_along(p)]
    for (seed in 1:400) {
      set.seed(seed)
      want <- match(sample(arms, size = 1, prob = p), arms)
      rng_want <- rng_now()
      set.seed(seed)
      got <- smp(p)
      rng_got <- rng_now()
      n_all <- n_all + 1L
      if (!identical(got, want) || !identical(rng_got, rng_want)) n_bad <- n_bad + 1L
    }
  }
  expect_identical(n_bad, 0L)
  expect_identical(n_all, length(cases) * 400L)
  # a single eligible arm still consumes exactly one uniform
  set.seed(7); u <- runif(2)
  set.seed(7); smp(1); expect_identical(runif(1), u[2])
})

# ---- chunk identity on R01 ------------------------------------------------------------------------------------------
test_that("C++ chunk is identical to the R body on R01 with return_variance TRUE and FALSE and with percentiles", {
  skip_if_not(file.exists(r01_path()), "step 2a capture args_R01.rds not available")
  d <- readRDS(r01_path())
  settings <- list(c(rv = TRUE, rp = FALSE), c(rv = FALSE, rp = FALSE), c(rv = FALSE, rp = TRUE), c(rv = TRUE, rp = TRUE))
  for (s in seq_along(d$scens)) for (st in settings) {
    env <- build_frame(d, s, 40L, return_variance = st[["rv"]], return_percentiles = st[["rp"]])
    expect_true(isTRUE(env$cfg$supported))
    out <- run_both(env, 3:30, 424242L + s)
    expect_identical(out$simulate_chunk_cpp$res, out$simulate_chunk_r$res)
    expect_identical(out$simulate_chunk_cpp$rng, out$simulate_chunk_r$rng)
    expect_identical(is.null(out$simulate_chunk_cpp$res$final_events_raw), !st[["rv"]])
  }
})

test_that("C++ chunk is identical to the R body on the sequential path (seed = NULL continues the caller's stream)", {
  skip_if_not(file.exists(r01_path()), "step 2a capture args_R01.rds not available")
  d <- readRDS(r01_path())
  env <- build_frame(d, 2L, 30L)
  set.seed(99); a <- env$simulate_chunk_r(seq_len(30), seed = NULL); ra <- rng_now()
  set.seed(99); b <- env$simulate_chunk_cpp(seq_len(30), seed = NULL); rb <- rng_now()
  expect_identical(b, a)
  expect_identical(rb, ra)
})

# ---- the supported-set predicate and the dispatcher -----------------------------------------------------------------
test_that("each violated predicate item sends the dispatcher to the R body and makes simulate_chunk_cpp stop", {
  skip_if_not(file.exists(r01_path()), "step 2a capture args_R01.rds not available")
  d <- readRDS(r01_path())
  env <- build_frame(d, 1L, 20L)
  set_both <- function(e, name, value) {   # the frame variable and its copy in args / base_args_for_interim
    if (is.null(value)) {
      assign(name, NULL, envir = e)
      if (name %in% names(e$args)) { e$args[name] <- list(NULL); e$base_args_for_interim[name] <- list(NULL) }
    } else {
      assign(name, value, envir = e)
      if (name %in% names(e$args)) { e$args[[name]] <- value; e$base_args_for_interim[[name]] <- value }
    }
  }
  set_args <- function(e, name, value) {
    if (is.null(value)) { e$base_args_for_interim[name] <- list(NULL) } else e$base_args_for_interim[[name]] <- value
  }
  mods <- list(
    "1"   = function(e) set_both(e, "compare_arms_option", TRUE),
    "2"   = function(e) set_both(e, "arm_names", c("Experimental", "Other")),
    "3"   = function(e) set_both(e, "max_total_patients_per_arm", unname(e$max_total_patients_per_arm)),
    "4"   = function(e) assign("looks_stop_at_final", FALSE, envir = e),
    "5"   = function(e) assign("pt_targets_abs", 100, envir = e),
    "6"   = function(e) assign("random_dropout", TRUE, envir = e),
    "7"   = function(e) set_both(e, "max_trial_time", Inf),
    "8"   = function(e) set_both(e, "rebalance_after_events", NULL),
    "9"   = function(e) set_both(e, "efficacy_stopping_rule_hc", TRUE),
    "10a" = function(e) set_both(e, "median_pfs_success_threshold_arms", NULL),
    "10b" = function(e) assign("median_pfs_futility_threshold_arms", c(Other = 4), envir = e),
    "10c" = function(e) { set_args(e, "futility_median_arms", NULL); set_args(e, "median_pfs_success_threshold_arms", NULL);
                          set_args(e, "null_median_arms", 3.8) },
    "12"  = function(e) assign("show_progress", TRUE, envir = e),
    "13a" = function(e) set_both(e, "min_person_time_frac_per_arm", 0.1),
    "13b" = function(e) assign("diagnostics", TRUE, envir = e),
    "13c" = function(e) set_both(e, "futility_threshold_hc_prob", NULL),
    "13d" = function(e) assign("final_success_posterior_prob_threshold", NA_real_, envir = e),
    "13e" = function(e) set_both(e, "num_posterior_draws", 0),
    "13f" = function(e) set_both(e, "min_events_hc", c(3L, 4L)),
    "13g" = function(e) set_both(e, "prior_alpha_params_model", e$prior_alpha_params_model[-1])
  )
  dispatch_hits <- function(e) {
    hits <- new.env(); hits$cpp <- 0L; hits$r <- 0L
    e$simulate_chunk_cpp <- function(...) { hits$cpp <- hits$cpp + 1L; NULL }
    e$simulate_chunk_r <- function(...) { hits$r <- hits$r + 1L; NULL }
    e$simulate_chunk(1:3, seed = 5L)
    c(cpp = hits$cpp, r = hits$r)
  }
  # the unmodified copy dispatches to C++ (default option) and to R with the option FALSE
  e0 <- copy_frame(env)
  expect_identical(dispatch_hits(e0), c(cpp = 1L, r = 0L))
  old <- options(evolveTrial.use_cpp_engine = FALSE)
  expect_identical(dispatch_hits(copy_frame(env)), c(cpp = 0L, r = 1L))
  options(old)
  for (item in names(mods)) {
    e <- copy_frame(env, mods[[item]])
    real_cpp <- e$simulate_chunk_cpp
    expect_false(isTRUE(e$cfg$supported), info = item)
    expect_false(ns(".chunk_cpp_supported")(e), info = item)
    set.seed(11); r0 <- rng_now()
    expect_error(real_cpp(1:3, seed = 5L), "outside the supported set", info = item)
    expect_identical(rng_now(), r0, info = item)   # refused before set.seed: stream untouched
    expect_identical(dispatch_hits(e), c(cpp = 0L, r = 1L), info = item)
  }
  # item 11: EVOLVETRIAL_USE_CPP is read at dispatch time
  old_env <- Sys.getenv("EVOLVETRIAL_USE_CPP", unset = NA)
  Sys.setenv(EVOLVETRIAL_USE_CPP = "FALSE")
  on.exit(if (is.na(old_env)) Sys.unsetenv("EVOLVETRIAL_USE_CPP") else Sys.setenv(EVOLVETRIAL_USE_CPP = old_env), add = TRUE)
  e <- copy_frame(env)
  real_cpp <- e$simulate_chunk_cpp
  expect_true(isTRUE(e$cfg$supported))
  expect_false(ns(".chunk_cpp_supported")(e))
  expect_error(real_cpp(1:3, seed = 5L), "EVOLVETRIAL_USE_CPP")
  expect_identical(dispatch_hits(e), c(cpp = 0L, r = 1L))
})

# Step 4: between-arm frames are now checked against the vs-ref items. This frame leaves compare_arms_hr_margin at its
# NULL default, so it fails item V6 (it failed item 1 in step 3) and must still fall back.
test_that("an unsupported between-arm frame (no HR margin) packs without error and falls back", {
  rsp <- ns("run_simulation_pure")
  old <- options(evolveTrial.chunk_probe = function(env) env, evolveTrial.looks_stop_at_final = TRUE)
  on.exit(options(old), add = TRUE)
  env <- rsp(num_simulations = 10L, arm_names = c("Doublet", "Triplet"), reference_arm_name = "Doublet",
             compare_arms_option = TRUE, weibull_shape_true_arms = c(Doublet = 1.3, Triplet = 1.3),
             weibull_median_true_arms = c(Doublet = 6, Triplet = 9), interval_cutpoints_sim = seq(0, 24, by = 6),
             max_follow_up_sim = 120, censor_max_time_sim = 120, prior_alpha_params_model = rep(0.3, 4),
             prior_beta_params_model = rep(0.3, 4), num_posterior_draws = 50, cohort_size_per_arm = 1L,
             max_total_patients_per_arm = c(Doublet = 20, Triplet = 20), overall_accrual_rate = 2,
             randomization_probs = c(Doublet = 0.5, Triplet = 0.5), max_trial_time = 72, use_ph_model_vs_ref = TRUE,
             efficacy_threshold_vs_ref_prob = 0.95, futility_threshold_vs_ref_prob = 0.4, rebalance_after_events = 50)
  expect_false(isTRUE(env$cfg$supported))
  expect_true(any(grepl("^V6 ", env$cfg$violations)))
  expect_error(env$simulate_chunk_cpp(1:2, seed = 3L), "outside the supported set")
  a <- env$simulate_chunk(1:4, seed = 3L); ra <- rng_now()
  b <- env$simulate_chunk_r(1:4, seed = 3L); rb <- rng_now()
  expect_identical(a, b)
  expect_identical(ra, rb)
})
