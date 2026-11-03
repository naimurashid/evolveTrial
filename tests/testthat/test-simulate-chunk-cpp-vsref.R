# Tests for the C++ chunk engine, between-arm (vs-reference, PH model) path (engine C++ port, step 4): the PH sampler
# port against the live R sampler, the restructured supported-set predicate, the dispatcher, and chunk identity with the
# R body on the R23 capture. The chunk tests need the step 2a captures dev/perf/evidence/args/args_R23.rds (and
# args_R01.rds for the predicate label check) and skip without them.

ns <- function(name) getFromNamespace(name, "evolveTrial")
rng_now <- function() get(".Random.seed", envir = globalenv())
set_rng <- function(s) assign(".Random.seed", s, envir = globalenv())

cap_path <- function(D) testthat::test_path("..", "..", "dev", "perf", "evidence", "args", paste0("args_", D, ".rds"))

# As in test-simulate-chunk-cpp.R: the run_simulation_pure frame for one scenario of a capture (run_scenarios' args_i).
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

run_both <- function(env, idx, seed) {
  out <- lapply(c("simulate_chunk_r", "simulate_chunk_cpp"), function(nm) {
    res <- get(nm, envir = env)(idx, seed = seed)
    list(res = res, rng = rng_now())
  })
  names(out) <- c("simulate_chunk_r", "simulate_chunk_cpp")
  out
}

# ---- the PH sampler port (plan P5, P15, P16) ------------------------------------------------------------------------
test_that("the PH sampler port matches sample_vs_ref_medians_ph on slices from R23 trials (interim and final calls)", {
  skip_if_not(file.exists(cap_path("R23")), "step 2a capture args_R23.rds not available")
  d <- readRDS(cap_path("R23"))
  calls <- list()
  rec <- function(slCtrl, slTrt, args, num_samples) {
    calls[[length(calls) + 1L]] <<- list(slCtrl = slCtrl, slTrt = slTrt, args = args, num_samples = num_samples,
                                         seed = rng_now())
  }
  nsenv <- asNamespace("evolveTrial")
  orig <- get("sample_vs_ref_medians_ph", envir = nsenv)
  wrapped <- function(slCtrl, slTrt, args, num_samples) {
    rec(slCtrl, slTrt, args, num_samples)
    orig(slCtrl, slTrt, args, num_samples)
  }
  unlockBinding("sample_vs_ref_medians_ph", nsenv)
  assign("sample_vs_ref_medians_ph", wrapped, envir = nsenv)
  on.exit({ assign("sample_vs_ref_medians_ph", orig, envir = nsenv); lockBinding("sample_vs_ref_medians_ph", nsenv) },
          add = TRUE)
  for (s in seq_along(d$scens)) {
    env <- build_frame(d, s, 20L)
    env$simulate_chunk_r(1:6, seed = 777L + s)
  }
  assign("sample_vs_ref_medians_ph", orig, envir = nsenv)
  expect_gt(length(calls), 20L)
  # the captured calls include rebalanced cutpoints (last cutpoint max_follow_up_sim) and the original ones
  last_cut <- vapply(calls, function(x) tail(x$args$interval_cutpoints_sim, 1), numeric(1))
  expect_true(any(last_cut == 120) && any(last_cut == 24))
  n_bad <- 0L
  for (x in calls) {
    set_rng(x$seed)
    want <- orig(x$slCtrl, x$slTrt, x$args, x$num_samples)
    rng_want <- rng_now()
    set_rng(x$seed)
    got <- ns(".chunk_ph_sample_r")(x$slCtrl, x$slTrt, x$args, x$num_samples)
    rng_got <- rng_now()
    if (!identical(got, want) || !identical(rng_got, rng_want)) n_bad <- n_bad + 1L
  }
  expect_identical(n_bad, 0L)
  # an empty control registry slice (double zero vectors, SM:118-130) and an empty treatment slice
  x <- calls[[1]]
  empty <- list(metrics = list(events_per_interval = rep(0, length(x$args$interval_cutpoints_sim) - 1),
                               person_time_per_interval = rep(0, length(x$args$interval_cutpoints_sim) - 1)))
  for (pair in list(list(empty, x$slTrt), list(x$slCtrl, empty), list(empty, empty))) {
    set.seed(31); want <- orig(pair[[1]], pair[[2]], x$args, 800); rw <- rng_now()
    set.seed(31); got <- ns(".chunk_ph_sample_r")(pair[[1]], pair[[2]], x$args, 800); rg <- rng_now()
    expect_identical(got, want)
    expect_identical(rg, rw)
  }
})

# ---- chunk identity on R23 ------------------------------------------------------------------------------------------
test_that("C++ chunk is identical to the R body on R23 with return_variance TRUE and FALSE and with percentiles", {
  skip_if_not(file.exists(cap_path("R23")), "step 2a capture args_R23.rds not available")
  d <- readRDS(cap_path("R23"))
  settings <- list(c(rv = TRUE, rp = FALSE), c(rv = FALSE, rp = FALSE), c(rv = FALSE, rp = TRUE), c(rv = TRUE, rp = TRUE))
  for (s in seq_along(d$scens)) for (st in settings) {
    env <- build_frame(d, s, 20L, return_variance = st[["rv"]], return_percentiles = st[["rp"]])
    expect_true(isTRUE(env$cfg$supported))
    expect_identical(env$cfg$mode, "vsref")
    out <- run_both(env, 3:12, 515151L + s)
    expect_identical(out$simulate_chunk_cpp$res, out$simulate_chunk_r$res)
    expect_identical(out$simulate_chunk_cpp$rng, out$simulate_chunk_r$rng)
    expect_identical(is.null(out$simulate_chunk_cpp$res$final_events_raw), !st[["rv"]])
  }
})

test_that("C++ chunk is identical to the R body on R23 on the sequential path (seed = NULL)", {
  skip_if_not(file.exists(cap_path("R23")), "step 2a capture args_R23.rds not available")
  d <- readRDS(cap_path("R23"))
  env <- build_frame(d, 2L, 10L)
  set.seed(99); a <- env$simulate_chunk_r(seq_len(10), seed = NULL); ra <- rng_now()
  set.seed(99); b <- env$simulate_chunk_cpp(seq_len(10), seed = NULL); rb <- rng_now()
  expect_identical(b, a)
  expect_identical(rb, ra)
})

test_that("C++ chunk is identical to the R body with the reference arm listed second", {
  skip_if_not(file.exists(cap_path("R23")), "step 2a capture args_R23.rds not available")
  d <- readRDS(cap_path("R23"))
  d$base_args$arm_names <- rev(d$base_args$arm_names)
  d$base_args$max_total_patients_per_arm <- rev(d$base_args$max_total_patients_per_arm)
  env <- build_frame(d, 2L, 10L)
  expect_true(isTRUE(env$cfg$supported))
  expect_identical(env$cfg$ref_idx, 1L)
  out <- run_both(env, 1:10, 2468L)
  expect_identical(out$simulate_chunk_cpp$res, out$simulate_chunk_r$res)
  expect_identical(out$simulate_chunk_cpp$rng, out$simulate_chunk_r$rng)
})

# ---- the restructured predicate -------------------------------------------------------------------------------------
test_that("the HC frame is checked against exactly the step 3 items and the vs-ref frame against common plus V items", {
  skip_if_not(file.exists(cap_path("R01")) && file.exists(cap_path("R23")), "step 2a captures not available")
  step3 <- c("1 compare_arms_option not TRUE", "2 one arm", "3 max_total_patients_per_arm named by arm_names",
             "4 looks_stop_at_final TRUE", "5 pt_targets_abs NULL", "6 random_dropout FALSE", "7 max_trial_time finite",
             "8 rebalance_after_events non-NULL", "9 efficacy_stopping_rule_hc not TRUE",
             "10 success and futility median thresholds named and finite", "12 show_progress FALSE",
             "13a min_person_time_frac_per_arm is the frozen unnamed 0 (2.49)", "13b diagnostics FALSE (3.5)",
             "13c interim futility_threshold_hc_prob one finite number (2.54)",
             "13d final posterior probability thresholds finite (2.25, 2.27)",
             "13e posterior draw counts at least 1 (2.56, 6.8)", "13f HC gate values are scalars (2.48)",
             "13g well-formed design inputs")
  env1 <- build_frame(readRDS(cap_path("R01")), 1L, 10L)
  expect_identical(names(ns(".chunk_cpp_static_checks")(env1)), step3)
  env23 <- build_frame(readRDS(cap_path("R23")), 1L, 10L)
  v_labels <- names(ns(".chunk_cpp_static_checks")(env23))
  expect_identical(v_labels[1:8], step3[c(4:8, 11, 13, 18)])
  expect_true(all(c(paste0("V", 1:10), "V11a") %in% sub(" .*", "", v_labels)))
  expect_identical(ns(".chunk_cpp_violations")(env23), character(0))
})

test_that("each violated vs-ref predicate item sends the dispatcher to the R body and makes simulate_chunk_cpp stop", {
  skip_if_not(file.exists(cap_path("R23")), "step 2a capture args_R23.rds not available")
  d <- readRDS(cap_path("R23"))
  env <- build_frame(d, 1L, 10L)
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
  # item label prefix -> modification
  mods <- list(
    "V1"   = function(e) set_args(e, "compare_arms_option", FALSE),
    "V2"   = function(e) set_args(e, "use_ph_model_vs_ref", FALSE),
    "V3"   = function(e) set_both(e, "reference_arm_name", "Other"),
    "V3"   = function(e) set_both(e, "arm_names", c("Doublet", "Triplet", "Third")),
    "V3"   = function(e) set_args(e, "reference_arm_name", "Triplet"),
    "V4"   = function(e) set_both(e, "max_total_patients_per_arm", unname(e$max_total_patients_per_arm)),
    "V4"   = function(e) set_both(e, "max_total_patients_per_arm", rev(e$max_total_patients_per_arm)),
    "V5"   = function(e) set_both(e, "efficacy_stopping_rule_vs_ref", TRUE),
    "V6"   = function(e) set_args(e, "compare_arms_hr_margin", NULL),
    "V6"   = function(e) set_args(e, "compare_arms_hr_eff_margin", 0.05),
    "V7"   = function(e) set_both(e, "futility_threshold_vs_ref_prob", NULL),
    "V7"   = function(e) assign("efficacy_threshold_vs_ref_prob", NA_real_, envir = e),
    "V8"   = function(e) set_args(e, "compare_arms_target_arm", "Triplet"),
    "V9"   = function(e) set_args(e, "min_events_per_arm", c(Doublet = 11, Triplet = 11)),
    "V9"   = function(e) set_args(e, "min_median_followup_per_arm", c(3, 3)),
    "V10"  = function(e) set_args(e, "ph_loghr_prior_sd", c(1, 2)),
    "13e"  = function(e) set_both(e, "num_posterior_draws", 0),
    "V11a" = function(e) set_both(e, "num_posterior_draws", 800.5),
    "4"    = function(e) assign("looks_stop_at_final", FALSE, envir = e),
    "5"    = function(e) assign("pt_targets_abs", 100, envir = e),
    "6"    = function(e) assign("random_dropout", TRUE, envir = e),
    "7"    = function(e) set_both(e, "max_trial_time", Inf),
    "8"    = function(e) set_both(e, "rebalance_after_events", NULL),
    "12"   = function(e) assign("show_progress", TRUE, envir = e),
    "13b"  = function(e) assign("diagnostics", TRUE, envir = e),
    "13g"  = function(e) set_both(e, "prior_alpha_params_model", e$prior_alpha_params_model[-1]),
    "13g"  = function(e) set_both(e, "weibull_scale_true_arms", c(Doublet = 1))
  )
  dispatch_hits <- function(e) {
    hits <- new.env(); hits$cpp <- 0L; hits$r <- 0L
    e$simulate_chunk_cpp <- function(...) { hits$cpp <- hits$cpp + 1L; NULL }
    e$simulate_chunk_r <- function(...) { hits$r <- hits$r + 1L; NULL }
    e$simulate_chunk(1:3, seed = 5L)
    c(cpp = hits$cpp, r = hits$r)
  }
  e0 <- copy_frame(env)
  expect_true(isTRUE(e0$cfg$supported))
  expect_identical(dispatch_hits(e0), c(cpp = 1L, r = 0L))
  old <- options(evolveTrial.use_cpp_engine = FALSE)
  expect_identical(dispatch_hits(copy_frame(env)), c(cpp = 0L, r = 1L))
  options(old)
  for (i in seq_along(mods)) {
    item <- names(mods)[i]
    e <- copy_frame(env, mods[[i]])
    info <- paste(item, i)
    v <- e$cfg$violations
    expect_true(any(sub(" .*", "", v) == item), info = info)
    real_cpp <- e$simulate_chunk_cpp
    expect_false(isTRUE(e$cfg$supported), info = info)
    expect_false(ns(".chunk_cpp_supported")(e), info = info)
    set.seed(11); r0 <- rng_now()
    expect_error(real_cpp(1:3, seed = 5L), "outside the supported set", info = info)
    expect_identical(rng_now(), r0, info = info)
    expect_identical(dispatch_hits(e), c(cpp = 0L, r = 1L), info = info)
  }
  # item 11: EVOLVETRIAL_USE_CPP is read at dispatch time (it also switches ph_beta_mode_var to its R fallback, P14)
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
