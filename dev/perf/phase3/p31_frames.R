# P3-1: acceptance frames for gates B, S and the G3 check. Built from the R01, R03 and R23 captures in
# tests/testthat/fixtures (base_args + scens). Each frame: design, a list of base_args modifications (NULL values set
# explicitly with keep.null), a vector of base_args names to drop (so run_simulation_pure's default applies), and the
# chronology (lsaf = looks_stop_at_final passed explicitly; NA = leave it to the option default, which is FALSE).
# Every frame runs with return_variance = TRUE.
p31_root <- function() normalizePath("~/research/evolveTrial-merge")
p31_fixture <- function(D) readRDS(file.path(p31_root(), "tests/testthat/fixtures", paste0("args_", D, ".rds")))

# The required formals of run_simulation_pure (no default). "unset" frames keep only these (plus return_variance).
P31_REQUIRED <- c("arm_names", "reference_arm_name", "compare_arms_option", "weibull_shape_true_arms",
                  "weibull_median_true_arms", "interval_cutpoints_sim", "max_follow_up_sim", "censor_max_time_sim",
                  "prior_alpha_params_model", "prior_beta_params_model", "num_posterior_draws", "cohort_size_per_arm",
                  "max_total_patients_per_arm", "overall_accrual_rate", "randomization_probs")

p31_frame_specs <- function() {
  half_cut <- function(D) { b <- p31_fixture(D)$base_args; 0.5 * sum(b$max_total_patients_per_arm) / b$overall_accrual_rate }
  thr_null <- list(median_pfs_success_threshold_arms = NULL, median_pfs_futility_threshold_arms = NULL)
  list(
    # ---- single-arm, R01 (32 patients at 2 per month, F = 6, max_trial_time 72, beat 1.29)
    B01_R01_off          = list(D = "R01", lsaf = FALSE, mod = list()),
    B02_R01_off_cut      = list(D = "R01", lsaf = FALSE, mod = list(max_trial_time = 10)),
    B03_R01_off_F0       = list(D = "R01", lsaf = FALSE, mod = list(min_follow_up_at_final = 0)),
    B04_R01_off_Fneg     = list(D = "R01", lsaf = FALSE, mod = list(min_follow_up_at_final = -1), expect_error = TRUE),
    B05_R01_off_inf      = list(D = "R01", lsaf = FALSE, mod = list(max_trial_time = Inf)),
    B06_R01_on_inf       = list(D = "R01", lsaf = TRUE,  mod = list(max_trial_time = Inf)),
    B07_R01_off_drop     = list(D = "R01", lsaf = FALSE, mod = list(random_dropout = TRUE, censor_max_time_sim = 30)),
    B08_R01_on_drop      = list(D = "R01", lsaf = TRUE,  mod = list(random_dropout = TRUE, censor_max_time_sim = 30)),
    # B09 has no rebalancing; B09c is the same frame with a threshold of 8 events, which does rebalance, so the pair shows
    # that the NULL flag switches off a branch the frame would otherwise take (R01 never reaches the captured 50 events)
    B09_R01_off_norebal  = list(D = "R01", lsaf = FALSE, mod = list(rebalance_after_events = NULL)),
    B09c_R01_off_rebal8  = list(D = "R01", lsaf = FALSE, mod = list(rebalance_after_events = 8)),
    B10_R01_off_succNULL = list(D = "R01", lsaf = FALSE, mod = list(median_pfs_success_threshold_arms = NULL)),
    # futility_threshold_hc_prob raised from 0.346 to 0.9: at 0.346 the interim looks stop every trial that would reach the
    # final futility test, so the -Inf fallback was hit 0 times (first gate B run); 0.9 lets trials reach it
    B11_R01_off_futNULL  = list(D = "R01", lsaf = FALSE, mod = list(median_pfs_futility_threshold_arms = NULL,
                                                                    futility_threshold_hc_prob = 0.9)),
    B12_R01_off_chain    = list(D = "R01", lsaf = FALSE, mod = c(thr_null, list(futility_median_arms = NULL)),
                                scen_mod = list(null_median_arms = c(Experimental = 4.2))),
    B13_R01_off_futprobNULL = list(D = "R01", lsaf = FALSE, mod = list(futility_threshold_hc_prob = NULL)),
    B14_R01_on_futprobNULL  = list(D = "R01", lsaf = TRUE,  mod = list(futility_threshold_hc_prob = NULL)),
    B15_R01_unset        = list(D = "R01", lsaf = NA, unset = TRUE),
    B16_R01_unset_thr    = list(D = "R01", lsaf = NA, unset = TRUE,
                                keep = c("median_pfs_success_threshold_arms", "median_pfs_futility_threshold_arms",
                                         "futility_threshold_hc_prob", "null_median_arms", "futility_median_arms")),
    # ---- single-arm, R03 (30 patients, beat 2.03, min_events_hc 11)
    B17_R03_off          = list(D = "R03", lsaf = FALSE, mod = list()),
    B18_R03_off_F0_cut   = list(D = "R03", lsaf = FALSE, mod = list(min_follow_up_at_final = 0, max_trial_time = 12)),
    B19_R03_unset_thr    = list(D = "R03", lsaf = NA, unset = TRUE,
                                keep = c("median_pfs_success_threshold_arms", "median_pfs_futility_threshold_arms",
                                         "futility_threshold_hc_prob", "null_median_arms", "futility_median_arms",
                                         "min_events_hc", "interim_calendar_beat")),
    # ---- between-arm, R23 (52 + 51 patients at 2 per month, F = 6, beat 1.01)
    B20_R23_off          = list(D = "R23", lsaf = FALSE, mod = list()),
    B21_R23_off_cut      = list(D = "R23", lsaf = FALSE, mod = list(max_trial_time = half_cut("R23"))),
    B22_R23_off_F0       = list(D = "R23", lsaf = FALSE, mod = list(min_follow_up_at_final = 0)),
    B23_R23_off_inf      = list(D = "R23", lsaf = FALSE, mod = list(max_trial_time = Inf)),
    B24_R23_on_inf       = list(D = "R23", lsaf = TRUE,  mod = list(max_trial_time = Inf)),
    B25_R23_off_drop     = list(D = "R23", lsaf = FALSE, mod = list(random_dropout = TRUE, censor_max_time_sim = 30)),
    B26_R23_on_drop      = list(D = "R23", lsaf = TRUE,  mod = list(random_dropout = TRUE, censor_max_time_sim = 30)),
    B27_R23_off_norebal  = list(D = "R23", lsaf = FALSE, mod = list(rebalance_after_events = NULL)),
    B28_R23_off_Fneg     = list(D = "R23", lsaf = FALSE, mod = list(min_follow_up_at_final = -1), expect_error = TRUE),
    B29_R23_unset_vsref  = list(D = "R23", lsaf = NA, unset = TRUE,
                                keep = c("use_ph_model_vs_ref", "compare_arms_hr_margin", "efficacy_threshold_vs_ref_prob",
                                         "futility_threshold_vs_ref_prob", "min_events_per_arm",
                                         "min_median_followup_per_arm"))
  )
}

# base_args and scens of a frame (the scens are the capture's, with scen_mod applied to each).
p31_frame_inputs <- function(spec) {
  d <- p31_fixture(spec$D)
  ba <- d$base_args
  ba$looks_stop_at_final <- NULL
  if (isTRUE(spec$unset)) {
    ba <- ba[intersect(names(ba), c(P31_REQUIRED, spec$keep %||% character(0)))]
  } else {
    ba <- utils::modifyList(ba, spec$mod %||% list(), keep.null = TRUE)
  }
  if (!is.na(spec$lsaf)) ba$looks_stop_at_final <- spec$lsaf
  scens <- d$scens
  if (isTRUE(spec$unset)) {
    keep <- c(P31_REQUIRED, spec$keep %||% character(0))
    scens <- lapply(scens, function(s) s[intersect(names(s), keep)])
  }
  if (!is.null(spec$scen_mod)) scens <- lapply(scens, function(s) utils::modifyList(s, spec$scen_mod, keep.null = TRUE))
  list(base_args = ba, scens = scens)
}

# The run_simulation_pure frame (chunk_probe seam) for one scenario, R trials, return_variance TRUE.
p31_chunk_frame <- function(spec, scenario, R) {
  inp <- p31_frame_inputs(spec)
  rsp <- get("run_simulation_pure", envir = asNamespace("evolveTrial"))
  args_i <- utils::modifyList(inp$base_args, inp$scens[[scenario]], keep.null = TRUE)
  args_i$return_percentiles <- FALSE
  args_i$return_variance <- TRUE
  args_i$parallel_replicates <- FALSE
  args_i$progress <- FALSE
  args_i <- args_i[intersect(names(args_i), names(formals(rsp)))]
  args_i$num_simulations <- as.integer(R)
  old <- options(evolveTrial.chunk_probe = function(env) env)
  on.exit(options(old), add = TRUE)
  env <- do.call(rsp, args_i)
  stopifnot(is.environment(env), is.function(env$simulate_chunk_r))
  env
}
`%||%` <- function(x, y) if (is.null(x)) y else x
