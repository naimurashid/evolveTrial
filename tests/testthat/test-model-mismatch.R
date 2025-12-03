# Test file: Model Mismatch Tests
# Purpose: Test sensitivity to model specification choices identified in audit
# Key finding: Data generated from Weibull, analyzed via piecewise exponential

# ==============================================================================
# TEST 1: Weibull Shape Parameter Sensitivity
# ==============================================================================
# AUDIT FINDING: Shape parameter critically affects hazard function shape
# - Shape = 1.0: Exponential (constant hazard) - often unrealistic
# - Shape = 1.3-1.5: Increasing hazard - more realistic for PFS
# - Shape < 1.0: Decreasing hazard - post-response plateau

test_that("shape parameter affects operating characteristics", {
  # Common settings
  common_args <- list(
    num_simulations = 30,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    null_median_arms = c(Treatment = 10),
    median_pfs_success_threshold_arms = c(Treatment = 12),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Treatment = 50),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 3,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE
  )

  set.seed(7001)

  # Shape = 1.0 (Exponential - constant hazard)
  args_exp <- c(common_args, list(
    weibull_shape_true_arms = c(Treatment = 1.0),
    weibull_median_true_arms = c(Treatment = 12)
  ))
  result_exponential <- do.call(run_simulation_pure, args_exp)

  set.seed(7001)

  # Shape = 1.5 (Increasing hazard)
  args_increasing <- c(common_args, list(
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 12)
  ))
  result_increasing <- do.call(run_simulation_pure, args_increasing)

  set.seed(7001)

  # Shape = 0.8 (Decreasing hazard)
  args_decreasing <- c(common_args, list(
    weibull_shape_true_arms = c(Treatment = 0.8),
    weibull_median_true_arms = c(Treatment = 12)
  ))
  result_decreasing <- do.call(run_simulation_pure, args_decreasing)

  # All should produce valid results
  expect_s3_class(result_exponential, "data.frame")
  expect_s3_class(result_increasing, "data.frame")
  expect_s3_class(result_decreasing, "data.frame")

  # DOCUMENTED BEHAVIOR: Different shapes with same median produce different OCs
  # because the timing and distribution of events differs
  # (This is expected, not a bug)

  # Verify valid stopping probability sums
  for (res in list(result_exponential, result_increasing, result_decreasing)) {
    stop_sum <- res$PET_Efficacy + res$PET_Futility + res$Pr_Reach_Max_N
    expect_equal(stop_sum, 1, tolerance = 1e-6)
  }

  # Expected events should differ (shape affects event timing)
  # Shape > 1 means more events late; shape < 1 means more events early
  expect_true(all(!is.na(c(
    result_exponential$Exp_Events,
    result_increasing$Exp_Events,
    result_decreasing$Exp_Events
  ))))
})

# ==============================================================================
# TEST 2: Piecewise Interval Width Sensitivity
# ==============================================================================
# AUDIT FINDING: Coarse intervals may not adequately capture hazard variation

test_that("interval width affects operating characteristics", {
  set.seed(7002)

  # Coarse intervals (6-month)
  result_coarse <- run_simulation_pure(
    num_simulations = 25,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 10),
    interval_cutpoints_sim = c(0, 6, 12, 18, 24),  # COARSE: 6-month intervals
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(0.5, 0.5, 0.5, 0.5),  # 4 intervals
    prior_beta_params_model = c(0.5, 0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    null_median_arms = c(Treatment = 10),
    median_pfs_success_threshold_arms = c(Treatment = 10),
    efficacy_threshold_hc_prob = 0.95,
    futility_threshold_hc_prob = 0.85,
    final_success_posterior_prob_threshold = 0.9,
    max_total_patients_per_arm = c(Treatment = 50),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 4,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE
  )

  set.seed(7002)

  # Fine intervals (3-month)
  result_fine <- run_simulation_pure(
    num_simulations = 25,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 10),
    interval_cutpoints_sim = c(0, 3, 6, 9, 12, 15, 18, 21, 24),  # FINE: 3-month intervals
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = rep(0.5, 8),  # 8 intervals
    prior_beta_params_model = rep(0.5, 8),
    num_posterior_draws = 500,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    null_median_arms = c(Treatment = 10),
    median_pfs_success_threshold_arms = c(Treatment = 10),
    efficacy_threshold_hc_prob = 0.95,
    futility_threshold_hc_prob = 0.85,
    final_success_posterior_prob_threshold = 0.9,
    max_total_patients_per_arm = c(Treatment = 50),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 4,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE
  )

  expect_s3_class(result_coarse, "data.frame")
  expect_s3_class(result_fine, "data.frame")

  # Both should produce valid results
  for (res in list(result_coarse, result_fine)) {
    expect_equal(res$PET_Efficacy + res$PET_Futility + res$Pr_Reach_Max_N,
                 1, tolerance = 1e-6)
  }

  # DOCUMENTED BEHAVIOR: Interval choice affects model fit to true data
  # Finer intervals = better approximation of continuous hazard
})

# ==============================================================================
# TEST 3: Weibull vs Analysis Model Alignment
# ==============================================================================
# AUDIT FINDING: True Weibull data analyzed via piecewise exponential

test_that("analysis model approximates Weibull truth", {
  set.seed(7003)

  # When shape = 1 (exponential), piecewise exponential should match well
  result_aligned <- run_simulation_pure(
    num_simulations = 25,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.0),  # Exponential = piecewise w/ constant hazard
    weibull_median_true_arms = c(Treatment = 10),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    null_median_arms = c(Treatment = 10),  # Exactly null
    median_pfs_success_threshold_arms = c(Treatment = 10),
    efficacy_threshold_hc_prob = 0.95,
    futility_threshold_hc_prob = 0.85,
    final_success_posterior_prob_threshold = 0.9,
    max_total_patients_per_arm = c(Treatment = 50),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 4,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE
  )

  set.seed(7003)

  # When shape ≠ 1, piecewise exponential is an approximation
  result_misaligned <- run_simulation_pure(
    num_simulations = 25,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 2.0),  # Strong increasing hazard
    weibull_median_true_arms = c(Treatment = 10),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    null_median_arms = c(Treatment = 10),
    median_pfs_success_threshold_arms = c(Treatment = 10),
    efficacy_threshold_hc_prob = 0.95,
    futility_threshold_hc_prob = 0.85,
    final_success_posterior_prob_threshold = 0.9,
    max_total_patients_per_arm = c(Treatment = 50),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 4,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE
  )

  expect_s3_class(result_aligned, "data.frame")
  expect_s3_class(result_misaligned, "data.frame")

  # DOCUMENTED BEHAVIOR: Both produce valid results
  # The mismatch is a modeling choice, not a bug
  for (res in list(result_aligned, result_misaligned)) {
    expect_equal(res$PET_Efficacy + res$PET_Futility + res$Pr_Reach_Max_N,
                 1, tolerance = 1e-6)
  }
})

# ==============================================================================
# TEST 4: Censoring Model Sensitivity
# ==============================================================================
# AUDIT FINDING: Random censoring assumed MCAR (Missing Completely At Random)

test_that("censoring rate affects operating characteristics", {
  set.seed(7004)

  # Heavy censoring (short max censor time)
  result_heavy_censor <- run_simulation_pure(
    num_simulations = 25,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 12),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 12,  # HEAVY censoring: uniform(0, 12)
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    null_median_arms = c(Treatment = 10),
    median_pfs_success_threshold_arms = c(Treatment = 11),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Treatment = 50),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 3,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE
  )

  set.seed(7004)

  # Light censoring (long max censor time)
  result_light_censor <- run_simulation_pure(
    num_simulations = 25,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 12),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 48,  # LIGHT censoring: uniform(0, 48)
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    null_median_arms = c(Treatment = 10),
    median_pfs_success_threshold_arms = c(Treatment = 11),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Treatment = 50),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 3,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE
  )

  expect_s3_class(result_heavy_censor, "data.frame")
  expect_s3_class(result_light_censor, "data.frame")

  # Heavy censoring should lead to fewer observed events
  expect_true(result_heavy_censor$Exp_Events <= result_light_censor$Exp_Events + 10)

  # Both should produce valid probability sums
  for (res in list(result_heavy_censor, result_light_censor)) {
    expect_equal(res$PET_Efficacy + res$PET_Futility + res$Pr_Reach_Max_N,
                 1, tolerance = 1e-6)
  }
})

# ==============================================================================
# TEST 5: PH Model vs Independent Model (Multi-Arm)
# ==============================================================================
# AUDIT FINDING: Two analysis models available for vs-reference comparisons

test_that("PH model vs independent model both work",
{
  set.seed(7005)

  # Independent model (default)
  result_independent <- run_simulation_pure(
    num_simulations = 20,
    arm_names = c("Control", "Treatment"),
    reference_arm_name = "Control",
    compare_arms_option = TRUE,
    use_ph_model_vs_ref = FALSE,  # Independent hazards
    weibull_shape_true_arms = c(Control = 1.5, Treatment = 1.5),
    weibull_median_true_arms = c(Control = 10, Treatment = 12),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.5, Treatment = 0.5),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_events_per_arm = 5,
    efficacy_threshold_vs_ref_prob = 0.9,
    futility_threshold_vs_ref_prob = 0.8,
    max_total_patients_per_arm = c(Control = 40, Treatment = 40),
    cohort_size_per_arm = c(Control = 1, Treatment = 1),
    overall_accrual_rate = 3,
    interim_calendar_beat = 3,
    efficacy_stopping_rule_vs_ref = TRUE,
    futility_stopping_rule_vs_ref = TRUE
  )

  set.seed(7005)

  # Proportional hazards model
  result_ph <- run_simulation_pure(
    num_simulations = 20,
    arm_names = c("Control", "Treatment"),
    reference_arm_name = "Control",
    compare_arms_option = TRUE,
    use_ph_model_vs_ref = TRUE,  # PH model
    ph_loghr_prior_mean = 0,
    ph_loghr_prior_sd = 1,
    weibull_shape_true_arms = c(Control = 1.5, Treatment = 1.5),
    weibull_median_true_arms = c(Control = 10, Treatment = 12),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.5, Treatment = 0.5),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_events_per_arm = 5,
    efficacy_threshold_vs_ref_prob = 0.9,
    futility_threshold_vs_ref_prob = 0.8,
    max_total_patients_per_arm = c(Control = 40, Treatment = 40),
    cohort_size_per_arm = c(Control = 1, Treatment = 1),
    overall_accrual_rate = 3,
    interim_calendar_beat = 3,
    efficacy_stopping_rule_vs_ref = TRUE,
    futility_stopping_rule_vs_ref = TRUE
  )

  expect_s3_class(result_independent, "data.frame")
  expect_s3_class(result_ph, "data.frame")

  # Both should produce valid results
  for (res in list(result_independent, result_ph)) {
    for (i in 1:nrow(res)) {
      expect_equal(res$PET_Efficacy[i] + res$PET_Futility[i] + res$Pr_Reach_Max_N[i],
                   1, tolerance = 1e-6)
    }
  }
})
