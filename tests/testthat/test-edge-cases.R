# Test file: Edge Cases
# Purpose: Test edge case handling identified during code audit

# ==============================================================================
# TEST 1: Zero Events in an Interval
# ==============================================================================
# AUDIT FINDING: Gamma posteriors with weak priors and zero events can be unstable

test_that("handles zero events in an interval gracefully", {
  set.seed(6001)

  # Very short follow-up with low event rate = likely zero events in later intervals
  result <- run_simulation_pure(
    num_simulations = 20,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 36),  # Very long median = few events
    interval_cutpoints_sim = c(0, 3, 6, 12, 24),  # Fine intervals
    max_follow_up_sim = 12,  # Short follow-up
    censor_max_time_sim = 15,
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(0.1, 0.1, 0.1, 0.1),  # Weak priors
    prior_beta_params_model = c(0.1, 0.1, 0.1, 0.1),
    num_posterior_draws = 500,
    min_patients_for_analysis = 5,
    min_events_hc = 1,  # Very low gate
    null_median_arms = c(Treatment = 20),
    median_pfs_success_threshold_arms = c(Treatment = 25),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Treatment = 30),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 3,
    interim_calendar_beat = 2,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)

  # Should produce valid OCs even with sparse data
  expect_true(all(!is.na(result$Type_I_Error_or_Power)))
  expect_true(all(!is.na(result$Exp_N)))
  expect_true(result$Exp_N > 0 && result$Exp_N <= 30)
})

# ==============================================================================
# TEST 2: Median Not Reached Within Follow-up
# ==============================================================================
# AUDIT FINDING: Tail extrapolation uses last interval hazard when median
# exceeds the piecewise grid.

test_that("handles median not reached within follow-up", {
  set.seed(6002)

  # Very long true median relative to follow-up
  result <- run_simulation_pure(
    num_simulations = 20,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 48),  # Much longer than 24 month FU
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_patients_for_analysis = 10,
    min_events_hc = 3,
    null_median_arms = c(Treatment = 30),
    median_pfs_success_threshold_arms = c(Treatment = 40),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Treatment = 40),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 4,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE
  )

  expect_s3_class(result, "data.frame")

  # Simulation should complete and produce valid results
  expect_true(all(!is.na(result$Type_I_Error_or_Power)))
  expect_true(all(result$Type_I_Error_or_Power >= 0 & result$Type_I_Error_or_Power <= 1))

  # Median extrapolation documented: uses last interval hazard
  # This is tested implicitly - no errors should occur
})

# ==============================================================================
# TEST 3: Very Extreme Gate Values
# ==============================================================================
# AUDIT FINDING: Extreme gate values can prevent any interim analyses

test_that("handles extreme gate values", {
  set.seed(6003)

  # Gates that can never be satisfied
  result_impossible <- run_simulation_pure(
    num_simulations = 15,
    arm_names = c("Control", "Treatment"),
    reference_arm_name = "Control",
    compare_arms_option = TRUE,
    weibull_shape_true_arms = c(Control = 1.5, Treatment = 1.5),
    weibull_median_true_arms = c(Control = 10, Treatment = 15),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.5, Treatment = 0.5),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_events_per_arm = 100,  # IMPOSSIBLE: more events than patients
    min_median_followup_per_arm = 0,
    efficacy_threshold_vs_ref_prob = 0.9,
    futility_threshold_vs_ref_prob = 0.8,
    max_total_patients_per_arm = c(Control = 30, Treatment = 30),
    cohort_size_per_arm = c(Control = 1, Treatment = 1),
    overall_accrual_rate = 3,
    interim_calendar_beat = 2,
    efficacy_stopping_rule_vs_ref = TRUE,
    futility_stopping_rule_vs_ref = TRUE
  )

  expect_s3_class(result_impossible, "data.frame")

  # With impossible gates, no early stopping should occur
  trt_row <- result_impossible[result_impossible$Arm_Name == "Treatment", ]
  expect_equal(trt_row$PET_Efficacy, 0, tolerance = 1e-6)
  expect_equal(trt_row$PET_Futility, 0, tolerance = 1e-6)
  expect_equal(trt_row$Pr_Reach_Max_N, 1, tolerance = 1e-6)
})

test_that("handles zero gates (always-open)", {
  set.seed(6004)

  # All gates set to zero - should allow earliest possible interims
  result_open <- run_simulation_pure(
    num_simulations = 20,
    arm_names = c("Control", "Treatment"),
    reference_arm_name = "Control",
    compare_arms_option = TRUE,
    weibull_shape_true_arms = c(Control = 1.5, Treatment = 1.5),
    weibull_median_true_arms = c(Control = 10, Treatment = 10),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.5, Treatment = 0.5),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_events_per_arm = 0,  # Open gate
    min_median_followup_per_arm = 0,  # Open gate
    min_person_time_frac_per_arm = 0,  # Open gate
    efficacy_threshold_vs_ref_prob = 0.9,
    futility_threshold_vs_ref_prob = 0.8,
    max_total_patients_per_arm = c(Control = 40, Treatment = 40),
    cohort_size_per_arm = c(Control = 1, Treatment = 1),
    overall_accrual_rate = 3,
    interim_calendar_beat = 2,
    efficacy_stopping_rule_vs_ref = TRUE,
    futility_stopping_rule_vs_ref = TRUE
  )

  expect_s3_class(result_open, "data.frame")

  # With open gates, early stopping is possible
  # (actual rates depend on the random data)
  for (i in 1:nrow(result_open)) {
    stop_sum <- result_open$PET_Efficacy[i] +
      result_open$PET_Futility[i] +
      result_open$Pr_Reach_Max_N[i]
    expect_equal(stop_sum, 1, tolerance = 1e-6)
  }
})

# ==============================================================================
# TEST 4: Single Patient / Minimal Sample Size
# ==============================================================================
# AUDIT FINDING: Very small samples can lead to unstable posteriors

test_that("handles minimal sample size", {
  set.seed(6005)

  result <- run_simulation_pure(
    num_simulations = 15,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 12),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_patients_for_analysis = 1,  # Allow analysis with 1 patient
    min_events_hc = 1,
    null_median_arms = c(Treatment = 10),
    median_pfs_success_threshold_arms = c(Treatment = 11),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Treatment = 5),  # Very small max N
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 1,
    interim_calendar_beat = 2,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE
  )

  expect_s3_class(result, "data.frame")
  expect_true(result$Exp_N <= 5)  # Should not exceed max
  expect_true(result$Exp_N >= 0)  # Should be non-negative
})

# ==============================================================================
# TEST 5: Unbalanced Randomization
# ==============================================================================
# AUDIT FINDING: Gate scaling is proportional to randomization probabilities

test_that("handles unbalanced randomization with gate scaling", {
  set.seed(6006)

  # Very unbalanced randomization
  result <- run_simulation_pure(
    num_simulations = 20,
    arm_names = c("Control", "Treatment"),
    reference_arm_name = "Control",
    compare_arms_option = TRUE,
    weibull_shape_true_arms = c(Control = 1.5, Treatment = 1.5),
    weibull_median_true_arms = c(Control = 10, Treatment = 10),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.2, Treatment = 0.8),  # 1:4 ratio
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_events_per_arm = 10,  # Will be scaled proportionally
    min_median_followup_per_arm = 3,
    min_person_time_frac_per_arm = 0.15,  # Will be scaled proportionally
    efficacy_threshold_vs_ref_prob = 0.9,
    futility_threshold_vs_ref_prob = 0.8,
    max_total_patients_per_arm = c(Control = 25, Treatment = 100),
    cohort_size_per_arm = c(Control = 1, Treatment = 1),
    overall_accrual_rate = 4,
    interim_calendar_beat = 3,
    efficacy_stopping_rule_vs_ref = TRUE,
    futility_stopping_rule_vs_ref = TRUE
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)

  # Check valid results for both arms
  for (i in 1:nrow(result)) {
    expect_true(!is.na(result$Exp_N[i]))
    expect_true(result$Exp_N[i] >= 0)
    stop_sum <- result$PET_Efficacy[i] + result$PET_Futility[i] + result$Pr_Reach_Max_N[i]
    expect_equal(stop_sum, 1, tolerance = 1e-6)
  }
})

# ==============================================================================
# TEST 6: Very High Posterior Thresholds
# ==============================================================================
# AUDIT FINDING: Threshold selection critically affects operating characteristics

test_that("handles very high posterior thresholds", {
  set.seed(6007)

  # Extremely conservative thresholds
  result <- run_simulation_pure(
    num_simulations = 20,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 15),  # Strong effect
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
    median_pfs_success_threshold_arms = c(Treatment = 11),
    efficacy_threshold_hc_prob = 0.999,  # VERY conservative
    futility_threshold_hc_prob = 0.999,  # VERY conservative
    final_success_posterior_prob_threshold = 0.999,
    max_total_patients_per_arm = c(Treatment = 50),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 3,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE
  )

  expect_s3_class(result, "data.frame")

  # With very conservative thresholds, early stopping should be rare
  # Most trials should reach max N
  expect_true(result$Pr_Reach_Max_N >= 0.5 ||
                (result$PET_Efficacy + result$PET_Futility) < 0.5)
})

# ==============================================================================
# TEST 7: Person-Time Milestone Scheduling
# ==============================================================================
# AUDIT FINDING: Person-time milestones are an alternative to calendar beats

test_that("person-time milestones work correctly", {
  set.seed(6008)

  result <- run_simulation_pure(
    num_simulations = 15,
    arm_names = c("Control", "Treatment"),
    reference_arm_name = "Control",
    compare_arms_option = TRUE,
    weibull_shape_true_arms = c(Control = 1.5, Treatment = 1.5),
    weibull_median_true_arms = c(Control = 10, Treatment = 10),
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
    person_time_milestones = c(0.25, 0.5, 0.75, 1.0),  # PT-based scheduling
    latest_calendar_look = 48,
    efficacy_stopping_rule_vs_ref = TRUE,
    futility_stopping_rule_vs_ref = TRUE
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)

  # Verify valid results
  for (i in 1:nrow(result)) {
    stop_sum <- result$PET_Efficacy[i] + result$PET_Futility[i] + result$Pr_Reach_Max_N[i]
    expect_equal(stop_sum, 1, tolerance = 1e-6)
  }
})
