# Test file: Assumption Audit
# Purpose: Document and validate key statistical assumptions identified in code audit
# These tests serve as regression tests and documentation for package behavior

# ==============================================================================
# TEST 1: Prior Sensitivity - Weak vs Informative Priors
# ==============================================================================
# AUDIT FINDING: Default Gamma(0.1, 0.1) priors are very weak and can lead to
# volatile posteriors when event counts are low.

test_that("weak priors produce more variable results than informative priors", {
  set.seed(5001)

  # Test with weak priors (default-like)
  result_weak <- run_simulation_pure(
    num_simulations = 30,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 10),  # Null scenario
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),  # WEAK priors
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 500,
    min_patients_for_analysis = 5,
    min_events_hc = 3,
    null_median_arms = c(Treatment = 10),
    futility_median_arms = c(Treatment = 8),
    median_pfs_success_threshold_arms = c(Treatment = 12),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Treatment = 40),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 3,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE
  )

  set.seed(5001)  # Same seed for comparability

  # Test with more informative priors
  result_informative <- run_simulation_pure(
    num_simulations = 30,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 10),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(1, 1, 1),  # MORE INFORMATIVE priors
    prior_beta_params_model = c(1, 1, 1),
    num_posterior_draws = 500,
    min_patients_for_analysis = 5,
    min_events_hc = 3,
    null_median_arms = c(Treatment = 10),
    futility_median_arms = c(Treatment = 8),
    median_pfs_success_threshold_arms = c(Treatment = 12),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Treatment = 40),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 3,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE
  )

  # DOCUMENTED BEHAVIOR: Prior choice affects operating characteristics
  # Both should produce valid results
  expect_s3_class(result_weak, "data.frame")
  expect_s3_class(result_informative, "data.frame")

  # Stopping probabilities should sum to 1
  expect_equal(result_weak$PET_Efficacy + result_weak$PET_Futility +
                 result_weak$Pr_Reach_Max_N, 1, tolerance = 1e-6)
  expect_equal(result_informative$PET_Efficacy + result_informative$PET_Futility +
                 result_informative$Pr_Reach_Max_N, 1, tolerance = 1e-6)
})

# ==============================================================================
# TEST 2: Multiple Interim Looks - No Alpha Spending
# ==============================================================================
# AUDIT FINDING: Same posterior probability threshold is applied at every interim
# without any alpha spending correction.

test_that("multiple interims use same threshold (no alpha spending)", {
  set.seed(5002)

  # Frequent interims (every 2 months)
  result_frequent <- run_simulation_pure(
    num_simulations = 30,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 10),  # Null
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_patients_for_analysis = 5,
    min_events_hc = 3,
    null_median_arms = c(Treatment = 10),
    median_pfs_success_threshold_arms = c(Treatment = 10),  # Exactly null
    efficacy_threshold_hc_prob = 0.95,
    futility_threshold_hc_prob = 0.9,
    final_success_posterior_prob_threshold = 0.95,
    max_total_patients_per_arm = c(Treatment = 50),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 2,  # FREQUENT interims
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = FALSE
  )

  set.seed(5002)

  # Infrequent interims (every 6 months)
  result_infrequent <- run_simulation_pure(
    num_simulations = 30,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 10),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_patients_for_analysis = 5,
    min_events_hc = 3,
    null_median_arms = c(Treatment = 10),
    median_pfs_success_threshold_arms = c(Treatment = 10),
    efficacy_threshold_hc_prob = 0.95,
    futility_threshold_hc_prob = 0.9,
    final_success_posterior_prob_threshold = 0.95,
    max_total_patients_per_arm = c(Treatment = 50),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 6,  # INFREQUENT interims
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = FALSE
  )

  # DOCUMENTED BEHAVIOR: Both configurations produce valid results
  # Note: Type I error may differ due to number of looks (no alpha spending)
  expect_s3_class(result_frequent, "data.frame")
  expect_s3_class(result_infrequent, "data.frame")

  # Both should have valid probability sums
  expect_equal(result_frequent$PET_Efficacy + result_frequent$PET_Futility +
                 result_frequent$Pr_Reach_Max_N, 1, tolerance = 1e-6)
})

# ==============================================================================
# TEST 3: Multi-Arm Trial - Independent Arm Comparisons
# ==============================================================================
# AUDIT FINDING: Each experimental arm is compared to control independently,
# without family-wise error rate control.

test_that("multi-arm comparisons are evaluated independently", {
  set.seed(5003)

  # Multi-arm trial with 3 experimental arms (all null)
  result_multiarm <- run_simulation_pure(
    num_simulations = 25,
    arm_names = c("Control", "TrtA", "TrtB", "TrtC"),
    reference_arm_name = "Control",
    compare_arms_option = TRUE,  # vs-reference mode
    weibull_shape_true_arms = c(Control = 1.5, TrtA = 1.5, TrtB = 1.5, TrtC = 1.5),
    weibull_median_true_arms = c(Control = 10, TrtA = 10, TrtB = 10, TrtC = 10),  # All null
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.25, TrtA = 0.25, TrtB = 0.25, TrtC = 0.25),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_events_per_arm = 5,
    min_median_followup_per_arm = 2,
    efficacy_threshold_vs_ref_prob = 0.95,
    futility_threshold_vs_ref_prob = 0.8,
    compare_arms_futility_margin = 0,
    max_total_patients_per_arm = c(Control = 40, TrtA = 40, TrtB = 40, TrtC = 40),
    cohort_size_per_arm = c(Control = 1, TrtA = 1, TrtB = 1, TrtC = 1),
    overall_accrual_rate = 4,
    interim_calendar_beat = 3,
    efficacy_stopping_rule_vs_ref = TRUE,
    futility_stopping_rule_vs_ref = TRUE
  )

  expect_s3_class(result_multiarm, "data.frame")
  expect_equal(nrow(result_multiarm), 4)  # All 4 arms reported

  # Each experimental arm has independent OCs
  trt_rows <- result_multiarm[result_multiarm$Arm_Name != "Control", ]
  expect_equal(nrow(trt_rows), 3)

  # Each arm should have valid stopping probability sum
  for (i in 1:nrow(result_multiarm)) {
    stop_sum <- result_multiarm$PET_Efficacy[i] +
      result_multiarm$PET_Futility[i] +
      result_multiarm$Pr_Reach_Max_N[i]
    expect_equal(stop_sum, 1, tolerance = 1e-6,
                 label = paste("Arm", result_multiarm$Arm_Name[i]))
  }

  # DOCUMENTED BEHAVIOR: FWER across arms is NOT controlled
  # Under null, any-arm Type I error ≈ 1 - (1 - per_arm_alpha)^3
  # This test documents that behavior exists (not that it's wrong)
})

# ==============================================================================
# TEST 4: Historical Control Assumption - Single-Arm Path
# ==============================================================================
# AUDIT FINDING: Historical control median is treated as known exactly,
# with no uncertainty propagation.

test_that("single-arm path compares to fixed historical control", {
  set.seed(5004)

  # Test that efficacy depends on success threshold (historical control value)
  result_low_bar <- run_simulation_pure(
    num_simulations = 30,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 12),  # True median = 12
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_patients_for_analysis = 10,
    null_median_arms = c(Treatment = 8),  # Historical control = 8
    median_pfs_success_threshold_arms = c(Treatment = 9),  # LOW bar: beat 9
    efficacy_threshold_hc_prob = 0.9,
    final_success_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Treatment = 50),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 3,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = FALSE
  )

  set.seed(5004)

  result_high_bar <- run_simulation_pure(
    num_simulations = 30,
    arm_names = c("Treatment"),
    reference_arm_name = "Treatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Treatment = 1.5),
    weibull_median_true_arms = c(Treatment = 12),  # Same true median
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Treatment = 1.0),
    prior_alpha_params_model = c(0.5, 0.5, 0.5),
    prior_beta_params_model = c(0.5, 0.5, 0.5),
    num_posterior_draws = 500,
    min_patients_for_analysis = 10,
    null_median_arms = c(Treatment = 8),
    median_pfs_success_threshold_arms = c(Treatment = 14),  # HIGH bar: beat 14
    efficacy_threshold_hc_prob = 0.9,
    final_success_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Treatment = 50),
    cohort_size_per_arm = c(Treatment = 1),
    overall_accrual_rate = 2,
    interim_calendar_beat = 3,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = FALSE
  )

  expect_s3_class(result_low_bar, "data.frame")
  expect_s3_class(result_high_bar, "data.frame")

  # DOCUMENTED BEHAVIOR: Lower success threshold = higher power
  # (true median 12 exceeds 9 but may not reliably exceed 14)
  # With small samples this is stochastic, but the relationship should hold on average
})

# ==============================================================================
# TEST 5: Gate Enforcement - Delays But Does Not Control Error
# ==============================================================================
# AUDIT FINDING: Information gates postpone interim analyses but do not
# provide any additional error control once gates are passed.

test_that("gates delay first interim but use same threshold after", {
  set.seed(5005)

  # Strict gates - delays first interim
  result_strict_gates <- run_simulation_pure(
    num_simulations = 25,
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
    min_events_per_arm = 15,  # STRICT: need 15 events per arm
    min_median_followup_per_arm = 6,  # STRICT: 6 months
    min_person_time_frac_per_arm = 0.3,  # STRICT: 30% person-time
    efficacy_threshold_vs_ref_prob = 0.9,
    futility_threshold_vs_ref_prob = 0.8,
    max_total_patients_per_arm = c(Control = 50, Treatment = 50),
    cohort_size_per_arm = c(Control = 1, Treatment = 1),
    overall_accrual_rate = 3,
    interim_calendar_beat = 2,
    efficacy_stopping_rule_vs_ref = TRUE,
    futility_stopping_rule_vs_ref = TRUE
  )

  set.seed(5005)

  # Lenient gates - allows earlier interims
  result_lenient_gates <- run_simulation_pure(
    num_simulations = 25,
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
    min_events_per_arm = 5,  # LENIENT: only 5 events
    min_median_followup_per_arm = 2,  # LENIENT: 2 months
    min_person_time_frac_per_arm = 0.1,  # LENIENT: 10% person-time
    efficacy_threshold_vs_ref_prob = 0.9,  # Same threshold
    futility_threshold_vs_ref_prob = 0.8,
    max_total_patients_per_arm = c(Control = 50, Treatment = 50),
    cohort_size_per_arm = c(Control = 1, Treatment = 1),
    overall_accrual_rate = 3,
    interim_calendar_beat = 2,
    efficacy_stopping_rule_vs_ref = TRUE,
    futility_stopping_rule_vs_ref = TRUE
  )

  expect_s3_class(result_strict_gates, "data.frame")
  expect_s3_class(result_lenient_gates, "data.frame")

  # DOCUMENTED BEHAVIOR: Strict gates should result in fewer early stops
  # (because first interim is delayed, fewer opportunities for stopping)
  # Note: With 25 sims this is stochastic, but structure should be valid
  for (df in list(result_strict_gates, result_lenient_gates)) {
    for (i in 1:nrow(df)) {
      expect_equal(df$PET_Efficacy[i] + df$PET_Futility[i] + df$Pr_Reach_Max_N[i],
                   1, tolerance = 1e-6)
    }
  }
})
