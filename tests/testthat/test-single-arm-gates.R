test_that("single-arm respects person-time fraction gates with proportional scaling", {
  set.seed(5555)

  # Test with unbalanced randomization - expect proportional scaling
  result <- run_simulation_pure(
    num_simulations = 30,
    arm_names = c("Control", "Treatment"),
    reference_arm_name = "Control",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Control = 1.5, Treatment = 1.5),
    weibull_median_true_arms = c(Control = 10, Treatment = 10),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.3, Treatment = 0.7),  # Unbalanced!
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    min_median_followup_hc = 2,
    min_person_time_frac_per_arm = 0.15,  # Should be scaled by randomization ratio
    null_median_arms = c(Control = 8, Treatment = 8),
    futility_median_arms = c(Control = 9, Treatment = 9),
    median_pfs_success_threshold_arms = c(Control = 11, Treatment = 11),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    final_futility_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Control = 50, Treatment = 50),
    cohort_size_per_arm = c(Control = 5, Treatment = 5),
    overall_accrual_rate = 2,
    min_follow_up_at_final = 3,
    interim_calendar_beat = 3,
    diagnostics = FALSE,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE,
    efficacy_stopping_rule_vs_ref = FALSE,
    futility_stopping_rule_vs_ref = FALSE,
    pred_success_pp_threshold_hc = 1,
    pred_futility_pp_threshold_hc = 0,
    num_posterior_draws_pred = 100
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(result$Exp_N > 0 & result$Exp_N <= 50))
})

test_that("single-arm works with rebalance_after_events", {
  set.seed(6666)

  result <- run_simulation_pure(
    num_simulations = 20,
    arm_names = c("Arm_A", "Arm_B"),
    reference_arm_name = "Arm_A",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Arm_A = 1.5, Arm_B = 1.5),
    weibull_median_true_arms = c(Arm_A = 12, Arm_B = 14),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Arm_A = 0.5, Arm_B = 0.5),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    min_median_followup_hc = 2,
    null_median_arms = c(Arm_A = 8, Arm_B = 8),
    futility_median_arms = c(Arm_A = 9, Arm_B = 9),
    median_pfs_success_threshold_arms = c(Arm_A = 11, Arm_B = 11),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    final_futility_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(Arm_A = 40, Arm_B = 40),
    cohort_size_per_arm = c(Arm_A = 5, Arm_B = 5),
    overall_accrual_rate = 2,
    min_follow_up_at_final = 3,
    interim_calendar_beat = 3,
    rebalance_after_events = 15,  # Rebalance intervals after 15 events
    diagnostics = FALSE,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE,
    efficacy_stopping_rule_vs_ref = FALSE,
    futility_stopping_rule_vs_ref = FALSE,
    pred_success_pp_threshold_hc = 1,
    pred_futility_pp_threshold_hc = 0,
    num_posterior_draws_pred = 100
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  # Check that rebalancing didn't break anything
  expect_true(all(!is.na(result$Exp_N)))
  expect_true(all(result$Exp_N > 0))
})

test_that("single-arm with multiple gates enforced correctly", {
  set.seed(7777)

  # Use strict gates to test they're all enforced
  result <- run_simulation_pure(
    num_simulations = 25,
    arm_names = c("Experimental"),
    reference_arm_name = "Experimental",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Experimental = 1.5),
    weibull_median_true_arms = c(Experimental = 12),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Experimental = 1.0),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    min_patients_for_analysis = 15,  # Higher gates
    min_events_hc = 10,
    min_median_followup_hc = 4,
    min_person_time_frac_per_arm = 0.25,
    null_median_arms = c(Experimental = 8),
    futility_median_arms = c(Experimental = 9),
    median_pfs_success_threshold_arms = c(Experimental = 11),
    efficacy_threshold_hc_prob = 0.95,
    futility_threshold_hc_prob = 0.85,
    final_success_posterior_prob_threshold = 0.9,
    final_futility_posterior_prob_threshold = 0.9,
    max_total_patients_per_arm = c(Experimental = 60),
    cohort_size_per_arm = c(Experimental = 5),
    overall_accrual_rate = 2,
    min_follow_up_at_final = 3,
    interim_calendar_beat = 2,
    diagnostics = FALSE,
    efficacy_stopping_rule_hc = TRUE,
    futility_stopping_rule_hc = TRUE,
    efficacy_stopping_rule_vs_ref = FALSE,
    futility_stopping_rule_vs_ref = FALSE,
    pred_success_pp_threshold_hc = 1,
    pred_futility_pp_threshold_hc = 0,
    num_posterior_draws_pred = 100
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$Arm_Name, "Experimental")
  # With high gates, should generally reach higher sample sizes
  expect_true(result$Pr_Reach_Max_N > 0.3)
})
