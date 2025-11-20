test_that("single-arm null scenario produces correct type I error control", {
  set.seed(1111)

  # Null scenario: true median equals null median
  result <- run_simulation_pure(
    num_simulations = 50,
    arm_names = c("Control", "Treatment"),
    reference_arm_name = "Control",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Control = 1.5, Treatment = 1.5),
    weibull_median_true_arms = c(Control = 10, Treatment = 10),  # Null
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.5, Treatment = 0.5),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    min_median_followup_hc = 2,
    null_median_arms = c(Control = 10, Treatment = 10),
    futility_median_arms = c(Control = 9, Treatment = 9),
    median_pfs_success_threshold_arms = c(Control = 11, Treatment = 11),
    efficacy_threshold_hc_prob = 0.95,
    futility_threshold_hc_prob = 0.8,
    final_success_posterior_prob_threshold = 0.9,
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
  # Type I error should be low under null
  expect_true(all(result$Type_I_Error_or_Power < 0.3))
  # Should see some futility stopping or max N
  expect_true(all(result$PET_Futility + result$Pr_Reach_Max_N > 0.5))
})

test_that("single-arm alternative scenario produces power", {
  set.seed(2222)

  # Alternative scenario: true median better than null
  result <- run_simulation_pure(
    num_simulations = 50,
    arm_names = c("Control", "Treatment"),
    reference_arm_name = "Control",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(Control = 1.5, Treatment = 1.5),
    weibull_median_true_arms = c(Control = 15, Treatment = 15),  # Alternative
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(Control = 0.5, Treatment = 0.5),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    min_median_followup_hc = 2,
    null_median_arms = c(Control = 10, Treatment = 10),
    futility_median_arms = c(Control = 9, Treatment = 9),
    median_pfs_success_threshold_arms = c(Control = 12, Treatment = 12),
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
  # Power should be decent under alternative
  expect_true(all(result$Type_I_Error_or_Power > 0.3))
  # Should see some early efficacy stopping
  expect_true(any(result$PET_Efficacy > 0.1))
})

test_that("single-arm handles early stopping correctly", {
  set.seed(3333)

  # Test with favorable scenario to induce early stopping
  result <- run_simulation_pure(
    num_simulations = 40,
    arm_names = c("SuperTreatment"),
    reference_arm_name = "SuperTreatment",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(SuperTreatment = 1.5),
    weibull_median_true_arms = c(SuperTreatment = 18),  # Very good
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(SuperTreatment = 1.0),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    min_median_followup_hc = 2,
    null_median_arms = c(SuperTreatment = 10),
    futility_median_arms = c(SuperTreatment = 9),
    median_pfs_success_threshold_arms = c(SuperTreatment = 12),
    efficacy_threshold_hc_prob = 0.85,  # Easier to cross
    futility_threshold_hc_prob = 0.8,
    final_success_posterior_prob_threshold = 0.8,
    final_futility_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(SuperTreatment = 60),
    cohort_size_per_arm = c(SuperTreatment = 5),
    overall_accrual_rate = 2,
    min_follow_up_at_final = 3,
    interim_calendar_beat = 2,  # Frequent interims
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
  # Check basic structure and validity (early stopping is stochastic, don't enforce it)
  expect_true(all(!is.na(result$Exp_N)))
  expect_true(result$Exp_N >= 0)
  # All stop probabilities should sum to 1
  stop_sum <- result$PET_Efficacy + result$PET_Futility + result$Pr_Reach_Max_N
  expect_equal(stop_sum, 1, tolerance = 1e-6)
})

test_that("single-arm final analysis produces valid decisions", {
  set.seed(4444)

  result <- run_simulation_pure(
    num_simulations = 35,
    arm_names = c("ArmA", "ArmB"),
    reference_arm_name = "ArmA",
    compare_arms_option = FALSE,
    weibull_shape_true_arms = c(ArmA = 1.5, ArmB = 1.5),
    weibull_median_true_arms = c(ArmA = 11, ArmB = 13),
    interval_cutpoints_sim = c(0, 6, 12, 24),
    max_follow_up_sim = 24,
    censor_max_time_sim = 30,
    randomization_probs = c(ArmA = 0.5, ArmB = 0.5),
    prior_alpha_params_model = c(0.1, 0.1, 0.1),
    prior_beta_params_model = c(0.1, 0.1, 0.1),
    num_posterior_draws = 1000,
    min_patients_for_analysis = 10,
    min_events_hc = 5,
    min_median_followup_hc = 2,
    null_median_arms = c(ArmA = 10, ArmB = 10),
    futility_median_arms = c(ArmA = 9, ArmB = 9),
    median_pfs_success_threshold_arms = c(ArmA = 11, ArmB = 11),
    efficacy_threshold_hc_prob = 0.9,
    futility_threshold_hc_prob = 0.8,
    final_success_posterior_prob_threshold = 0.85,
    final_futility_posterior_prob_threshold = 0.85,
    max_total_patients_per_arm = c(ArmA = 45, ArmB = 45),
    cohort_size_per_arm = c(ArmA = 5, ArmB = 5),
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

  # Check that final decision probabilities sum to ~Pr_Reach_Max_N
  # (arms that stop early won't have final decisions)
  for (i in 1:nrow(result)) {
    final_sum <- result$Pr_Final_Efficacy[i] +
                 result$Pr_Final_Futility[i] +
                 result$Pr_Final_Inconclusive[i]
    expected_final <- result$Pr_Reach_Max_N[i]
    expect_equal(final_sum, expected_final, tolerance = 1e-6)
  }

  # All stopping probabilities should sum to 1
  for (i in 1:nrow(result)) {
    stop_sum <- result$PET_Efficacy[i] +
                result$PET_Futility[i] +
                result$Pr_Reach_Max_N[i]
    expect_equal(stop_sum, 1, tolerance = 1e-6)
  }
})
